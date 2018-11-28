#defines the size of each simulation chunk
const ASPU_NCHUNK = 10^4

#returns SPU(gamma) for all gammas. fast
function getspu!(spu, pows, z, n)
  @inbounds for i in eachindex(pows)
    @fastmath pows[i] != 0 && (spu[i] = z[1]^(pows[i]))
    pows[i] == 0 && (spu[i] = abs(z[1]))
  end
  @inbounds for j = 2:n
    for i in eachindex(pows)
      @fastmath (pows[i] != 0) && (spu[i] += z[j]^(pows[i]))
      (pows[i] == 0) && (abs(z[j]) > spu[i]) && (spu[i] = abs(z[j]))
    end
  end
  @inbounds for i = eachindex(pows)
    spu[i] = abs(spu[i])
  end
end

function getspu(pows::Array{Int64, 1}, z::Vector{T}, n::Int64) where {T<:Real}
  tmpspu = Array{T}(undef, length(pows))
  getspu!(tmpspu,pows,z,n)
  tmpspu
end

#get highest SPU(gamma) rank across gamma
function rank_spus!(rnk::AbstractArray{Int64, 2}, zb::Array{T,2}, B = size(zb, 2)) where {T<:Real}
  rnk1_v = view(rnk,1,1:B)
  for i in 1:B
    rnk1_v[i] = i
  end
  quicksort!(zb[1,1:B], rnk1_v)
  for (i, val) in enumerate(rnk1_v)
    rnk[2,val] = i
    rnk1_v[i] = i
  end
  for i in 2:(size(zb,1)-1)
    quicksort!(zb[i,1:B], rnk1_v)
    for (j, val) in enumerate(rnk1_v)
      rnk[2,val] < j && (rnk[2,val] = j)
      rnk1_v[j] = j
    end
  end
  quicksort!(zb[size(rnk,1),1:B], rnk1_v)
  for (j, val) in enumerate(rnk1_v)
    rnk[2,val] < j && (rnk[2,val] = j)
  end
  0
end


#add values exceeding threshold to array
function create_arref(mvn, thresh, pows)
    nchunk = ASPU_NCHUNK

    tmp = zeros(Float64, length(pows))
    np = length(pows)
    ntraits = length(mvn)
    out = zeros(Float64, length(pows), nchunk)
    n = 0
    narr = zeros(Int, length(pows))
    ran = rand(mvn, nchunk)
    for i in 1:nchunk
        getspu!(tmp, pows, ran[:,i], ntraits)
        topind = tmp .> thresh
        for p in eachindex(pows)[topind]
            narr[p] += 1
        end
        if sum(topind) > 0
            n += 1
            out[:, n] = tmp
        end
    end
    n, narr, out[:, 1:n]
end

#Start workers for initialization
function do_initwork(vars, jobs, results)
    pows, mvn = take!(vars)
    while true
        thresh = take!(jobs)
        out = create_arref(mvn, thresh, pows)
        put!(results, out)
    end
end

#function simulating SPUs in parallel
function init_aspu_par(pows::Vector{Int64}, mvn, maxiter::Int64; verbose = true)
    nchunk = ASPU_NCHUNK

    maxchunks = Int(maxiter/nchunk)
    ntraits = length(mvn)
    ranspu = zeros(length(pows), nchunk)
    tmp = zeros(length(pows))
    maxin = zeros(Int, ceil(Int, 1+log10(maxiter/nchunk)))
    maxin_arr = zeros(Int, length(pows), ceil(Int, 1+log10(maxiter/nchunk)))

    #for parallel; define channels and start workers
    thresh::Vector{Float64} = fill(0.0, length(pows))
    jobs = RemoteChannel(()->Channel{typeof(thresh)}(maxchunks))
    results = RemoteChannel(()->Channel{Any}(maxchunks))
    init_vars = RemoteChannel(()->Channel{Any}(length(workers())))

    for p in workers()
        put!(init_vars, (pows, mvn))
        remote_do(do_initwork, p, init_vars, jobs, results)
    end

    #big data structure to hold simulations
    allvals = [zeros(length(pows), nchunk*10) for i in 1:ceil(Int, 1+log10(maxiter/nchunk))]

    #fill first bucket to the brim
    verbose && println_timestamp("Generating SPUs for 1e-$(Int(log10(nchunk)))")
    ran = rand(mvn, nchunk)
    for i in 1:nchunk
        getspu!(view(allvals[1],:,i), pows, ran[:,i], ntraits)
    end
    maxin[1] = nchunk
    fill!(view(maxin_arr,:,1), nchunk)

    #pass each bucket of values to workers
    for i in eachindex(allvals)[2:end]
        verbose && println_timestamp("Generating SPUs for 1e-$(Int(log10(nchunk))+i-1)")
        iternow = nchunk*10^(i-1)
        thresh = [partialsort(view(allvals[i-1],j,:), Int(nchunk*0.10); rev = true) for j in eachindex(pows)]
        fullchunks = floor(Int, iternow/nchunk)

        for chunk in 1:fullchunks
            put!(jobs, thresh)
        end
        for chunk in 1:fullchunks
            tn, tnarr, tout = take!(results)
            maxin_arr[:,i] = maxin_arr[:,i] .+ tnarr
            allvals[i][:, (maxin[i]+1) : (maxin[i]+tn)] = tout
            maxin[i] += tn
        end
    end

    allranks = [ zeros(Int, 2, n) for n in maxin ]

    [ rank_spus!(allranks[i], allvals[i], maxin[i]) for i in eachindex(allvals) ]
    allsorted = [ [ sort(allvals[i][g,:], rev=true)[1:maxin_arr[g,i]] for g in eachindex(pows) ] for i in eachindex(allvals) ]

    verbose && println_timestamp("Simulations initialized")
    allsorted, allranks, maxin_arr, maxin
end

#function to process a single z-scores vector
function getaspu(z, allsorted, allranks, maxin_arr, maxin, pows)
    spu = getspu(pows, z, length(z))
    pval = zeros(Int, size(allsorted[1],1))
    ind_p = fill(length(allsorted), length(pows))
    i = 0
    B = maxin[1]
    while minimum(pval) < 850 && i < length(allsorted)
        i += 1
        for k in eachindex(pows)
            pval[k] = sum( spu[k] .< allsorted[i][k][1:maxin_arr[k, i]] )
            (ind_p[k] > i && pval[k] > 850) && (ind_p[k] = i)
        end
    end
    minp, gamma = findmin(pval)
    aspu_n = count(x->(x > maxin[i] - minp), allranks[i][2,:])
    aspu_p = (aspu_n + 1) / (B*10^(i-1) + 1)
    p_out = (pval .+ 1) ./ (B .* 10 .^ (ind_p .- 1) .+ 1)

    aspu_p, p_out, gamma
end

#Arguments are passed once through the channel, then workers are started
function do_aspuwork(vars, jobs, results)
    aspu_obj, pows, delim, trans, na = take!(vars)
    allsorted, allranks, maxin_arr, maxin = aspu_obj
    while true
        line = take!(jobs)
        ls = split(line, delim)
        if in(na, ls)
            put!(results, (ls, fill(na, length(pows)+2)))
        else
            z = trans*parse.(Float64, ls[2:end])
            out = getaspu(z, allsorted, allranks, maxin_arr, maxin, pows)
            put!(results, (ls, out))
        end
    end
end


invsd(mat::Matrix) = sqrt(inv(Diagonal(mat)))
cov2cor(mat::Matrix) = Matrix(Hermitian( invsd(mat) * mat * invsd(mat )))

function aspu(
    filein::AbstractString, niter::Int64;
    pows::Vector{Int64} = collect(0:8), invR_trans::Bool = false,
    covfile::AbstractString = "", plim::Float64 = 1e-4,
    delim::Char = '\t', noheader::Bool = false, skip::Int64 = 1, na::String = "NA",
    out::AbstractString = "", verbose::Bool = true, nosavecov::Bool = false,
    outtest::Real = Inf
    )

    nchunk = ASPU_NCHUNK

    maxiter = Int(10^(ceil(Int, log10(niter))))
    maxiter > niter && println("\nIterations were rounded up to the next power of 10: $maxiter \n")

    #Create output file path
    defdir = string("aspu_results_", dtnow())
    outdir = out == "" ? defdir : (isdir(out) ? out : dirname(abspath(out)))
    isdir(outdir) || mkdir(outdir)
    defname = string("aspu_results_1e", ceil(Int, log10(maxiter)), "_", basename(filein))
    outname = (out == "" || isdir(out)) ? defname : basename(out)
    fileout = string(outdir, "/", outname)

    #Calculate R
    Σ::Matrix{Float64} = covfile == "" ?
        cov_io(filein; delim = delim, plim = plim, skip = skip, header = !noheader) :
        readdlm(covfile, ',')
    R = cov2cor(Σ)

    #Open input and output files
    fout = fileout == "-" ? stdout : open(fileout, "w")
    fin = filein == "-" ? stdin : open(filein, "r")

    #Write header
    line1 = noheader ? join(["snpid",[string("z",i) for i in 1:ntraits]...], delim) : readline(fin)
    join(fout, vcat(line1, "aspu_p", map(*, fill("p_spu_",length(pows)), string.(pows)), "gamma", '\n'), delim, "")

    #Show and save R
    if verbose
        println("Covariance matrix computed")
        display(R); println("")
    end
    outcov = string(outdir, "/aspu_z_covariance_", basename(filein))
    nosavecov || ( writedlm(outcov, R) )

    #Create MVM distribution from R
    mvn = invR_trans ? MvNormal(inv(R)) : MvNormal(R)
    trans = invR_trans ? inv(R) : one(R)
    ntraits = length(mvn)

    #Run simulations, and store forever
    aspu_obj = init_aspu_par(pows, mvn, maxiter; verbose=verbose)

    #Setup parallel channels
    buffer_s = min(10*nworkers(), outtest)
    jobs = RemoteChannel(()->Channel{String}(buffer_s))
    results = RemoteChannel(()->Channel{Any}(buffer_s))
    vars = RemoteChannel(()->Channel{Any}(length(workers())))

    #Start workers
    for p in workers()
        put!(vars, (aspu_obj, pows, delim, trans, na))
        remote_do(do_aspuwork, p, vars, jobs, results)
    end

    #Give workers a head start
    verbose && println_timestamp("Processing file...")
    buffer_n = 0
    for i in 1:buffer_s
        eof(fin) && break
        line = readline(fin)
        put!(jobs, line)
        buffer_n += 1
    end

    #loop through the whole file
    outtest2 = outtest - length(buffer_s)
    for (n, line) in enumerate(eachline(fin))
        n > outtest2 && break
        put!(jobs, line)
        out = take!(results)
        join(fout, vcat(out[1], out[2]..., '\n'), delim, "")
    end

    #clear out buffer
    for i in 1:buffer_n
        out = take!(results)
        join(fout, vcat(out[1], out[2]..., '\n'), delim, "")
    end

    #close files
    verbose && println_timestamp("All done.\n")
    close(fout)
    close(fin)
end
