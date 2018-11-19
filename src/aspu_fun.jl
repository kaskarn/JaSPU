#fast spu(gamma)
function getspu!(spu, pows, z, n)
  for i in eachindex(pows)
    @inbounds pows[i] < 9 && (spu[i] = z[1]^(pows[i]))
    @inbounds pows[i] == 9 && (spu[i] = abs(z[1]))
  end
  for j = 2:n
    for i in eachindex(pows)
      @inbounds (pows[i] < 9) && (spu[i] += z[j]^(pows[i]))
      @inbounds (pows[i] == 9) && (abs(z[j]) > spu[i]) && (spu[i] = abs(z[j]))
    end
  end
  for i = eachindex(pows)
    @inbounds spu[i] = abs(spu[i])
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



function init_aspu(pows, mvn, ntest, maxiter; trans = Matrix{Float64}(I,length(mvn),length(mvn)))

    ntraits = length(mvn)
    ran = zeros(ntraits, ntest)
    ranspu = zeros(length(pows), ntest)
    tmp = zeros(length(pows))

    allvals = [zeros(max(ntraits, length(pows)), ntest*10) for i in 0:ceil(log10(maxiter/ntest))]

    #fill first bucket
    rand!(mvn, ran)
    for i in 1:ntest
        getspu!(view(allvals[1],:,i), pows, trans*ran[:,i], ntraits)
    end

    maxin = zeros(Int, ceil(Int, 1+log10(maxiter/ntest)))
    maxin_arr = zeros(Int, length(pows), ceil(Int, 1+log10(maxiter/ntest)))
    maxin[1] = ntest
    fill!(view(maxin_arr,:,1), ntest)

    for i in eachindex(allvals)[2:end]
        iternow = ntest*10^(i-1)
        toppows = [partialsort(view(allvals[i-1],j,:), Int(ntest*0.10); rev = true) for j in eachindex(pows)]
        x = create_arref!(ran, tmp, ntest, iternow, mvn, toppows, allvals, i, pows)
        maxin_arr[:,i] = x[2]
        maxin[i] = x[1]
    end

    allranks = [ zeros(Int, 2, n) for n in maxin ]
    allranks_f = [ zeros(Float64, 2, n) for n in maxin ]

    [ rank_spus!(allranks[i], allvals[i], maxin[i]) for i in eachindex(allvals)]
    allsorted = [[ sort(allvals[i][g,:],rev=true)[1:maxin_arr[g,i]] for g in eachindex(pows) ] for i in eachindex(allvals) ]
    allsorted_rank = [ sort(1 .- ((allranks[i][2,:] .+ 1) ./ maxin[i])) for i in eachindex(allvals) ]

    allranks_f = [ ((1 + maxin[i]) .- allranks[i][2,:])./maxin[i] for i in eachindex(allvals)]
    # allsorted, allsorted_rank, allranks_f, allranks, maxin_arr, maxin
    allsorted, allranks, maxin_arr, maxin
end


function getaspu(z, allsorted, allranks, maxin_arr, maxin, pows)
    spu = getspu(pows, z, length(z))
    minp = 0.0
    minp_f = 0.0
    gamma = 0
    pval = fill(0.0, size(allsorted[1],1))
    pval2 = copy(pval)
    ind_p = zeros(Float64, length(pows))
    pval_f = fill(0.0, size(allsorted[1],1))
    aspu_p = 1
    i = 0
    B = maxin[1]
    while minp_f < 0.085 && i < length(allsorted)
        i += 1
        for k in eachindex(pows)
            @inbounds pval[k] = sum( spu[k] .< allsorted[i][k][1:maxin_arr[k, i]] ) + 1
            pval_f[k] < 0.085 && (pval2[k] = pval[k])
        end
        minp, gamma = findmin(pval)
        pval_f = pval ./ (B + 1)
        minp_f = (minp)/(B + 1)
        for k in eachindex(ind_p)
            (ind_p[k] == 0 && pval_f[k] > 0.085) && (ind_p[k] = i)
        end
    end
    for j in allranks[i][2,:]
        @inbounds (maxin[i] - j + 1)/B < minp_f && (aspu_p += 1)
    end
    for k in eachindex(ind_p)
        ind_p[k] > 0 && (pval2[k] /= 10^(ind_p[k]-1)*B)
        ind_p[k] == 0 && (pval2[k] /= 10^(i-1)*B)
    end

    # aspu_p/(B+1)/10^(i-1), minp_f/10^(i-1), pval, pval_f, gamma, ind_p
    aspu_p/(B+1)/10^(i-1), pval2, gamma
end

function aspu(
    filename, outfile=""; covfile="", maxiter=Int(1e7), ntest=Int(1e4), pows=[1, 2, 9],
    plim = 1e-5, header = true, skip = 1, invcor=false, delim = '\t'
    )
    Σ, R = cor_io(filename; delim = delim)
    mvn = invcor ? MvNormal(inv(Σ)) : MvNormal(Σ)
    trans = invcor ? inv(Σ) : one(Σ)
    ntraits = length(mvn)

    allsorted, allranks, maxin_arr, maxin = init_aspu(pows, mvn, ntest, maxiter)

    fout = outfile == "" ? stdout : open(outfile, "w")
    f = open(filename, "r")

    #assumes header
    write(fout, readline(f))
    write(fout, delim)
    join(fout, vcat("aspu_p", map(*, fill("p_spu_",length(pows)), string.(pows)), "gamma"), delim)
    write(fout, '\n')

    for line in eachline(f)
        ls = split(line, delim)
        z = parse.(Float64, ls[2:end])
        aspu_out = getaspu(trans*z, allsorted, allranks, maxin_arr, maxin, pows)
        join(fout, vcat(ls, aspu_out..., '\n'), delim, "")
    end
    close(fout)
    close(f)
end
