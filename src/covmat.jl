#computes X Xt without loading file in memory
function xxt_io(filename::AbstractString; plim = 1e-5, delim::Char = '\t', skip = 1, header = true, na = "NA")
    f = open(filename, "r")
    xlen = length(split(readline(f), delim)) - skip
    header || seek(f, 0)
    xxt = zeros(Float64, xlen, xlen)
    tot = zeros(Float64, xlen)
    n, s = 0, 0, 0
    ztresh = abs(quantile(Normal(0,1), plim/2))
    for line in eachline(f)
        ls = split(line, delim)[(1+skip):end]
        in(na, ls) && continue
        lp = parse.(Float64, ls)
        if maximum(abs.(lp)) > ztresh
            s += 1
            continue
        end
        n += 1
        for i in eachindex(lp)
            tot[i] += lp[i]
            for j in 1:i
                xxt[i,j] += lp[i]*lp[j]
            end
        end
    end
    for j in 1:xlen
        for i in (j+1):xlen
            xxt[j,i] = xxt[i,j]
        end
    end

    close(f)
    xxt, tot, n, s
end

#compute cov(X) from sums of squares computed with xxt() -- without loading file in memory
function cov_io(fin; plim = 1e-5, delim::Char = '\t', skip = 1, header = true)
    xxt, tot, n, s = xxt_io(fin; plim = plim, delim=delim, header = header, skip = skip)
    μ = tot./n
    Σ = (xxt .- μ*tot' .- tot*μ' .+  n*μ*μ')/(n-1)
end
