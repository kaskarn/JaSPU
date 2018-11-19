#computes X Xt without loading file in memory
function xxt_io(filename::AbstractString; plim = 1e-5, delim::Char = '\t', skip = 1, header = true)
    print(delim)
    f = open(filename, "r")
    xlen = length(split(readline(f), delim)) - skip
    header || seek(f, 0)
    xxt = zeros(Float64, xlen, xlen)
    n, s = 0, 0
    ztresh = abs(quantile(Normal(0,1), plim/2))

    for line in eachline(f)
        lp = parse.(Float64, split(line, delim)[2:end])
        if maximum(abs.(lp)) > ztresh
            s += 1
            continue
        end
        for i in eachindex(lp)
            for j in 1:i
                xxt[i,j] += lp[i]*lp[j]
            end
        end
        n += 1
    end
    for j in 1:xlen
        for i in (j+1):xlen
            xxt[j,i] = xxt[i,j]
        end
    end

    close(f)
    xxt, n, s
end

#compute cov(X) and cor(X) without loading file in memory
function cor_io(filename::AbstractString; delim::Char = '\t', covfile="")
    covfile=="" || return readdlm(open(covfile,"r"), delim)
    xxt, n, s = xxt_io(filename; delim=delim)
    Σ = xxt ./ n
    R = Σ * inv(Diagonal(Σ))

    Σ
end
