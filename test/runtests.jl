println("Testing...")
using Distributed
# addprocs()

# @everywhere begin
#     using JaSPU
#     using Distributions
#     using CSV
#     using DelimitedFiles
#     using Random
#     using LinearAlgebra
# end

@everywhere using JaSPU


@time aspu("data/smalldat", "data/tmptest"; maxiter=Int(1e4))
@time aspu("data/smalldat", "data/tmptest"; maxiter=Int(1e5), verbose = false, plim = 1e-4)
@time aspu("data/smalldat", "data/tmptest"; maxiter=Int(1e6), verbose = false)
@time aspu("data/smalldat", maxiter=Int(1e6))
# @time aspu("data/smalldat", "data/tmptest"; maxiter=Int(1e7), verbose = false)
# @time aspu("data/testdat", "data/tmptest"; maxiter=Int(1e8), outtest = 10) #82 sec
# @time aspu("data/testdat", "data/tmptest"; maxiter=Int(1e9), outtest = 10) #800sec
# @time aspu("data/smalldat", "data/tmpsmall"; maxiter=Int(1e5))


# println()
# rmprocs(workers())
# aspu("data/testdat", "data/tmptest"; maxiter=Int(1e5), outtest = 10)
# @time aspu("data/testdat", "data/tmptest"; maxiter=Int(1e5), outtest = 10)
# @time aspu("data/testdat", "data/tmptest"; maxiter=Int(1e6), outtest = 10)
# @time aspu("data/testdat", "data/tmptest"; maxiter=Int(1e7), outtest = 10)


# @time aspu("data/testdat", "data/tmptest"; maxiter=Int(1e8), outtest = 10)
# @time aspu("data/testdat", "data/tmpout"; maxiter=Int(1e5), invcor=true, covfile="data/covmat")
# @time aspu("data/tmp", "data/tmpout"; maxiter=Int(1e5))
# @time aspu("data/tmp", "data/tmpout"; maxiter=Int(1e5))


# @time aspu("tmp", "tmpout_1e9", maxiter=Int(1e9))
