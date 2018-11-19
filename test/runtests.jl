println("Testing...")

using JaSPU

using Distributions
using Distributed
using CSV
using DelimitedFiles
using Random
using LinearAlgebra

@time aspu("data/testdat", "data/tmpout"; maxiter=Int(1e5), invcor=true, covfile="data/covmat")
@time aspu("tmp", "data/tmpout"; maxiter=Int(1e5))
@time aspu("data/tmp", "data/tmpout"; maxiter=Int(1e5))


# @time aspu("tmp", "tmpout_1e9", maxiter=Int(1e9))
