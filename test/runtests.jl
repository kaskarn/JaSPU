println("Testing...")

using JaSPU

using Distributions
using Distributed
using CSV
using DelimitedFiles
using Random
using LinearAlgebra

@time aspu("data/testdat", "tmpout", maxiter=Int(1e5), invcor=true, covfile="data/covmat")
@time aspu("tmp", "tmpout", maxiter=Int(1e5))
@time aspu("data/tmp", "data/tmpout", maxiter=Int(1e5))


# @time aspu("tmp", "tmpout_1e9", maxiter=Int(1e9))

#
# addprocs(2)
# @everywhere using Distributions
# @everywhere include("C:/Users/balda/Documents/GitHub/JaSPU/src/types.jl")
# include("mkcor.jl")
#
# a = @everywhere pwd()
