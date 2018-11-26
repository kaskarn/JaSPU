module JaSPU

using Distributions
using Distributed
using DelimitedFiles
using Random
using LinearAlgebra
using Dates

export
    cov_io,
    aspu

include("utils.jl")
include("covmat.jl")
include("aspu_fun.jl")


end # module
