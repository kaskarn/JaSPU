println("Testing...")
using Distributed

@everywhere using JaSPU

aspu("data/smalldat", "data/tmptest"; maxiter=Int(1e5))
aspu("data/smalldat", maxiter=Int(1e6))

println("\n\n")

@time aspu("data/smalldat", "data/tmptest"; maxiter=Int(1e5), verbose = false, plim = 1e-4)
@time aspu("data/smalldat", "data/tmptest"; maxiter=Int(1e6), verbose = false)
