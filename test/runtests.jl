println("Testing...")
using Test
using Distributed
using Distributions
using Random

@everywhere using JaSPU

#generic tests and benchmarks
aspu("data/smalldat", "results/tmptest"; maxiter=10^5)
aspu("data/smalldat", "results", maxiter=10^6, plim = 10^-4)

@time aspu("data/smalldat", "results"; maxiter=10^5, verbose = false)
@time aspu("data/smalldat", "results"; maxiter=10^6, verbose = false)


#test values
Random.seed!(0809)
pows = collect(0:8)
R = [1 0.2; 0.2 1]
z = [2, -1]
mvn = MvNormal(R)
aspuobj = JaSPU.init_aspu_par(pows, mvn, 10^4, 10^5; verbose=false)
out = JaSPU.getaspu(z, aspuobj..., pows)


@test out[1] == 0.12266877331226687
