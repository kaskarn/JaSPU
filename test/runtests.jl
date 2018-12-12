println("Testing...")
using Test
using Distributed
using Distributions
using Random

@everywhere using JaSPU

#generic tests and benchmarks
aspu("data/smalldat", 10^5; out = "results/tmptest")
aspu("data/smalldat", 10^6; out = "results", plim = 10^-4)

@time aspu("data/smalldat", 10^5; out = "results")
@time aspu("data/smalldat", 10^6; out = "results", verbose = false)

#test values
Random.seed!(0809)
pows = collect(0:8)
R = [1 0.2; 0.2 1]
z = [2, -1]
mvn = MvNormal(R)
aspuobj = JaSPU.init_aspu_par(pows, mvn, 10^5; verbose=false)
out = JaSPU.getaspu(z, aspuobj..., pows)

# z = [2, -1, 3, 8]
# JaSPU.getspu(pows, z)
# getspu(pows, z, 4)

@test out[1] == 0.12266877331226687
