# Adaptive Sum of Powered Tests (aSPU) in Julia

## Installation

The package may be installed from github:

```julia
(v1.0) pkg> add https://github.com/kaskarn/JaSPU
```
or
```julia
julia> using Pkg
julia> Pkg.add("https://github.com/kaskarn/JaSPU")
```

## The Adaptive Sum of Powered Tests

The Adaptive Sum of Powered Tests (aSPU) is used in genome-wide association
settings to evaluate the effect of SNPs across k traits. To do so, it
uses k z-scores from previous regression analyses, and aggregates them as
multiple sums of powered scores:

<a href="https://www.codecogs.com/eqnedit.php?latex=SPU(\gamma)=\sum_k&space;|S_k^\gamma&space;|" target="_blank"><img src="https://latex.codecogs.com/gif.latex?SPU(\gamma)=\sum_k&space;|S_k^\gamma&space;|" title="SPU(\gamma)=\sum_k |S_k^\gamma |" /></a>

Where S may either be untransformed z-scores (the default), or z-scores transformed with the inverse of R:
<a href="https://www.codecogs.com/eqnedit.php?latex=S=R^{-1}Z" target="_blank"><img src="https://latex.codecogs.com/gif.latex?S=R^{-1}Z" title="S=R^{-1}Z" /></a>, and gamma takes on integer values, by default:
<a href="https://www.codecogs.com/eqnedit.php?latex=\gamma&space;=&space;1,2,\ldots,\infty" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\gamma&space;=&space;1,2,\ldots,\infty" title="\gamma = 1,2,\ldots,\infty" /></a>
with <a href="https://www.codecogs.com/eqnedit.php?latex=SPU(\infty)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?SPU(\infty)" title="SPU(\infty)" /></a>
equivalently computed as <a href="https://www.codecogs.com/eqnedit.php?latex=SPU(\infty)=\max_k&space;|S_k|" target="_blank"><img src="https://latex.codecogs.com/gif.latex?SPU(\infty)=\max_k&space;|S_k|" title="SPU(\infty)=\max_k |S_k|" /></a>

aSPU adaptively selects the SPU with the greatest power, and performs Monte-Carlo simulations to produce p-values,
using the empirical multivariate-normal distribution of null z-scores.

## Package description

The package provides the function `aspu` to compute SPU and aSPU p-values for each SNP contained in an input file. This implementation
fully uses available processors (which must be previously added with `addprocs()`, or `ClusterManagers` functions), has a minimal
memory footprint, and is orders of magnitude faster than R implementations.

### Input file
The input file must include one SNP per line, with the first column containing SNP names, and, by default, each subsequent column
containing  a z-score value; this default behavior can be changed using the option skip = n, so that z-scores start at the nth
column, instead of the second one. Lines with missing values ("NA" by default, or specified with the option `na = string`) are
permitted, and will appear in the output with missing entries.

### Output file
The output file uses the same delimiter as the input file. It replicates the input files, and adds columns for the aspu p-value,
the p-values for each SPU score (each gamma), and the best-powered gamma value (in the case of ties, the higher gamma value is
returned).

By default, the output file is named following the pattern "aspu_results_1eN_filein", where N is the number of iterations,
and filein is the input filename. By default, the output file is placed in a timestamped folder created in the working directory.
These behaviors can be overriden by specifying `out = path`, where `path` can be a directory (to override the location of results), or
a file name.

### Usage

The `aspu` function has two required arguments: a path to the input file, and the number of iterations used to calculate p-values.
Since Monte-Carlo simulations are used to compute p-values, the minimum achievable aSPU p-value using N iterations will be <a href="https://www.codecogs.com/eqnedit.php?latex=p_{min}=(N&plus;1)^{-1}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?p_{min}=(N&plus;1)^{-1}" title="p_{min}=(N+1)^{-1}" /></a>

*The number of iterations should be provided as a power of 10, any other number will be rounded up to the next power of
10 (e.g. 200,000 -> 1,000,000)*

```julia
   aspu(
    filein::AbstractString, maxiter::Int64;                                   #required options
    pows::Vector{Int64} = collect(0:8), invR_trans::Bool = false,             #key aSPU parameters
    covfile::AbstractString = "", plim::Float64 = 1e-4,                       #R estimation options
    delim::Char = '\t', noheader::Bool = false, skip::Int64 = 1,              #input file options
    out::AbstractString = "", verbose::Bool = true, nosavecov::Bool = false,  #output options
    outtest::Real = Inf                                                       #testing/development
    )
```

A simple, common example would be to compute aSPU with default gammas, allowing p-values to well exceed the genome-wide significance threshold of 5e-8:

```julia
   aspu("myzscores.txt", 10^9) #use default output
   aspu("myzscores.txt", 10^9, out = "aspu_results/myresults.txt") #save output to specific destination
   aspu("myzscores.txt", 10^9, plim = 1e-5) #use a lower threshold for null SNPs when computing Z correlation
```

Users may choose a more succinct set of gamma values, since added gains from large gammas are uncertain. A reasonable alternative set may be 1, 2, 3, and infinity. We use 0 to represent infinity, and the `aspu` call can be made as:

```julia
   aspu("myzscores.txt", 10^9, pows = [0, 1, 2, 3])
```

## Performance
Below are rough estimates for the CPU-hours spent on Monte-Carlo simulations, estimated on the UNC-Chapel hill high-performance computing cluster (longleaf):

| Iterations | minimum p-value | CPU-hours |
| ---------- | --------------- | --------- |
| 10^9       | 1E-9            | 2.8       |
| 10^10      | 1E-10           | 24        |
| 10^11      | 1E-11           | 264       |

These do not count the time taken to read and process SNPs, which depends on the number of SNPs. In our case, it took 5.5 hours to run `aspu` on 18M SNPs, using 10^11 iterations, and 60 CPUs.

### Memory usage
The `aspu` function does not create particularly large objects
