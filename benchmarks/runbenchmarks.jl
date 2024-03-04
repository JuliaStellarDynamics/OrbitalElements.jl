
using BenchmarkTools
using Dates
using LibGit2
using OrbitalElements

function bench()
    # Date of the benchmark and configuration
    println("----------------------------")
    println("time (UTC): "*string(now(UTC)))
    println("git commit: "*LibGit2.head((@__DIR__)*"/.."))

    # Considered potential (isochrone)
    pot = IsochronePotential()
    params = OrbitalParameters(Ω₀=Ω₀(pot))
    # Considered locations
    a, e = 1.0, 0.05
    u = 0.4
    Ω1, Ω2 = ComputeFrequenciesAE(pot,a,e)

    println("ANALYTICAL COMPUTATIONS:")
    println("Frequencies computation Benchmark")
    @btime ComputeFrequenciesAE($pot,$a,$e)
    println("Frequencies inversion Benchmark")
    @btime ComputeAEFromFrequencies($pot,$Ω1,$Ω2)
    println("NUMERICAL COMPUTATIONS:")
    println("Frequencies integrand computation Benchmark")
    @btime OrbitalElements.ΘAE($pot,$u,$a,$e,$params)
    println("Frequencies computation Benchmark")
    @btime ComputeFrequenciesAE($pot,$a,$e,$params)
    println("Frequencies inversion Benchmark")
    @btime ComputeAEFromFrequencies($pot,$Ω1,$Ω2,$params)
    _, _, niter, _ = OrbitalElements.AEFromΩ1Ω2Brute(Ω1,Ω2,pot,params)
    println("Number of iterations : ",niter)
end

bench()