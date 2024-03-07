
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
    anapot = AnalyticIsochrone()
    numpot = NumericalIsochrone()
    params = OrbitalParameters(Ω₀=Ω₀(numpot))
    # Considered locations
    a, e = 1.0, 0.05
    u = 0.4
    Ω1, Ω2 = frequencies_from_ae(a,e,anapot)

    println("ANALYTICAL COMPUTATIONS:")
    println("Frequencies computation Benchmark")
    @btime frequencies_from_ae($a,$e,$anapot)
    println("Frequencies inversion Benchmark")
    @btime ae_from_frequencies($Ω1,$Ω2,$anapot)
    println("NUMERICAL COMPUTATIONS:")
    println("Frequencies integrand computation Benchmark")
    @btime OrbitalElements.Θ($u,$a,$e,$numpot,$params)
    println("Frequencies computation Benchmark")
    @btime frequencies_from_ae($a,$e,$numpot,$params)
    println("Frequencies inversion Benchmark")
    @btime ae_from_frequencies($Ω1,$Ω2,$numpot,$params)
end

bench()