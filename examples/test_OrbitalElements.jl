########################################
# Installing the necessary libraries
########################################
using Pkg
Pkg.add(["Plots","LaTeXStrings"])

########################################
# Importing the necessary libraries
########################################
using OrbitalElements
using Plots
using LaTeXStrings

########################################
# Handling directory location
########################################
cd(@__DIR__()) # Setting pwd to the directory where the example file is located

########################################
# Simple Frequencies
########################################
# Creating the potential
println("Creating the potential ... ")
G, M, bc = 1.0, 1.0, 1.0
model = NumericalIsochrone(G=G, M=M, bc=bc)
analytic_model = AnalyticIsochrone(G=G, M=M, bc=bc)

# Creating the parameter structure
params = OrbitalParameters(rc=radial_scale(model))

# Points (line of constant semi-major axis)
println("Compute frequencies ... ")
a = 1.0
ne = 200
emin, emax = 0., 1.
tabe = collect(LinRange(emin,emax,ne))

# Compute the values of the frequencies and store them
tabΩ = zeros(2,ne)      # Storage for the numerical frequencies
tabΩIso = zeros(2,ne)   # Storage for the analytic frequencies
for j = 1:ne
    # Compute the frequencies at (a,e)
    Ω1, Ω2 = frequencies_from_ae(a, tabe[j], model, params)
    Ω1true, Ω2true = frequencies_from_ae(a, tabe[j], analytic_model)
    # Store them
    tabΩ[1,j], tabΩ[2,j] = Ω1, Ω2
    tabΩIso[1,j], tabΩIso[2,j] = Ω1true, Ω2true
end

# Plot the curves
println("Plotting ... ")
p1 = plot(tabe, tabΩ[1,:], xlabel=L"e", ylabel=L"\Omega_1", label="Numerical")
plot!(tabe, tabΩIso[1,:], label="True")
p2 = plot(tabe, tabΩ[2,:], xlabel=L"e", ylabel=L"\Omega_2", label=false)
plot!(tabe, tabΩIso[2,:], label=false)
plΩ=plot(
    p1, p2;
    layout=(1,2), 
    size=(1000,400),
    plot_title="Isochrone frequencies",
    left_margin = 5Plots.mm,
    bottom_margin = 5Plots.mm
)
savefig(plΩ,"IsochroneFrequencies.png")
println("The plot has been saved in the same folder as this example script 
    under the name 'IsochroneFrequencies.png'.")

########################################
# Error on forward+backward mapping
########################################

# Semi-major axis and eccentricity domains
println("Computing errors on forward/backward mappings to frequencies...")
na, ne = 100, 101
aminexp, amaxexp = -2, 2
emin, emax = 0.0, 1.0
tabaexp = collect(LinRange(aminexp, amaxexp, na))
taba = 10 .^tabaexp
tabe = collect(LinRange(emin, emax, ne))

errtabAE = zeros(na, ne)

for ka = 1:na
    loca = taba[ka] 
    for ke = 1:ne
        e = tabe[ke]
        # Forward mapping
        Ω1, Ω2 = frequencies_from_ae(loca, e, model, params)
        # Backward mapping
        aback, eback = ae_from_frequencies(Ω1, Ω2, model, params)
        # Error
        errtabAE[ka,ke] = abs((aback - loca) / loca) + abs(eback - e)
    end
end 

# Error made on (a,e)->mapped variables->(a,e) in log scale
println("Plotting ... ")
plerr = heatmap(
    taba, tabe, transpose(log10.(errtabAE));
    xscale=:log10,
    colorbar_ticks=(-14:2:-4,string.(-14:2:-4)),
    colormap=:viridis,xlabel=L"a",
    ylabel=L"e",clim=(-14,-4),
    colorbar_title=" \n"*L"\log_{10} (\mathrm{error})",
    title="Error on successive \n forward/backward mappings",
    size=(400,400),
    right_margin=7Plots.mm
)
# colorbar_ticks is not available for GR (default) backend for now
# However it does not prevent the figure plotting.
savefig(plerr,"ForwardBackwardErrors.png")
println("The plot has been saved in the same folder as this example script 
    under the name 'ForwardBackwardErrors.png'.")