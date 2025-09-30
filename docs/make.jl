# Adding the package src to the load path
push!(LOAD_PATH,"../src/")

using Documenter
using OrbitalElements

# Use the readme as the installation page
# Open  readme page
readme = open(dirname(@__FILE__) * "/../README.md") do io
      read(io, String)
end
# Remove the images by which are not under docs/ folder
# and fail to render
ref1="![`Isochrone frequencies`](examples/IsochroneFrequencies_original.png)"
ref2="![`Forward+backward errors`](examples/ForwardBackwardErrors_original.png)"
readme = replace(readme, ref1 => s"")
readme = replace(readme, ref2 => s"")
# Same with the notebook reference
readme = replace(readme, "[`examples/test_OrbitalElements.ipynb`](examples/test_OrbitalElements.ipynb)" => "examples/test_OrbitalElements.ipynb")
# Paste it in installation page
open(dirname(@__FILE__) * "/src/installation.md", "w") do io
      write(io, readme)
end

makedocs(sitename = "OrbitalElements.jl",
         pages=[
                "Home" => "index.md",
                "Quick start" => "installation.md",
                "Manual" => "manual.md",
                "Formulae" => "formulae.md",
                "References" => "functions.md"
               ],
         format = Documenter.HTML(prettyurls=false))

deploydocs(repo="github.com/JuliaStellarDynamics/OrbitalElements.jl")