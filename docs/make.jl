# Adding the package src to the load path
push!(LOAD_PATH,"../src/")

using Documenter, OrbitalElements

makedocs(sitename = "OrbitalElements.jl",
         pages=[
                "Home" => "index.md",
                "Functions" => "functions.md"
               ],
         format = Documenter.HTML(prettyurls=false))
deploydocs(repo="github.com/michael-petersen/OrbitalElements.jl")