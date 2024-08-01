using Documenter, Muscade

repo = normpath((@__DIR__)*"/..")
push!(LOAD_PATH,joinpath(repo,"src"))
cp(joinpath(repo,"LICENSE.md"),joinpath(repo,"docs/src/LICENSE.md"),force=true)
makedocs(sitename ="Muscade.jl",
        format    = Documenter.HTML(prettyurls = false,sidebar_sitename = false),
        pages     = ["index.md",
                     "Theory.md",
                     "Modelling.md",
                     "Solvers.md",
                     "Elements.md",
                     "ElementLibrary.md",
                     "TypeStable.md",
                     "Memory.md",
                     "Adiff.md",
                     "reference.md",
                     "LICENSE.md"],
                     source  = "src",
                     build   = "build"   
        )

#deploydocs(repo = "github.com/SINTEF/Muscade.jl.git",target="build",devbranch="philippe")

# https://sintef.github.io/Muscade.jl/v0.3.6/index.html
# https://sintef.github.io/Muscade.jl/stable/index.html

