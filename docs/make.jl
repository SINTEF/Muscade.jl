using Documenter, Muscade

repo = normpath((@__DIR__)*"/..")
push!(LOAD_PATH,joinpath(repo,"src"))
cp(joinpath(repo,"LICENSE.md"),joinpath(repo,"docs/src/LICENSE.md"),force=true)
makedocs(sitename ="Muscade.jl",
        format    = Documenter.HTML(prettyurls = false,sidebar_sitename = false),
        pages     = ["index.md",
                     "README.md",
                     "Theory.md",
                     "Modelling.md",
                     "Elements.md",
                     "Solvers.md",
                     "Builtin_els.md",
                     "Builtin_sol.md",
                     "TypeStable.md",
                     "Memory.md",
                     "Adiff.md",
                     "reference.md",
                     "LICENSE.md"],
                     source  = "src",
                     build   = "build"   
        )

deploydocs(repo = "github.com/SINTEF/Muscade.jl.git",target="build",devbranch="philippe")

# for f ∈ readdir(joinpath(repo,"docs/build")) # needed so that the doc will deploy on 
# mv(joinpath(repo,"docs/build/",f), joinpath(repo,"docs",f) , force=true)        
# end       
# deploydocs(repo = "github.com/SINTEF/Muscade.jl.git")

# https://sintef.github.io/muscade.jl/dev    → philippe
# https://sintef.github.io/muscade.jl/stable → main 
# https://sintef.github.io/muscade.jl/v0.3.3 → ? 
# https://sintef.github.io/Muscade.jl/       → README.html