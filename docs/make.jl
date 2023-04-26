using Documenter, Muscade

push!(LOAD_PATH,"../src")
cp("./LICENSE.md","./docs/src/LICENSE.md",force=true)
makedocs(sitename ="Muscade.jl",
        format    = Documenter.HTML(prettyurls = false),
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
                     "reference.md",
                     "LICENSE.md"],
        source  = "src",
        build   = "build"   
        )

#here = @__DIR__
# mv(here*"\\build\\index.html",here*"\\index.html", force=true)        
# mv(here*"\\build\\search.html",here*"\\search.html", force=true)        
# mv(here*"\\build\\search_index.js",here*"\\search_index.js", force=true)        
# mv(here*"\\build\\assets",here*"\\assets", force=true)        
#deploydocs(repo = "github.com/PhilippeMaincon/Muscade.jl.git")

