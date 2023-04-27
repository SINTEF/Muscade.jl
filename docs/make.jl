using Documenter, Muscade

push!(LOAD_PATH,"../src")
cp("./LICENSE.md","./docs/src/LICENSE.md",force=true)
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

# mv("./build/index.html"     ,"./index.html"     , force=true)        
# mv("./build/search.html"    ,"./search.html"    , force=true)        
# mv("./build/search_index.js","./search_index.js", force=true)        
# mv("./build/assets"         ,"./assets"         , force=true)        
deploydocs(repo = "github.com/SINTEF/Muscade.jl.git")
