using Documenter, EspyInsideFunction

push!(LOAD_PATH,"../src")
makedocs(sitename ="Muscade.jl",
#        modules   = [EspyInsideFunction],
        format    = Documenter.HTML(prettyurls = false),
        pages     = ["index.md"],
        source  = "src",
        build   = "build"   
        )

here = @__DIR__
mv(here*"\\build\\index.html",here*"\\index.html", force=true)        
mv(here*"\\build\\search.html",here*"\\search.html", force=true)        
mv(here*"\\build\\search_index.js",here*"\\search_index.js", force=true)        
mv(here*"\\build\\assets",here*"\\assets", force=true)        
#deploydocs(repo = "github.com/PhilippeMaincon/Muscade.jl.git")

