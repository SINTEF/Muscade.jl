docs    = @__DIR__
muscade = normpath(joinpath(docs,".."))
docsrc  = joinpath(docs,"src")

using Pkg

cd(muscade)
Pkg.activate(".") 
using Muscade

cd(docs)
Pkg.activate(".") 
using Documenter, Literate, DocumenterCitations

cp(joinpath(muscade,"LICENSE.md"),joinpath(docsrc,"LICENSE.md"),force=true)

Literate.markdown(joinpath(docsrc,"tutorial1.jl"),docsrc)

bib = CitationBibliography(
    joinpath(docsrc, "ref.bib");
    style=:authoryear #:numeric
)

makedocs(sitename ="Muscade.jl",
        # modules = [Elements],
        format    = Documenter.HTML(prettyurls = false,sidebar_sitename = false),
        plugins=[bib],
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
                     "litterature.md",
                     "Demo tutorial" => "tutorial1.md",
                     "LICENSE.md"],
                     source  = "src",
                     build   = "build"                 
        )

#deploydocs(muscade = "github.com/SINTEF/Muscade.jl.git",target="build",devbranch="dev")

# https://sintef.github.io/Muscade.jl/v0.3.6/index.html
# https://sintef.github.io/Muscade.jl/stable/index.html

cd(muscade)
Pkg.activate(".")
