nodocstr(str) =  replace(str, r"(*ANYCRLF)^\"\"\"$.*?^\"\"\"$"ms => "")

docs    = @__DIR__
#muscade = dirname(pathof(Muscade))
muscade = normpath(joinpath(docs,".."))
docsrc  = joinpath(docs,"src")

Pkg.activate(docs) 

using Muscade, Documenter, Literate, DocumenterCitations,Printf

cp(joinpath(muscade,"LICENSE.md"),joinpath(docsrc,"LICENSE.md"),force=true)

#Literate.markdown(joinpath(docsrc,"tutorial1.jl"),docsrc)
els = ["DryFriction","BeamElement"]
for el ∈ els
        Literate.markdown(joinpath(muscade,"src","Elements",@sprintf("%s.jl",el)),docsrc,
                       execute=false,codefence= "````julia" => "````", # prevent execution of the code by Literate and Documenter
                       preprocess = nodocstr) 
end

bib = CitationBibliography(joinpath(docsrc, "ref.bib"); style=:authoryear) #:numeric

makedocs(sitename ="Muscade.jl",
        # modules = [Elements],
        format    = Documenter.HTML(prettyurls = false,sidebar_sitename = false),
        plugins   = [bib],
        pages     = ["index.md",
                     "Theory.md",
                     "User manual" => ["Modelling.md",
                                        "Solvers.md",
                                        "Elements.md"],
                     "reference.md",
                     "Appendix" => [
                                            "Element code examples" => [@sprintf("%s.md",el) for el∈els],
                                            "TypeStable.md",
                                            "Memory.md",
                                            "Adiff.md",
                                            "litterature.md"],
#                     "Demo tutorial" => "tutorial1.md",  # source for this is Muscade/docs/src/tutorial1.jl NB: *.jl
                     "LICENSE.md"],
                     source  = "src",
                     build   = "build"                 
        )

#deploydocs(muscade = "github.com/SINTEF/Muscade.jl.git",target="build",devbranch="dev")

# https://sintef.github.io/Muscade.jl/v0.3.6/index.html
# https://sintef.github.io/Muscade.jl/stable/index.html

Pkg.activate(muscade)
