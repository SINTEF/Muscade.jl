nodocstr(str) =  replace(str, r"(*ANYCRLF)^\"\"\"$.*?^\"\"\"$"ms => "")

#Assumes that pwd() is docs
docs    = @__DIR__
muscade = normpath(joinpath(docs,".."))
docsrc  = joinpath(docs,"src")
examplessrc = normpath(joinpath(docs,"..","examples"))

using Pkg
Pkg.activate(docs) 

using Muscade, Muscade.BeamElements, Muscade.SdofElements
using Documenter, Literate, DocumenterCitations,Printf

cp(joinpath(muscade,"LICENSE.md"),joinpath(docsrc,"LICENSE.md"),force=true)

examples = ["StaticAnalysisBeam","DecayAnalysis"]
for ex ∈ examples
        Literate.markdown(joinpath(examplessrc,@sprintf("%s.jl",ex)),docsrc)
end

els = ["SdofElements","BeamElements"]
for el ∈ els
        Literate.markdown(joinpath(muscade,"src","Elements",@sprintf("%s.jl",el)),docsrc,
                       execute=false,codefence= "````julia" => "````", # prevent execution of the code by Literate and Documenter
                       preprocess = nodocstr) 
end

bib = CitationBibliography(joinpath(docsrc, "ref.bib"); style=:authoryear) #:numeric

makedocs(sitename ="Muscade.jl",
        modules   = [Muscade, Muscade.BeamElements, Muscade.SdofElements],
        format    = Documenter.HTML(    prettyurls          = false,
                                        sidebar_sitename    = false,
                                        size_threshold_warn = 256*1024,
                                        # example_size_threshold = nothing,                 
                                        size_threshold      = nothing),
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
                    "Examples" => [@sprintf("%s.md",ex) for ex∈examples], 
                     "LICENSE.md"],
                     source  = "src",
                     build   = "build"                 
        )

#deploydocs(muscade = "github.com/SINTEF/Muscade.jl.git",target="build",devbranch="dev")

# https://sintef.github.io/Muscade.jl/v0.3.6/index.html
# https://sintef.github.io/Muscade.jl/stable/index.html

Pkg.activate(muscade)
