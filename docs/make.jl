nodocstr(str) =  replace(str, r"(*ANYCRLF)^\"\"\"$.*?^\"\"\"$"ms => "")

#Assumes that pwd() is docs
docs    = @__DIR__
muscade = normpath(joinpath(docs,".."))
docsrc  = joinpath(docs,"src")
examplessrc = normpath(joinpath(docs,"..","examples"))

using Pkg
Pkg.activate(docs) 

using Muscade
using Documenter, Literate, DocumenterCitations,Printf

cp(joinpath(muscade,"LICENSE.md"),joinpath(docsrc,"LICENSE.md"),force=true)

## Literate

examples = ["StaticBeamAnalysis","DecayAnalysis","DryFriction"]

function replace_includes(str)
        included = ["BeamElements.jl","Rotations.jl"] # in this order
        path = "examples/"
        for ex in included
                content = read(path*ex, String)
                str     = replace(str, "include(\"$(ex)\")" => content)
        end
        return nodocstr(str)
end

for ex ∈ examples
        Literate.markdown(joinpath(examplessrc,@sprintf("%s.jl",ex)),docsrc, preprocess = replace_includes)
end

# els = ["DryFriction","BeamElements"]
# for el ∈ els
#         Literate.markdown(joinpath(muscade,"examples",@sprintf("%s.jl",el)),docsrc,
#                        execute=false,codefence= "````julia" => "````", # prevent execution of the code by Literate and Documenter
#                        preprocess = replace_includes) 
# end

## DocumenterCitations

bib = CitationBibliography(joinpath(docsrc, "ref.bib"); style=:authoryear) #:numeric

## Documenter

makedocs(sitename ="Muscade.jl",
        modules   = [Muscade],
        format    = Documenter.HTML(    prettyurls          = false,
                                        sidebar_sitename    = false,
                                        size_threshold_warn = 256*1024,
                                        # example_size_threshold = nothing,                 
                                        size_threshold      = nothing),
        plugins   = [bib],
                pages= ["index.md",
                        "Theory.md",
                        "User manual" => [      "Modelling.md",
                                                "Solvers.md",
                                                "Elements.md"],
                        "Examples" => [@sprintf("%s.md",ex) for ex∈examples], 
                        "reference.md",
                        "Appendix" => [
                                                "TypeStable.md",
                                                "Memory.md",
                                                "Adiff.md",
                                                "litterature.md"],
                        "LICENSE.md"],
                        source  = "src",
                        build   = "build"                 
        )

#deploydocs(muscade = "github.com/SINTEF/Muscade.jl.git",target="build",devbranch="dev")

# https://sintef.github.io/Muscade.jl/v0.3.6/index.html
# https://sintef.github.io/Muscade.jl/stable/index.html

Pkg.activate(muscade)
