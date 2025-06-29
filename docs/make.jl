
#Assumes that pwd() is docs
docs        = @__DIR__
muscade     = normpath(joinpath(docs,".."))
docsrc      = joinpath(docs,"src")
examplesrc(ex) = normpath(joinpath(docs,"..","examples",ex))
examples    = ["BeamElements","StaticBeamAnalysis","DynamicBeamAnalysis","DecayAnalysis","DryFriction"]

requiredIncludeFiles = ["BeamElements.jl","Rotations.jl"]
 for includeFile ∈ requiredIncludeFiles
        cp(joinpath(muscade,"examples",includeFile),joinpath(muscade,"docs","src",includeFile),force=true)             
 end


using Pkg
Pkg.activate(docs) 

using Muscade
using Documenter, Literate, DocumenterCitations,Printf

cp(joinpath(muscade,"LICENSE.md"),joinpath(docsrc,"LICENSE.md"),force=true)

## Literate

nodocstr(str) =  replace(str, r"(*ANYCRLF)^\"\"\"$.*?^\"\"\"$"ms => "") # Take """ somedocstring """ out of str

@printf "\nLiterate.markdown: *.jl → *.md\n\n"
for ex ∈ examples
        Literate.markdown(examplesrc(@sprintf("%s.jl",ex)), docsrc, preprocess = nodocstr)
end

## DocumenterCitations

bib = CitationBibliography(joinpath(docsrc, "ref.bib"); style=:authoryear) #:numeric

## Documenter

@printf "\nDocumenter.makedocs: *.md → *.html\n\n"
makedocs(sitename ="Muscade.jl",
        modules   = [Muscade],
        doctest   = false, # we do not use doctest, we run Literate.jl on mydemo.jl files that are also included in unit test files
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


deploydocs(repo = "github.com/SINTEF/Muscade.jl.git",devbranch="dev",   devurl="dev",   versions = ["stable" => "stable", "dev" => "dev"])
deploydocs(repo = "github.com/SINTEF/Muscade.jl.git",devbranch="main",  devurl="stable",versions = ["stable" => "stable", "dev" => "dev"])

# https://sintef.github.io/Muscade.jl/dev
# https://sintef.github.io/Muscade.jl/stable
