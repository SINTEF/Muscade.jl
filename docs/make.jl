
#Assumes that pwd() is docs
docs        = @__DIR__
muscade     = normpath(joinpath(docs,".."))
docsrc      = joinpath(docs,"src")
examplesrc(ex) = normpath(joinpath(docs,"..","examples",ex))
examples    = ["StaticBeamAnalysis","ModalBeamAnalysis","DynamicBeamAnalysis","DecayAnalysis","DryFriction"]

requiredIncludeFiles = ["toolbox" "BeamElement.jl"; 
                        "toolbox" "Rotations.jl";  
                        "examples" "SCR.csv"] 
 for idx ∈ 1:size(requiredIncludeFiles,1)
        cp(     joinpath(muscade,requiredIncludeFiles[idx,1],requiredIncludeFiles[idx,2]),
                joinpath(muscade,"docs","src",requiredIncludeFiles[idx,2]),force=true)             
 end

using Pkg
Pkg.activate(docs)

using Muscade, Muscade.Toolbox
using Documenter, Literate, DocumenterCitations, Printf, Interpolations

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
        modules   = [Muscade, Muscade.Toolbox],
        doctest   = false,
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
                        "Diagnostic.md",
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
