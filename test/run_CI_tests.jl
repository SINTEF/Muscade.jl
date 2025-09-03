# Runs tests, examples and documentation generation 
# from a clean environment, using the version of Muscade
# located in "muscadeDir". Useful to run prior to 
# a merge.

muscadeDir = joinpath(@__DIR__,"..")

using Pkg

function devMuscadeIn(workingDir,muscadeDir)
	cd(workingDir)
	isfile("Manifest.toml") && rm("Manifest.toml")    
	Pkg.activate(".")
    try Pkg.rm("Muscade") catch err end
    Pkg.develop(path=muscadeDir)
    Pkg.instantiate()
end

# run all Muscade tests, 
workingDir = joinpath(muscadeDir,"test")
devMuscadeIn(workingDir,muscadeDir)
include(joinpath(workingDir,"runtests.jl"));

# run all examples
workingDir = joinpath(muscadeDir,"examples")
devMuscadeIn(workingDir,muscadeDir)
include(joinpath(workingDir,"StaticBeamAnalysis.jl"))
include(joinpath(workingDir,"DynamicBeamAnalysis.jl"))
include(joinpath(workingDir,"DecayAnalysis.jl"))
include(joinpath(workingDir,"DryFriction.jl"))

# generate documentation
workingDir = joinpath(muscadeDir,"docs")
devMuscadeIn(workingDir,muscadeDir)
include(joinpath(workingDir,"make.jl"))
