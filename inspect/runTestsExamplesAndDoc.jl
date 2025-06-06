# Runs tests, examples and documentation generation 
# from a clean environment, using the version of Muscade
# located in "whereIsMuscade". Useful to run prior to 
# a merge.

whereIsMuscadeDev = "C:\\Users\\thsa\\code\\Muscade.jl\\"
# whereIsMuscadeDev = "/home/thomash/Muscade.jl/"

using Pkg

function devMuscadeIn(workingDir,whereIsMuscadeDev)
	cd(workingDir)
	rm("Manifest.toml")    
	Pkg.activate(".")
    try Pkg.rm("Muscade") catch err end
    Pkg.develop(path=whereIsMuscadeDev)
    Pkg.instantiate()
end

# run all Muscade tests, 
workingDir = joinpath(whereIsMuscadeDev,"test")
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(joinpath(workingDir,"runtests.jl"));

# run all examples
workingDir = joinpath(whereIsMuscadeDev,"examples")
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(joinpath(workingDir,"StaticBeamAnalysis.jl"))
include(joinpath(workingDir,"DynamicBeamAnalysis.jl"))
include(joinpath(workingDir,"DecayAnalysis.jl"))
include(joinpath(workingDir,"DryFriction.jl"))

# generate documentation
workingDir = joinpath(whereIsMuscadeDev,"docs")
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(joinpath(workingDir,"make.jl"))
