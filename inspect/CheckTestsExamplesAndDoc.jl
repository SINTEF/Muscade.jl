whereIsMuscadeDev = "C:\\Users\\thsa\\code\\Muscade.jl\\"

using Pkg

function devMuscadeIn(workingDir,whereIsMuscadeDev)
    cd(workingDir)
    Pkg.activate(".")
    try Pkg.rm("Muscade") catch err end
    Pkg.develop(path=whereIsMuscadeDev)
    Pkg.instantiate()
end

# run all Muscade tests, 
workingDir = whereIsMuscadeDev*"test\\"
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(workingDir*"runtests.jl");

# run all examples
workingDir = whereIsMuscadeDev*"examples\\"
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(workingDir*"StaticBeamAnalysis.jl")
include(workingDir*"DynamicBeamAnalysis.jl")
include(workingDir*"DecayAnalysis.jl")
include(workingDir*"DryFriction.jl")

# generate documentation
workingDir = whereIsMuscadeDev*"docs\\"
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(workingDir*"make.jl")
