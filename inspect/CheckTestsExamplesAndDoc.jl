using Pkg

function devMuscadeIn(workingDir,whereIsMuscadeDev)
    cd(workingDir)
    # versioninfo();
    # @show pwd();
    try Pkg.rm("Muscade") catch err end
    Pkg.develop(path=whereIsMuscadeDev)
    Pkg.activate(".")
    Pkg.instantiate()
    # Pkg.status()
end

whereIsMuscadeDev = "C:\\Users\\thsa\\code\\Muscade.jl\\"

# run all Muscade tests, 
workingDir = whereIsMuscadeDev*"test\\"
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(workingDir*"runtests.jl");

# run all examples
workingDir = whereIsMuscadeDev*"examples\\"
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(workingDir*"StaticBeamAnalysis.jl")
include(workingDir*"DecayAnalysis.jl")
include(workingDir*"DryFriction.jl")

# generate documentation
workingDir = whereIsMuscadeDev*"docs\\"
devMuscadeIn(workingDir,whereIsMuscadeDev)
include(workingDir*"make.jl")
