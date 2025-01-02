test    = @__DIR__
muscade = normpath(joinpath(test,".."))
#docs    = normpath(joinpath(test,"../docs"))
using Pkg
Pkg.activate(test)
#using Muscade # seems necessary for doc test to work on a cold start
include(normpath(joinpath(test,"runtests.jl")))
Pkg.activate(muscade) 
