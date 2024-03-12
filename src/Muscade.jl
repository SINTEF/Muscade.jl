module Muscade
    using  Printf,SparseArrays,StaticArrays,LinearAlgebra

    include("Dialect.jl")
    export ℝ,ℤ,𝕣,𝕫,𝔹,𝕓
    export ℝ1,ℤ1,𝕣1,𝕫1,𝔹1,𝕓1
    export ℝ2,ℤ2,𝕣2,𝕫2,𝔹2,𝕓2
    export ℝ11,ℤ11,𝕣11,𝕫11,𝔹11,𝕓11
    export toggle,@once,default

    include("Adiff.jl")
    export  ∂ℝ #\partial \bbR
    export  variate,δ,directional # \delta
    export  value,VALUE,∂,value_∂ # \partial, \nabla
    export  constants,precedence,npartial,norm
    export  ≗ #\circeq

    include("Dots.jl")
    export dots,∘₀,∘₁,∘₂,⊗

    include("Espy.jl") 
    export @request
    export @espy,@espydbg

    include("Exceptions.jl")
    export muscadeerror


    include("ModelDescription.jl")
    export AbstractElement
    export Model,addnode!,addelement!,setscale!,initialize!
    export Node
    export getndof

    include("ElementAPI.jl")
    export lagrangian,residual,espyable
    export coord,∂0,∂1,∂2,getsomedofs
    export doflist
    export noχ,noFB

    include("BasicElements.jl")
    export off,equal,positive
    export DofCost,SingleDofCost,ElementCost
    export DofConstraint,Hold,ElementConstraint
    export QuickFix,DofLoad

    include("Assemble.jl")
    export Assembly
    
    include("Solve.jl")
    export solve

    include("BlockSparse.jl")
    export prepare,cat!,addin!,zero!,getblock

    include("StaticX.jl")
    export StaticX
    
    include("StaticXUA.jl")
    export StaticXUA

    include("Diagnostic.jl")
    export studyscale,studysingular,describe

    include("Output.jl")
    export setdof!,getdof,getresult,findlastassigned

    include("SelfDraw.jl")
    export draw,request2draw

    
    module Unit
        include("Unit.jl")
        export unit,←,→,convert
    end
    module ElTest
        using Muscade
        include("ElTest.jl")
        export test_static_element,gradient
    end
end
