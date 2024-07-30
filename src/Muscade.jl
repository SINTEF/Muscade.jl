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

    include("Multiplex.jl")
    
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
    export noFB

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

    include("SweepX.jl")
    export SweepX

    include("StaticXUA.jl")
    export StaticXUA

    include("Diagnostic.jl")
    export studyscale,studysingular,describe

    include("Output.jl")
    export setdof!,getdof,getresult,findlastassigned

    include("SelfDraw.jl")
    export draw,request2draw

    export Unit
    module Unit  # using Muscade.Unit
        include("Unit.jl")
        export unit,←,→,convert
    end

    export ElementTestTools
    module ElementTestTools # using Muscade.ElementTestTools
        using Muscade
        include("ElementTestTools.jl")
        export test_static_element,gradient
    end

    export Elements
    module Elements  # using Muscade.Elements: EulerBeam3D
        using Muscade
        include("Elements/BeamElement.jl")
        include("Elements/DryFriction.jl")
    end
end
