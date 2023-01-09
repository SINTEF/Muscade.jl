module Muscade
    using  EspyInsideFunction
    export EspyInsideFunction,@request,makekey,forloop,scalar,@espy,@espydbg
 
    using  Printf,SparseArrays,StaticArrays,LinearAlgebra
    # using  Base.Threads
    # import Base.Threads.@spawn, Base.Threads.nthreads

    include("core/Dialect.jl")
    export ℝ,ℤ,𝕣,𝕫,𝔹,𝕓
    export ℝ1,ℤ1,𝕣1,𝕫1,𝔹1,𝕓1
    export ℝ2,ℤ2,𝕣2,𝕫2,𝔹2,𝕓2
    export ℝ11,ℤ11,𝕣11,𝕫11,𝔹11,𝕓11
    export toggle,@once

    include("core/Adiff.jl")
    export  ∂ℝ #\partial \bbR
    export  variate,δ # \delta
    export  value,VALUE,∂ # \partial, \nabla
    export  constants,precedence,npartial,norm
    export  ≗ #\circeq

    include("core/Dots.jl")
    export dots,∘₀,∘₁,∘₂,⊗

    include("core/Exceptions.jl")
    export muscadeerror

    include("core/ElementAPI.jl")
    export AbstractElement
    export lagrangian,residual,espyable,draw,request2draw 
    export ∂0,∂1,∂2
    export doflist

    include("core/ModelDescription.jl")
    export Model,addnode!,addelement!,coord,Node,DofID,NodID,EleID,setscale! 
    export getndof,getnele,getdofID,describe

    include("core/BasicElements.jl")
    export XdofCost,UdofCost,AdofCost,DofLoad,DofHold,Spring

    include("core/Assemble.jl")
 
    include("core/Solve.jl")
    export solve

    include("core/StaticX.jl")
    export staticX
    
    include("core/StaticXUA.jl")
    export staticXUA

    include("core/Output.jl")
    export getdof,getresult

    include("core/SelfDraw.jl")
    export draw
    
    module Unit
        include("core/Unit.jl")
        export unit,←,→,convert
    end
    module ElTest
        include("core/ElTest.jl")
        export test_static_element#,nodes_for_element_test
    end
end
