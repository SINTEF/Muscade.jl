module Muscade
    using  EspyInsideFunction
    export EspyInsideFunction,@request,makekey,forloop,scalar,@espy,@espydbg
 
    using  Printf,SparseArrays,StaticArrays,LinearAlgebra
    # using  Base.Threads
    # import Base.Threads.@spawn, Base.Threads.nthreads


    ## TEMPORARY STUFF
    using StaticArrays    
    struct Node
        coords :: SVector{3,Float64}
    end
    coords(n)= SMatrix{1,3}(n[i].coords[j] for i∈eachindex(n), j∈1:3)
    export Node,coords    
    ##

    include("core/Dialect.jl")
    export ℝ,ℤ,𝕣,𝕫,𝔹,𝕓
    export ℝ1,ℤ1,𝕣1,𝕫1,𝔹1,𝕓1
    export ℝ2,ℤ2,𝕣2,𝕫2,𝔹2,𝕓2
    export ℝ11,ℤ11,𝕣11,𝕫11,𝔹11,𝕓11
    export toggle

    include("core/Dots.jl")
    export dots,∘₀ ,∘₁,∘₂,⊗

    include("core/Exceptions.jl")
    export muscadeerror

    include("core/ElementAPI.jl")
    export Element
    export initχ,lagrangian,residual,espyable,draw,request2draw # element API
    export ∂0,∂1,∂2
    export Xdofid,Udofid,Adofid,dofid,neldof

    module Unit
        include("Core/Unit.jl")
        export unit,←,→,convert
    end
    module ElTest
        include("Core/ElTest.jl")
        export testStaticElement,nodesforelementtest
    end
end
