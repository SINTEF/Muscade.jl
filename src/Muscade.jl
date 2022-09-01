module Muscade
    using  EspyInsideFunction
    export EspyInsideFunction,@request,makekey,forloop,scalar,@espy,@espydbg
 
    using  Printf,SparseArrays,StaticArrays,LinearAlgebra
    # using  Base.Threads
    # import Base.Threads.@spawn, Base.Threads.nthreads


    using StaticArrays    
    struct Node
        coords :: SVector{3,Float64}
    end
    coords(n)= SMatrix{1,3}(n[i].coords[j] for i∈eachindex(n), j∈1:3)
    export Node,coords    


    # export muscadeerror
    # include(core/Exceptions.jl)
    include("core/Dialect.jl")
    include("core/Dots.jl")
    include("Core/Unit.jl")


    include("core/ElementAPI.jl")
    export AbstractElement,nostate
    export initstate,lagrangian,espyable,draw,request2draw # element API
    export ∂0,∂1,∂2
    export Xdofid,Udofid,Adofid,dofid,neldof
    # export dofid,neldof

    include("Core/ElementTestBench.jl")
    export testStaticElement,nodesforelementtest
    

    module Tools    
        module Dialect
            include("core/Dialect.jl")
            export ℝ,ℤ,𝕣,𝕫
            export ℝ1,ℤ1,𝕣1,𝕫1
            export ℝ2,ℤ2,𝕣2,𝕫2
            export ℝ11,ℤ11,𝕣11,𝕫11
            export toggle
        end
        module Dots
           include("core/Dots.jl")
            export dots,∘₀ ,∘₁,∘₂,⊗
        end
        module Unit
            include("Core/Unit.jl")
            export unit,←,→,convert
        end
#        module ElementTestBench
#        end
    end 
end
