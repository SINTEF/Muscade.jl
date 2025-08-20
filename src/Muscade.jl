module Muscade
    using Printf
    using LinearAlgebra
    using SparseArrays
    using StaticArrays
    using SpecialFunctions
    using KrylovKit: KrylovKit,eigsolve
    using MacroTools
    using MacroTools: postwalk,gensym_ids,rmlines,unblock 
    using Base.Cartesian
    using GLMakie

    include("Dialect.jl")
    export ℝ,ℤ,𝕣,𝕫,𝔹,𝕓,ℂ
    export ℝ1,ℤ1,𝕣1,𝕫1,𝔹1,𝕓1
    export ℝ2,ℤ2,𝕣2,𝕫2,𝔹2,𝕓2
    export ℝ11,ℤ11,𝕣11,𝕫11,𝔹11,𝕓11
    export toggle,default,@once,mod_onebased

    include("Adiff.jl")
    export  ∂ℝ #\partial \bbR
    export  variate,δ,directional # \delta
    export  value,VALUE,∂,value_∂ # \partial, \nabla
    export  constants,precedence,npartial,norm

    include("Taylor.jl")
    export  motion,motion⁻¹,revariate,compose,fast,apply,justinvoke,composevalue,composeJacobian 

    include("Functors.jl")
    export QuadraticFunction,FunctionFromVector 

    include("Dots.jl")
    export dots,∘₀,∘₁,∘₂,⊗

    include("Espy.jl") 
    export @request, mergerequest
    export @espy,@espydbg

    include("Exceptions.jl")
    export muscadeerror

    include("ModelDescription.jl")
    export AbstractElement
    export Model,addnode!,addelement!,setscale!,initialize!
    export Node
    export getndof

    include("ElementAPI.jl")
    export coord,∂0,∂1,∂2,getsomedofs
    export noFB

    include("ElementTestTools.jl")
    export diffed_residual,diffed_lagrangian,print_element_array, @typeof

    include("BasicElements.jl")
    export off,equal,positive
    export DofCost,SingleDofCost,SingleUdof,ElementCost,Acost,SingleAcost
    export DofConstraint,Hold,ElementConstraint
    export QuickFix,DofLoad

    include("Assemble.jl")
    export Assembly

    include("Solve.jl")
    export solve

    include("SparseTools.jl")
    export prepare,cat!,addin!,zero!,getblock
    
    include("FiniteDifferences.jl")

    include("SweepX.jl")
    export SweepX

    include("DirectXUA.jl")
    export DirectXUA

    include("EigX.jl")
    export EigX,increment

    include("EigXU.jl")
    export EigXU,GUI

    include("FreqXU.jl")
    export FreqXU

    include("Diagnostic.jl")
    export studyscale,studysingular,describe

    include("Output.jl")
    export setdof!,getdof,getresult,findlastassigned,eletyp

    include("SelfDraw.jl")
    export draw!,request2draw

    include("Unit.jl")
    export ←,→

    include("FFT.jl")
    #export getδf,getδt(n3,δf3′),𝔉𝕣(g.(t3),δt3),𝔉𝕣⁻¹(X3′′,δf3)

    include("Eigenmodes.jl")
end
