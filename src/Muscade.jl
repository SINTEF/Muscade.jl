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
    export default
    public toggle,mod_onebased

    include("OffsetVector.jl")

    include("Adiff.jl")
    export  ∂ℝ #\partial\BbbR
    export  value,VALUE,∂,value_∂ # \partial, \nabla
    public  precedence,npartial #,norm

    include("Taylor.jl")
    export motion,motion⁻¹,variate,variate0,revariate,chainrule,apply,AllElements
    export chainrule_Jacobian 
    public Taylor,McLaurin
    
    include("Functors.jl")
    export @functor, Functor
    public FunctionFromVector

    include("Dots.jl")
    export ∘₀,∘₁,∘₂,⊗
    public dots

    include("Espy.jl") 
    export @request,@espy,mergerequest
    public @espydbg

    include("Exceptions.jl")
    export muscadeerror

    include("ModelDescription.jl")
    export AbstractElement
    export Model,addnode!,addelement!,setscale!,initialize!
    export Node
    export getndof  
    public get # TODO does this get doc'ed?

    include("ElementAPI.jl")
    export coord,∂0,∂1,∂2,getsomedofs
    export noFB
    public no_second_order
    public allocate_drawing, update_drawing, display_drawing!
    public residual,lagrangian,doflist

    include("BasicElements.jl")
    public off,equal,positive
    export DofCost,SingleDofCost,SingleUdof,ElementCost,Acost,SingleAcost
    export DofConstraint,Hold,ElementConstraint
    export DofLoad
    public QuickFix

    include("Assemble.jl")

    include("Solve.jl")
    export solve

    include("SparseTools.jl")
    public prepare,addin!
    
    include("FiniteDifferences.jl")

    include("SweepX.jl")
    export SweepX

    include("SweepXA.jl")
    export SweepXA

    include("DirectXUA.jl")
    export DirectXUA

    include("EigX.jl")
    export EigX,increment

    include("EigXU.jl")
    export EigXU

    include("FreqXU.jl")
    export FreqXU

    include("Diagnostic.jl")
    export describe  
    public study_scale,study_singular,plot_matrix_sparsity
    public plot_block_matrix_sparsity,print_nz,Monitor,@typeof,print_element_array,diffed_lagrangian,diffed_residual
    public SpyAxis

    include("Output.jl")
    export setdof!,getdof,getresult,findlastassigned

    include("SelfDraw.jl")
    export draw!,request2draw
    public GUI

    include("FFT.jl")
    public getδω,getδt,𝔉,𝔉⁻¹

    include("Eigenmodes.jl")

    public Toolbox
    module Toolbox
        include("../toolbox/Basics.jl")
        export clutch        

        include("../toolbox/Rotations.jl")
        public scac, sinc1, sinc1′,sinc1″, sinc1‴, sinc1⁗
        public spin,spin⁻¹,Rodrigues, Rodrigues⁻¹, adjust , intrinsicrotationrates
        include("../toolbox/BarElement.jl")
        export Bar3D, AxisymmetricBarCrossSection
        
        include("../toolbox/BeamElement.jl")
        export EulerBeam3D, BeamCrossSection
        
        include("../toolbox/StrainGaugeOnBeamElement.jl")
        export EulerBeam3DwithStrainGauge
        
        include("../toolbox/PositionElement.jl")
        export Position3D  
        
        include("../toolbox/SoilContact.jl")
        export SoilContact

        include("../toolbox/RigidBodyKinematics.jl")
        export ExcentricRigidConnection

        include("../toolbox/MeshLine.jl")
        export MeshLine!

        include("../toolbox/Unit.jl") 
        export ←,→
    end

end
