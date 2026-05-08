function makeXUnorm!(vec,dofgr,σ,def=(1.,∞,∞))
    # in.class.field[ider] = σ
    getider(a::NamedTuple,ider,defᵢ) = map(aⱼ->getider(aⱼ,ider,defᵢ),a)
    getider(a::NTuple    ,ider,defᵢ) = a[ider]
    getider(a::𝕣         ,ider,defᵢ) = ider==1 ? a : defᵢ
    for ider ∈ eachindex(vec)
        vec[ider] .= def[ider]
        makevecfromfields!(vec[ider],dofgr,getider(σ,ider,def[ider]))
        vec[ider] .= vec[ider].^(-2)
    end 
end

function make_λxu_sparsepattern(out) 
    L2(α,β) = out.L2[α,β][1,1]
    α       = [2,3,1,2,3,1,2,3]  #   [0 . .]
    β       = [1,1,2,2,2,3,3,3]  #   [. . .]
    return sparse(α,β,L2.(α,β))  # = [. . .]
end
function assemblebigmat!(L2::Vector{Sparse𝕣2},L2bigasm::SparseMatrixCSC,asm,model,dis,out::AssemblyDirect,dbg) 
    # does not call assemble!: solve has previously called assemble! to prepare bigasm, so out.L2 is already set,
    for L2ᵢ∈L2
        zero!(L2ᵢ)
    end
    for     α ∈ λxu 
        for β ∈ λxu
            Lαβ = out.L2[α,β]
            for     αder = 1:size(Lαβ,1)
                for βder = 1:size(Lαβ,2)
                    ider =  αder+βder-1   
                    if isassigned(Lαβ,αder,βder)
                        sgn  = isodd(αder) ? +1 : -1 # TODO Antisymmetry for odd derivatives? conjugation? Check theory.  See also DirectXUA
                        addin!(L2bigasm,L2[ider],Lαβ[αder,βder],α,β,sgn) 
                    end
                end
            end
        end
    end
end
function assemblebigvec!(L1,L1bigasm::𝕫1,asm,model,dis,out::AssemblyDirect,state,Δt,dbg) 
    zero!.(L1)
    assemble!{:vectors}(out,asm,dis,model,state,Δt,(dbg...,asm=:assemblebigvec!)) # first assemble model vectors
    for β ∈ λxu                                                                # then collate them into
        Lβ = out.L1[β]
        for βder = 1:size(Lβ,1)
            addin!(L1bigasm,L1[βder],Lβ[βder],β,idmult) 
        end
    end
end

struct EigXUincrement{Tω}
    nmod  :: 𝕫
    dofgr :: DofGroup
    ω     ::Tω           # [iω] range
    ncv   :: 𝕫1          # [iω]
    λ     :: 𝕣11         # [iω][imod]
    nor   :: 𝕣11         # [iω][imod]
    ΔΛXU  :: Vector{𝕣11} # [iω][imod][idof]
end



"""
	EigXU{OX,OU}

Study the combinations of load and response that are least detected by sensor systems.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
eiginc          = solve(EigXU{OX,OU};Δω, p, nmod,initialstate)
```

The solver linearises the problem (computes the Hessian of the Lagrangian) at `initialstate` and solves 
the ΛXU-eigenvalue problem at frequencies ωᵢ = Δω*i with i∈{0,...,2ᵖ-1}.


# Parameters
- `OX`                0 for static analysis
                      1 for first OX problems in time (viscosity, friction, measurement of velocity)
                      2 for second OX problems in time (inertia, measurement of acceleration) 
- `OU`                0 for white noise prior to the unknown load process
                      2 otherwise

# Named arguments
- `dbg=(;)`           a named tuple to trace the call tree (for debugging).
- `verbose=true`      set to false to suppress printed output (for testing).
- `initialstate`      a `State`, at which the problem is linearized.
- `nmod`              the number of eigen-modes to identusy
- `Δω`                frequency step
- `p`                 `2^p` steps will be analysed.      
- `droptol=1e-10`     set to zero terms in the incremental matrices that are smaller than `droptol` in absolute value.                      

# Output
- an object of type `EigXUincrement` for use with [`increment`](@ref) to create a snapshot of the
  oscillating system.

See also: [`increment`](@ref),[`EigXU`](@ref), [`solve`](@ref), [`initialize!`](@ref), [`study_singular`](@ref), [`SweepX`](@ref), [`DirectXUA`](@ref)
"""
struct EigXU{OX,OU} <: AbstractSolver end 

function solve(::Type{EigXU{OX,OU}},pstate,verbose::𝕓,dbg;
    Δω::𝕣, p::𝕫, 
    initialstate::State,
    droptol::𝕣=1e-10,
    nmod::𝕫=5,
    σₓᵤ,
    kwargs...) where{OX,OU}

    #  Mostly constants
    local LU
    model,dis             = initialstate.model, initialstate.dis
    nω                    = 2^p
    IA                    = 0
    # State storage
    S                     = State{1,OX+1,OU+1,Nothing}
    pstate[] = state      = Vector{S}(undef,nω)                                                                           
    state₀                = State{1,OX+1,OU+1}(copy(initialstate))   

    verbose && @printf("    Preparing assembler\n")
    wanted                = Wanted{1,OX+1,OU+1,IA}(:all,:all) # TODO refine
    out,asm,dofgr         = prepare(AssemblyDirect,model,dis,wanted)   # model assembler for all arrays   

    verbose && @printf("    Computing matrices\n")
    assemble!{:matrices}(out,asm,dis,model,state₀,idmult,(dbg...,solver=:EigXU,phase=:matrices))            # assemble all model matrices - in class-blocks
    pattern               = make_λxu_sparsepattern(out)
    L2                    = Vector{Sparse𝕣2}(undef,5)
    L2[1],L2bigasm,L1bigasm,Ldis  = prepare(pattern)  
    λxu_dofgr             = allΛXUdofs(model,dis)                                        # NB same ordering of dofs in rhs as implied by pattern                                          
    for ider              = 2:5
        L2[ider]          = copy(L2[1])
    end    
    assemblebigmat!(L2,L2bigasm,asm,model,dis,out,(dbg...,solver=:EigXU))              # assemble all complete model matrices into L2
    nXdof,nUdof           = getndof(model,(:X,:U))
    ixu                   = (nXdof+1):(2nXdof+nUdof)
    B                     = sparse(ixu,ixu,𝕣1(undef,nXdof+nUdof)) # ndof×ndof
    N                     = [𝕣1(undef,nXdof+nUdof) for ider = 1:3]
    xu_dofgr              = allXUdofs(model,dis)                                        # NB same ordering of dofs in rhs as implied by pattern                                          
    makeXUnorm!(N,xu_dofgr,σₓᵤ)   

    verbose && @printf("    Improving sparsity ")    
    keep                  = sparser!(L2,droptol)
    verbose && @printf("from %i to %i nz terms\n",length(keep),sum(keep))    

    verbose && @printf("    Solving XU-eigenproblem for all ω\n")
    L2₁                   = L2[1]
    ndof                  = 2nXdof+nUdof
    A                     = Sparse𝕣2(ndof,ndof,L2₁.colptr,L2₁.rowval,𝕣1(undef,length(L2₁.nzval))) # but could be complex
    ΔΛXU                  = Vector{𝕣11}(undef,nω) # ΔΛXU[iω][imod][idof]
    λ                     = 𝕣11(undef,nω)         # λ[   iω][imod]
    nor                   = 𝕣11(undef,nω)         # nor[ iω][imod] 
    ncv                   = 𝕫1(undef,nω)          # ncv[ iω]
    wrk                   = zeros(ndof)           # wrk[ndof]

    ω                     = range(start=0.,step=Δω,length=nω) 
    for (iω,ωᵢ)           = enumerate(ω)
        B.nzval          .= N[1]        + ωᵢ^2*N[2]        + ωᵢ^4*N[3]     
        A.nzval          .= L2[1].nzval + ωᵢ^2*L2[3].nzval + ωᵢ^4*L2[5].nzval     # complex if exponents 1 and 3 included
        try 
            if iω==1 LU   = lu(A) 
            else     lu!(LU ,A)
            end 
            λ⁻¹, ΔΛXU[iω], ncv[iω] = geneig{:symmetric}(A,B,nmod;normalize=false,kwargs...)
            nor[iω]                = 𝕣1(undef,ncv[iω])
            λ[iω]                  = 1 ./λ⁻¹
            for imod               = 1:ncv[iω]
                Δ                  = ΔΛXU[iω][imod]
                wrk[ixu]          .= view(Δ,ixu)                   # this copy can be optimised by viewing the classes in Δ, operating on out.L2[α,β][αder,βder], and combining over derivatives.  Is it worth the effort?   
                Anorm              = √(ℜ(dot(wrk,A,wrk))/2)  # ΔΛXU is real, A is complex Hermitian, so square norm is real: (imag part is zero to machine precision)
                if iω>1  &&  imod≤nmod  &&  sum(Δ[idof]*ΔΛXU[iω-1][imod][idof] for idof∈λxu_dofgr.jX)<0
                    Anorm          = -Anorm
                end
                Δ                .*= 2.575829303549/Anorm          # corresponds to a probability of exceedance of 0.01                        
                nor[iω][imod]      = √(ℜ(dot(Δ,B,Δ))/2) 
            end
        catch 
            muscadewarning(@sprintf("Factorization of matrix A failed for ω=%f",ωᵢ));
        end
    end    
    any(ncv.<nmod) && verbose && muscadewarning("Some eigensolutions did not converge",4)
    pstate[] = EigXUincrement(nmod,allΛXUdofs(model,dis),ω,ncv,λ,nor,ΔΛXU)
    verbose && @printf("\n")
    return
end
"""
    state = increment{OX}(initialstate,eiginc,iω,imod,A)

Starting from `initalstate` for which an `EigX` analysis has been carried out, and using the output
`eiginc` of that analysis, construct new `State`s representing the instantaneous state of the 
vibrating structure
    
# Input
- `OX` the number of time derivatives to be computed.  `increment(initialstate,eiginc,imod,A)` defaults to `OX=2`
- `initialstate` the same initial `State` provided to `EigXU` to compute `eiginc`
- `eiginc` obtained from `EigXU`
- `iω`, the number of the frequency to consider. `ω=iω*Δω` where `Δω` is an input to [`EigXU`](@ref). 
- `imod`, an `AbstractVector` of integer mode numbers
- `A`, an `AbstractVector` of same length as `imod`, containing real or complex 
  amplitudes associated to the modes

# Output
- `state` a snapshot of the vibrating system

See also: [`EigXU`](@ref)
"""
function increment{OX}(initialstate,eiginc::EigXUincrement,iω::𝕫,imod::AbstractVector{𝕫},amplitude::AbstractVector) where{OX} 
    state       = State{1,OX+1,1}(copy(initialstate)) 
    ω, ΔΛXU     = eiginc.ω[iω], eiginc.ΔΛXU[iω]
    maximum(imod)≤length(eiginc.λ) || muscadeerror(@sprintf("eiginc only has %n modes for iω=%i.",length(ω),iω))
    for (i,imodᵢ)∈enumerate(imod)  
        for iOX = 0:OX
            increment!(state,iOX+1,ℜ.(ω^iOX*amplitude[i]*ΔΛXU[imodᵢ]),eiginc.dofgr)
        end
    end
    return state
end

# scales Λ,X and U differently (so is no longer a solution state) and only updates 0th derivatives. For graphical outputs
function visualincrement(initialstate,eiginc::EigXUincrement,iω::𝕫,imod::𝕫;Λscale::𝕣=1.,Xscale::𝕣=1.,Uscale::𝕣=1.) 
    imod≤length(eiginc.λ) || muscadeerror(@sprintf("eiginc only has %n modes for iω=%i.",length(ω),iω))
    model,dis             = initialstate.model, initialstate.dis
    gr, ΔΛXU              = eiginc.dofgr, eiginc.ΔΛXU[iω][imod]
    state                 = State{1,1,1}(copy(initialstate)) 
    for i ∈ eachindex(gr.iΛ); state.Λ[1][gr.iΛ[i]] += ΔΛXU[gr.jΛ[i]] * gr.scaleΛ[i] *Λscale end
    for i ∈ eachindex(gr.iX); state.X[1][gr.iX[i]] += ΔΛXU[gr.jX[i]] * gr.scaleX[i] *Xscale end
    for i ∈ eachindex(gr.iU); state.U[1][gr.iU[i]] += ΔΛXU[gr.jU[i]] * gr.scaleU[i] *Uscale end
    return state
end
"""

    Muscade.GUI(eiginc,initialstate;[draw_shadow=true],[shadow=...],[model=...])

Taking the output `eiginc` obtained from an `EigXU`, and the state `initstate` provided to `EigXU`, provide
a GUI to explore the results.

Optional keyword arguements are
- `draw_shadow` whether to superimpose a drawing of `initstate`
- `shadow` a `NamedTuple` with any arguments to be passed to `draw!` `initstate`
- `model` a `NamedTuple` with any arguments to be passed to `draw!` the `EigXU` modes.

See also [`EigXU`](@ref)
"""
function GUI(initialstate,eiginc::EigXUincrement,;kwargs...)
    args = default(kwargs, (draw_shadow=true, shadow=(;),model=(;)))

    ## Organize the window

    fig             = Figure(size = (1500,900))
    GLMakie.activate!( title = "Muscade.jl" )
    display(fig) # open interactive window (gets closed down by "save")
    panelFreqs      = fig[1,1]        
    panelNorm       = panelFreqs[1,1] 
    axisNorm        = Axis(panelNorm,xlabel="ω [rad/s]",ylabel="magnitude of error",yscale=log10)
    panelSlide      = panelFreqs[2,1] 
    panelModel      = fig[1,2:3]        
    Box(panelModel, cornerradius = 20,z=1., color = :transparent)
    axisModel       = Axis3(panelModel,title="EigXU mode shape",aspect=:data,viewmode=:free,perspectiveness=.5,clip=false)

    ## sliders

    ω0 = eiginc.ω[div(length(eiginc.ω),3)]
    sg = SliderGrid(panelSlide,
                    (label="ω"      , range = eiginc.ω        , startvalue = ω0,snap=true,update_while_dragging=true,format = "{:.1f} rad/s"),
                    (label="mode"   , range = 1:eiginc.nmod   , startvalue = 1 ,snap=true,update_while_dragging=true                        ),
                    (label="X scale", range = -5:0.01:5       , startvalue = 0 ,snap=true,update_while_dragging=true,format = "10^{:.1f}"   ),
                    (label="U scale", range = -5:0.01:5       , startvalue = 0 ,snap=true,update_while_dragging=true,format = "10^{:.1f}"   ))
    obs = (ω      = sg.sliders[1].value,
           imode  = sg.sliders[2].value,
           Xscale = sg.sliders[3].value,
           Uscale = sg.sliders[4].value)
    iωs  = map(obs.ω) do ω
        round(Int64,ω/step(eiginc.ω))+1
    end       
    nors = map(obs.imode,iωs) do imode,iω 
        eiginc.nor[iω][imode]
    end

    ## Model

    args.draw_shadow && draw!(axisModel,initialstate;args.shadow...); # draw initial state once to keep on screen
    
    graphic = draw!(axisModel,initialstate;args.model...);                             # and twice to start the pump
    _ = map(iωs,obs.imode,obs.Xscale,obs.Uscale) do iω,imod,Xscale,Uscale                                    # Then observe the sliders
        state = Muscade.visualincrement(initialstate,eiginc,iω,imod;Xscale=exp10(Xscale),Uscale=exp10(Uscale))
        draw!(graphic,state;args.model...);
    end

    ## norm spectre

    nω  = length(eiginc.ω)
    nor = 𝕣1(undef,nω)
    #λ   = 𝕣1(undef,nω)
    for imod = 1:maximum(eiginc.ncv)
        for iω= 1:nω
            if imod≤eiginc.ncv[iω]
                nor[iω] = eiginc.nor[iω][imod]
                #λ[  iω] = eiginc.λ[  iω][imod]
            else
                nor[iω] = NaN
                #λ[  iω] = NaN
            end
        end
        scatter!(axisNorm,eiginc.ω,nor,markersize=1,color=:black)
    end
    scatter!(axisNorm,obs.ω,nors,color=:red,markersize=10)

end
