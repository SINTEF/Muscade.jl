"""
	eiginc = solve(EigX{ℝ};state=initialstate,nmod)
	eiginc = solve(EigX{ℂ};state=initialstate,nmod)

Given an initial (typicaly static) state `initialstate`, computes the lowest `nmod` eigenmodes of a system. 
`EigX{ℝ}` computes real eigenmodes not accounting for damping. `EigX{ℝ}` computes complex eigenmodes accounting for damping.
The data structure `eiginc` can be passed to `increment` to obtain dynamic states superimposing mode shapes.

# Input
- `initialstate` - a `State`
- `nmod=5` - the number of eigenmodes to identify
- `droptol=1e-9` - in the stiffness and mass matrix, the magnitude of a term relative to the largest term in the matrix
        under which the term is set to zero.
- Further named arguments: see the optional keyword arguments to `geneig`.         

# Output
- an object of type `EigXℝincrement` or `EigXℂincrement`, for use with [`increment`](@ref) to create a snapshot of the
  oscillating system.

See also: [`solve`](@ref), [`initialize!`](@ref), [`increment`](@ref), [`geneig`](@ref)
"""
struct EigX{T} <: AbstractSolver end

struct EigXℝincrement
    dofgr :: DofGroup
    ω     :: 𝕣1
    Δx    :: 𝕣11
end
function solve(::Type{EigX{ℝ}},pstate,verbose,dbg; 
                   state::State, nmod::𝕫=5,droptol::𝕣=1e-9,kwargs...) 
    OX,OU,IA         = 2,0,0
    model,dis        = state.model,state.dis

    verbose && @printf("\n    Assembling\n")
    wantK,wantM      = (ind.Λ,ind.X,1,1),(ind.Λ,ind.X,1,3)
    wanted           = Wanted{1,OX+1,OU+1,IA}(:all,(wantK,wantM))
    out,asm,dofgr    = prepare(AssemblyDirect,model,dis,wanted)  
    nXdof            = getndof.(dofgr)[ind.X]
    state₀           = State{1,OX+1,OU+1}(copy(state)) 
    assemble!{:matrices}(out,asm,dis,model,state₀,idmult,(dbg...,solver=:EigXℝ))
    K                = out.L2[ind.Λ,ind.X][1,1]
    M                = out.L2[ind.Λ,ind.X][1,3]
    sparser!(K,droptol)
    sparser!(M,droptol)

    verbose && @printf("    Solving Eigenvalues\n\n")
    ω², Δx, ncv = geneig{:symmetric}(K,M,nmod;kwargs...)

    ncv≥nmod||muscadeerror(dbg,@sprintf("eigensolver only converged for %i out of %i modes",ncv,nmod))
    pstate[] = EigXℝincrement(dofgr[ind.X],sqrt.(ω²),Δx)
    verbose && @printf("\n")
    return 
end
"""
    state = increment{OX}(initialstate,eiginc,imod,A)

Starting from `initalstate` for which an `EigX` analysis has been carried out, and using the output
`eiginc` of that analysis, construct new `State`s representing the instantaneous state of the 
vibrating structure
    
# Input
- `OX` the number of time derivatives to be computed.  `increment(initialstate,eiginc,imod,A)` defaults to `OX=2`
- `initialstate` the same initial `State` provided to `EigX` to compute `eiginc`
- `eiginc` obtained from `EigX`
- `imod`, an `AbstractVector` of integer mode numbers
- `A`, an `AbstractVector` of same length as `imod`, containing real or complex 
  amplitudes associated to the modes

# Output
- `state` a snapshot of the vibrating system

See also: [`EigX`](@ref)
"""
struct increment{OX} end
increment(args...) = increment{2}(args...)
function increment{OX}(initialstate,eiginc::EigXℝincrement,imod::AbstractVector{𝕫},A::AbstractVector) where{OX}
    length(imod)==length(A)|| muscadeerror("imod and A must be of same length.")
    state            = State{1,OX+1,1}(copy(initialstate)) 
    maximum(imod)≤length(eiginc.ω) || muscadeerror(@sprintf("eiginc only has %n modes.",length(ω)))
    for i∈eachindex(imod)  
        ωᵢ,Δxᵢ = eiginc.ω[imod[i]],eiginc.Δx[imod[i]]
        for iOX     = 0:OX
            increment!(state,iOX+1,ℜ.((𝑖*ωᵢ)^iOX*A[i]*Δxᵢ),eiginc.dofgr)
        end
    end
    return state
end

struct EigXℂincrement
    dofgr :: DofGroup
    p     :: 𝕔1
    Δx    :: 𝕔11
end
function blkasm(iA,jA,vA)
    pattern = sparse(iA,jA,vA)
    A,Aasm,Vasm,Vdis = prepare(pattern)
    zero!(A)
    for i ∈ eachindex(vA)
        addin!(Aasm,A,vA[i],iA[i],jA[i],1.)
    end
    return A,Vasm,Vdis
end 
function solve(::Type{EigX{ℂ}},pstate,verbose,dbg; 
    state::State, nmod::𝕫=5,droptol::𝕣=1e-9,kwargs...) 
    OX,OU,IA         = 2,0,0
    model,dis        = state.model,state.dis

    verbose && @printf("\n    Assembing\n")
    wantK,wantC,wantM= (ind.Λ,ind.X,1,1),(ind.Λ,ind.X,1,2),(ind.Λ,ind.X,1,3)
    wanted           = Wanted{1,OX+1,OU+1,IA}(:all,(wantK,wantC,wantM))
    out,asm,dofgr    = prepare(AssemblyDirect,model,dis,wanted)  
    nXdof            = getndof.(dofgr)[ind.X]
    state₀           = State{1,OX+1,OU+1}(copy(state))   
    assemble!{:matrices}(out,asm,dis,model,state₀,idmult,(dbg...,solver=:EigXℂ))
    K                = out.L2[ind.Λ,ind.X][1,1]
    C                = out.L2[ind.Λ,ind.X][1,2]
    M                = out.L2[ind.Λ,ind.X][1,3]
    I                = spdiagm(ones(nXdof))
    sparser!(K,droptol)
    sparser!(C,droptol)
    sparser!(M,droptol)
    A,_   ,_        = blkasm([1,2  ],[1,2  ],[ I,K  ])  
    B,_   ,_        = blkasm([2,1,2],[1,2,2],[-I,M,C])

    verbose && @printf("\n    Solving Eigenvalues\n")
    p, vec, ncv    = geneig{:complex}(A,B,nmod;kwargs...)
    ncv≥nmod||muscadeerror(dbg,@sprintf("eigensolver only converged for %i out of %i modes",ncv,nmod))
    Δx             = [view(vecᵢ,nXdof+1:2nXdof) for vecᵢ ∈ vec]

    pstate[] = EigXℂincrement(dofgr[ind.X],p,Δx)
    return 
end
function increment{OX}(initialstate,eiginc::EigXℂincrement,imod::AbstractVector{𝕫},A::AbstractVector) where{OX} 
    state            = State{1,OX+1,1}(copy(initialstate)) 
    maximum(imod)≤length(eiginc.p)|| muscadeerror(@sprintf("eiginc only has %n modes.",length(ω)))
    for i∈eachindex(imod)  
        pᵢ,Δxᵢ = eiginc.p[imod[i]],eiginc.Δx[imod[i]]
        for iOX     = 0:OX
            increment!(state,iOX+1,ℜ.((pᵢ)^iOX*A[i]*Δxᵢ),eiginc.dofgr)
        end
    end
    return state
end

