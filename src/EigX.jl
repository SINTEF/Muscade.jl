"""
	res = solve(EigX{ℝ};state=initialstate,nmod)
	res = solve(EigX{ℂ};state=initialstate,nmod)

Given an initial (typicaly static) state `initialstate`, computes the lowest `nmod` eigenmodes of a system. 
`EigX{ℝ}` computes real eigenmodes not accounting for damping. `EigX{ℝ}` computes complex eigenmodes accounting for damping.
The data structure `res` can be passed to `increment` to obtain dynamic states superimposing mode shapes.

# Input
- `initialstate` - a `State`
- `nmod=5` - the number of eigenmodes to identify
- `droptol=1e-9` - in the stiffness and mass matrix, the magnitude of a term relative to the largest term in the matrix
        under which the term is set to zero.
- Further named arguments: see the optional keyword arguments to `geneig`.         

# Output
A `NamedTuple` with fields
- `ω` (for `EigX{ℝ}`) a `Vector` of length `nmod` or longer containing circular frequencies associated to each mode.
- `p` (for `EigX{ℂ}`) a `Vector` of length `nmod` or longer containing complex exponents `p=τ+𝑖ω` associated to each mode.
- `vecs` is an internal.
- `dofgr` is an internal.

See also: [`solve`](@ref), [`initialize!`](@ref), [`increment`](@ref), [`geneig`](@ref)
"""
struct        EigX{T} <: AbstractSolver end
function solve(::Type{EigX{ℝ}},pstate,verbose,dbg; 
                   state::State, nmod::𝕫=5,droptol::𝕣=1e-9,kwargs...) 
    OX,OU,IA         = 2,0,0
    model,dis        = state.model,state.dis

    verbose && @printf("\n    Assembing\n")
    out,asm,dofgr    = prepare(AssemblyDirect{OX,OU,IA},model,dis)  
    nXdof            = getndof.(dofgr)[ind.X]
    state            = State{1,OX+1,OU+1}(copy(state))   
    assemble!(out,asm,dis,model,state,(dbg...,solver=:EigXℝ))
    K                = out.L2[ind.Λ,ind.X][1,1]
    M                = out.L2[ind.Λ,ind.X][1,3]
    sparser!(K,droptol)
    sparser!(M,droptol)

    verbose && @printf("\n    Solving Eigenvalues\n")
    ω², vecs, ncv = geneig{:SDP}(K,M,nmod;kwargs...)
    ncv≥nmod||muscadeerror(dbg,@sprintf("eigensolver only converged for %i out of %i modes",ncv,nmod))
    pstate[] = (dofgr=dofgr[ind.X],ω=sqrt.(ω²),v=vecs)
    return 
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
    out,asm,dofgr    = prepare(AssemblyDirect{OX,OU,IA},model,dis)  
    nXdof            = getndof.(dofgr)[ind.X]
    state            = State{1,OX+1,OU+1}(copy(state))   
    assemble!(out,asm,dis,model,state,(dbg...,solver=:EigXℂ))
    K                = out.L2[ind.Λ,ind.X][1,1]
    C                = out.L2[ind.Λ,ind.X][1,2]
    M                = out.L2[ind.Λ,ind.X][1,3]
    I                = spdiagm(ones(nXdof))
    sparser!(K,droptol)
    sparser!(C,droptol)
    sparser!(M,droptol)
    A,_   ,_        = blkasm([1,2  ],[1,2  ],[ I,K  ])  # something wrong in my blkasm process, possibly initialisation
    B,_   ,_        = blkasm([2,1,2],[1,2,2],[-I,M,C])

    verbose && @printf("\n    Solving Eigenvalues\n")
    p, vecs, ncv    = geneig{:Complex}(A,B,nmod;kwargs...)
    ncv≥nmod||muscadeerror(dbg,@sprintf("eigensolver only converged for %i out of %i modes",ncv,nmod))
    v               = [view(vecsᵢ,nXdof+1:2nXdof) for vecsᵢ∈vecs]

    pstate[] = (dofgr=dofgr[ind.X],p=p,v=v)
    return 
end

"""
    state = increment(initialstate,res,imod,A;order=2)

Starting from `initalstate` for which an `EigX` analysis has been carried out, and using the output
`res` of that analysis, construct new `State`s representing the instantaneous state of the 
vibrating structure
    
# Input
- `initialstate` the same initial `State` provided to `EigX` to compute `res`
- `res` obtained from `EigX`
- `imod`, an `AbstractVector` of integer mode numbers
- `A`, an `AbstractVector` of same length as `imod`, containing real or complex 
  amplitudes associated to the modes
- `order=2` the number of time derivatives to be computed

# Output
- `state` a snapshot of the vibrating system

See also: [`EigX`](@ref)

"""
function increment(initialstate,res::Tres,imod::AbstractVector{𝕫},A::AbstractVector;order=2) where{Tres}
    length(imod)==length(A)|| muscadeerror("imod and A must be of same length.")
    state            = State{1,order+1,1}(copy(initialstate)) 
    if hasfield(Tres,:ω)
        maximum(imod)≤length(res.ω)|| muscadeerror(@sprintf("res only has %n modes.",length(ω)))
        for i∈eachindex(imod)  
            ωᵢ,vᵢ = res.ω[imod[i]],res.v[imod[i]]
            for n     = 0:order
                increment!(state,n+1,ℜ.((𝑖*ωᵢ)^n*A[i]*vᵢ),res.dofgr)
            end
        end
    elseif hasfield(Tres,:p)
        maximum(imod)≤length(res.p)|| muscadeerror(@sprintf("res only has %n modes.",length(ω)))
        for i∈eachindex(imod)  
            pᵢ,vᵢ = res.p[imod[i]],res.v[imod[i]]
            for n     = 0:order
                increment!(state,n+1,ℜ.((pᵢ)^n*A[i]*vᵢ),res.dofgr)
            end
        end
    else
        muscadeerror("To show results from eigenvalue analysis, expected field ω or p")
    end
    return state
end
