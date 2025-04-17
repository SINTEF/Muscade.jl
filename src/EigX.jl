"""
	res = solve(EigX;state=initialstate,nmod)

Given an initial (typicaly static) state `initialstate`, computes the lowest `nmod` eigenmodes. 
The data structure `res` can be passed to `increment` to obtain dynamic states superimposing
mode shapes.

# Input
- `initialstate` - a `State`
- `nmod=5` - the number of eigenmodes to identify
- `fastresidual=true` - limit automatic differentiation of `residual` to first order
- `droptol=1e-9` - in the stiffness and mass matrix, the magnitude of a term relative to the largest term in the matrix
        under which the term is set to zero.
- Further named arguments: see the optional keyword arguments to `geneig`.         

# Output
A `NamedTuple` with fields
- `ω` a `Vector` of length `nmod` or longer containing circular frequencies associated to each mode
- `vecs` is an internal
- `dofgr` is an internal

See also: [`solve`](@ref), [`initialize!`](@ref), [`increment`](@ref), [`geneig`](@ref)
"""
struct        EigX <: AbstractSolver end
function solve(TS::Type{EigX},pstate,verbose,dbg; 
                   state::State, nmod::𝕫=5,fastresidual::𝕓=true,droptol::𝕣=1e-9,kwargs...) 
    OX,OU,IA         = 2,0,0
    model,dis        = state.model,state.dis

    verbose && @printf("\n    Preparing assembler\n")
    out,asm,dofgr    = prepare(AssemblyDirect{OX,OU,IA},model,dis;fastresidual)  
    nXdof            = getndof.(dofgr)[ind.X]
    state            = State{1,OX+1,OU+1}(copy(state))   
    assemble!(out,asm,dis,model,state,(dbg...,solver=:EigX))
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
function increment(initialstate,res,imod::AbstractVector{𝕫},A::AbstractVector;order=2)
    length(imod)==length(A)|| muscadeerror("imod and A must be of same length.")
    maximum(imod)≤length(res.ω)|| muscadeerror(@sprintf("res only has %n modes.",length(res.ω)))
    state            = State{1,order+1,1}(copy(initialstate)) 
    for i∈eachindex(imod)  
        ω,v = res.ω[imod[i]],res.v[imod[i]]
        for n     = 0:order
            increment!(state,n+1,ℜ.((𝑖*ω)^n*A[i]*v),res.dofgr)
        end
    end
    return state
end
