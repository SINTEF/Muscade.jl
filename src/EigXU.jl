function make_λxu_sparsepattern(out) 
    L2(α,β) = out.L2[α,β][1,1]
    α       = [2,3,1,2,3,1,2,3]  #   [0 . .]
    β       = [1,1,2,2,2,3,3,3]  #   [. . .]
    return sparse(α,β,L2.(α,β))  # = [. . .]
end

function assemblebigmat!(L2::Vector{Sparse𝕣2},L2bigasm,asm,model,dis,out::AssemblyDirect{OX,OU,0},dbg) where{OX,OU}
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
                    sgn  = isodd(αder) ? +1 : -1 
                    if α==ind.Λ && β==ind.U
                        addin!(L2bigasm,L2[ider],Lαβ[αder,βder],α,β,sgn) 
                    else
                        addin!(L2bigasm,L2[ider],Lαβ[αder,βder],α,β,sgn) 
                    end
                end
            end
        end
    end
end
function assemblebigvec!(L1,L1bigasm,asm,model,dis,out::AssemblyDirect{OX,OU,0},state,dbg) where{OX,OU}
    zero!.(L1)
    out.matrices = false
    assemble!(out,asm,dis,model,state,(dbg...,asm=:assemblebigvec!))
    for β ∈ λxu
        Lβ = out.L1[β]
        for βder = 1:size(Lβ,1)
            addin!(L1bigasm,L1[βder],Lβ[βder],β,1) 
        end
    end
end

struct EigXUincrement
    dofgr :: DofGroup
    λ     :: 𝕣11         # [iω][imod]
    Δz    :: Vector{𝕣11} # [iω][imod][idof]
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
- `initialstate`      a `State`.
- `nmod`              the number of eigen-modes to identusy
- `Δω`                frequency step
- `p`                 `2^p` steps will be analysed.      
- `droptol=1e-10`     set to zero terms in the incremental matrices that are smaller than `droptol` in absolute value.                      

# Output
- an object of type `EigXUincrement` for use with [`increment`](@ref) to create a snapshot of the
  oscillating system.

See also: [`solve`](@ref), [`initialize!`](@ref), [`studysingular`](@ref), [`SweepX`](@ref), [`DirectXUA`](@ref)
"""
struct EigXU{OX,OU} <: AbstractSolver end 
function solve(::Type{EigXU{OX,OU}},pstate,verbose::𝕓,dbg;
    Δω::𝕣, p::𝕫, 
    initialstate::State,
    droptol::𝕣=1e-10,
    nmod::𝕫=5,
    kwargs...) where{OX,OU}

    #  Mostly constants
    local LU
    model,dis             = initialstate.model, initialstate.dis
    nω                    = 2^(p-1)
    IA                    = 0

    # State storage
    S                     = State{1,3,3,Nothing}
    pstate[] = state      = Vector{S}(undef,nω)                                                                           
    state₀                = State{1,3,3}(copy(initialstate))   

    verbose && @printf("    Preparing assembler\n")
    out,asm,dofgr         = prepare(AssemblyDirect{OX,OU,IA},model,dis)   # model assembler for all arrays   

    verbose && @printf("    Computing matrices\n")
    out.matrices          = true
    assemble!(out,asm,dis,model,state₀,(dbg...,solver=:EigXU,phase=:matrices))            # assemble all model matrices - in class-blocks
    pattern               = make_λxu_sparsepattern(out)
    L2                    = Vector{Sparse𝕣2}(undef,5)
    L2[1],L2bigasm,L1bigasm,Ldis  = prepare(pattern)  
    λxu_dofgr             = allΛXUdofs(model,dis)                                        # NB same ordering of dofs in rhs as implied by pattern                                          
    for ider              = 2:5
        L2[ider]          = copy(L2[1])
    end    
    assemblebigmat!(L2,L2bigasm,asm,model,dis,out,(dbg...,solver=:EigXU))              # assemble all complete model matrices into L2
    sparser!(L2,droptol)
    nXdof,nUdof = getndof(model,(:X,:U))
    ixu = (nXdof+1):(2nXdof+nUdof)
    N   = sparse(ixu,ixu,ones(nXdof+nUdof))

    verbose && @printf("    Improving sparsity ")    
    keep = [any(abs(L2ⱼ.nzval[i])>droptol for L2ⱼ∈L2) for i∈eachindex(L2[1].nzval)]
    for L2ⱼ∈L2
        sparser!(L2ⱼ,j->keep[j])
    end
    verbose && @printf("from %i to %i nz terms\n",length(keep),sum(keep))    

    verbose && @printf("    Solving XU-eigenproblem for all ω\n")
    L2₁  = L2[1]
    ndof = 2nXdof+nUdof
    M    = Sparse𝕔2(ndof,ndof,L2₁.colptr,L2₁.rowval,𝕔1(undef,length(L2₁.nzval)))
    Δz   = Vector{𝕣11}(undef,nω) # Δz[iω][imod][idof]
    λ    = 𝕣11(undef,nω)         # λ[iω][imod]

    ω   = range(start=0.,step=Δω,length=nω)
    for (iω,ωᵢ) = enumerate(ω)
        for inz ∈eachindex(M.nzval)
            M.nzval[inz] = L2[1].nzval[inz] 
        end
        for j = 1:4
            𝑖ωᵢʲ  = (𝑖*ωᵢ)^j
            for inz ∈eachindex(M.nzval)
                M.nzval[inz] += 𝑖ωᵢʲ *L2[j+1].nzval[inz]
            end
        end
        try 
            if iω==1 LU = lu(M) 
            else     lu!(LU ,M)
            end 
        catch 
            verbose && @printf("\n")
            muscadeerror(@sprintf("M matrix factorization failed for ω=%f",ωᵢ));
        end

        λ[iω], Δz[iω], ncv = geneig{:Hermitian}(LU,N,nmod;kwargs...)
        # error message if ncv < nω?
    end    
    pstate[] = EigXUincrement(allΛXUdofs(model,dis),λ,Δz)
    verbose && @printf("\n\n")
    return
end

function increment{OX}(initialstate,eiginc::EigXUincrement,iω::𝕫,imod::AbstractVector{𝕫},A::AbstractVector) where{OX} 
    state            = State{1,OX+1,1}(copy(initialstate)) 
    maximum(imod)≤length(eiginc.λ)|| muscadeerror(@sprintf("eiginc only has %n modes.",length(ω)))
    for i∈eachindex(imod)  
        λᵢ,Δzᵢ = eiginc.λ[iω][imod[i]],eiginc.Δz[iω][imod[i]]
        for iOX = 0:OX
            increment!(state,iOX+1,ℜ.((λᵢ)^iOX*A[i]*Δzᵢ),eiginc.dofgr)
        end
    end
    return state
end