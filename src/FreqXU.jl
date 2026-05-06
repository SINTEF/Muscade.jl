
#= 

TODO 
Does FFT take significant time? If so:

Avoid FFT of zeros, and addition of zeros
    find L1ᵢ that are all zero

Avoid FFT of zeros    
    in the non-all-zero L1ᵢ find dofs whose duals are zero over time
=#


"""
	FreqXU{OX,OU}

A linear frequency domain solver for optimisation FEM.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
stateXU         = solve(FreqXU{OX,OU};Δt, p, t₀,tᵣ,initialstate)
```

The solver linearises the problem (computes the Hessian of the Lagrangian) at `initialstate` with time `tᵣ`, and solves
it at times `t=range(start=t₀,step=Δt,length=2^p)`. The return


# Parameters
- `OX`                0 for static analysis
                      1 for first order problems in time (viscosity, friction, measurement of velocity)
                      2 for second order problems in time (inertia, measurement of acceleration) 
- `OU`                0 for white noise prior to the unknown load process
                      2 otherwise

# Named arguments
- `dbg=(;)`           a named tuple to trace the call tree (for debugging).
- `verbose=true`      set to false to suppress printed output (for testing).
- `silenterror=false` set to true to suppress print out of error (for testing) .
- `initialstate`      a `State`.
- `t₀=0.`             time of first step.                      
- `Δt`                time step.
- `p`                 `2^p` steps will be analysed.      
- `tᵣ=t₀`             reference time for linearisation.
- `droptol=1e-10`     set to zero terms in the incremental matrices that are smaller than `droptol` in absolute value.                      

# Output

A vector of length `2^p` containing the state of the model at each of these steps.                       

See also: [`solve`](@ref), [`initialize!`](@ref), [`study_singular`](@ref), [`SweepX`](@ref), [`DirectXUA`](@ref)
"""
struct FreqXU{OX,OU} <: AbstractSolver end 
function solve(::Type{FreqXU{OX,OU}},pstate,verbose::𝕓,dbg;
    Δt::𝕣, p::𝕫, t₀::𝕣=0.,tᵣ::𝕣=t₀, 
    initialstate::State,
    droptol::𝕣=1e-10,
    wantedhessians::Matrix{𝕓2}=want_all_hessians(1,OX+1,OU+1,0)) where{OX,OU}

    #  Mostly constants
    local LU
    model,dis             = initialstate.model, initialstate.dis
    nω                    = 2^(p-1)
    nstep                 = 2nω
    time                  = range(;start=t₀,step=Δt,length=nstep)
    IA                    = 0

    # State storage
    S                     = State{1,3,3,Nothing}
    pstate[] = state      = Vector{S}(undef,nstep)                                                                           
    stateᵣ                = State{1,3,3}(copy(initialstate,time=tᵣ))   

    # Prepare assembler
    verbose && @printf("    Preparing assembler\n")
    out,asm,dofgr         = prepare(AssemblyDirect{OX,OU,IA},model,dis,wantedhessians)   # model assembler for all arrays   

    verbose && @printf("    Computing matrices\n")
    assemble!{:matrices}(out,asm,dis,model,stateᵣ,Δt,(dbg...,solver=:FreqXU,phase=:matrices))            # assemble all model matrices - in class-blocks
    pattern               = make_λxu_sparsepattern(out)
    L2                    = Vector{Sparse𝕣2}(undef,5)
    L2[1],L2bigasm,L1bigasm,Ldis  = prepare(pattern)  
    λxu_dofgr             = allΛXUdofs(model,dis)                                        # NB same ordering of dofs in rhs as implied by pattern                                          
    for ider              = 2:5
        L2[ider]          = copy(L2[1])
    end    
    assemblebigmat!(L2,L2bigasm,asm,model,dis,out,(dbg...,solver=:FreqXU))              # assemble all complete model matrices into L2

    verbose && @printf("    Improving sparsity ")    
    keep = sparser!(L2,droptol)
    verbose && @printf("from %i to %i nz terms\n",length(keep),sum(keep))    

    verbose && @printf("    Computing rhs\n")
    ndof                  = size(L2[1],1)
    L1𝕔                   = ntuple(ider->𝕔2(undef,nω,ndof)       ,3)
    L1𝕣                   = ntuple(ider->reinterpret(𝕣,L1𝕔[ider]),3)
    #TODO Multithread
    for (step,timeᵢ)      = enumerate(time)
        L1ᵢ               = ntuple(ider->view(L1𝕣[ider],step,:),3)
        state[step]       = State(timeᵢ,deepcopy(stateᵣ.Λ),deepcopy(stateᵣ.X),deepcopy(stateᵣ.U),stateᵣ.A,nothing,stateᵣ.model,stateᵣ.dis)
        assemblebigvec!(L1ᵢ,L1bigasm,asm,model,dis,out,state[step],Δt,dbg)
    end
  
    verbose && @printf("    Fourier transform of rhs\n")
    for L1ᵢ∈ L1𝕔
        𝔉!(L1ᵢ,Δt)
    end
    Δω  = getδω(nstep,Δt)
    ω   = range(start=0.,step=Δω,length=nω)

    verbose && @printf("    Solving equations for all ω\n")
    local LU
    x   = L2[1]
    M   = Sparse𝕔2(ndof,ndof,x.colptr,x.rowval,𝕔1(undef,length(x.nzval)))
    rhs = 𝕔1(undef,ndof)
    Δz  = 𝕔1(undef,ndof)

    # TODO multithread
    for (iω,ωᵢ) = enumerate(ω)
        for inz ∈eachindex(M.nzval)
            M.nzval[inz] = L2[1].nzval[inz] 
        end
        for j = 1:4
            a = (𝑖*ωᵢ)^j
            for inz ∈eachindex(M.nzval)
                M.nzval[inz] += a *L2[j+1].nzval[inz]
            end
        end
        for idof ∈eachindex(rhs)
            rhs[idof]   = L1𝕔[1][iω,idof]
        end
        for j = 1:2
            a = (-𝑖*ωᵢ)^j
            for idof ∈eachindex(rhs)
                rhs[idof] += a *L1𝕔[j+1][iω,idof]
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
        ldiv!(Δz,LU,rhs)
        for (ider,L1ᵢ) ∈ enumerate(L1𝕔)
            L1ᵢ[iω,:] .= Δz * (𝑖*ωᵢ)^(ider-1)
        end
    end    

    verbose && @printf("    Inverse Fourier transform of solution and its time derivatives\n")
    for L1ᵢ∈ L1𝕔
        𝔉⁻¹!(L1ᵢ,Δω)
    end

    verbose && @printf("    Updating the states\n")
    # TODO multithread
    for (step,stateᵢ) = enumerate(state)
        for ider = 1:3
            decrement!(stateᵢ,ider,view(L1𝕣[ider],step,:),λxu_dofgr)
        end
    end
    verbose && @printf("\n\n")
    return
end


