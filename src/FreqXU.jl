
#= 

TODO 

Faster LU decomposition
    sparsify L2 (find entries that are zero in L2ᵢ ∀i, and update the sparsity structure)

Avoid FFT of zeros, and addition of zeros
    find L1ᵢ that are all zero

Avoid FFT of zeros    
    in the non-all-zero L1ᵢ find dofs whose duals are zero over time
=#


function makepattern(out) 
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


"""
	FreqXU{OX,OU}

"""
struct FreqXU{OX,OU} <: AbstractSolver end 
using GLMakie
function solve(::Type{FreqXU{OX,OU}},pstate,verbose::𝕓,dbg;
    Δt::𝕣, p::𝕫, t₀::𝕣=0.,tᵣ::𝕣=t₀, 
    initialstate::State,
    fastresidual::𝔹=false,
    kwargs...) where{OX,OU}

    #  Mostly constants
    local LU
    model,dis             = initialstate.model, initialstate.dis
    nω                    = 2^(p-1)
    nstep                 = 2*nω
    time                  = range(;start=t₀,step=Δt,length=nstep)
    IA                    = 0

    # State storage
    S                     = State{1,3,3,Nothing}
    pstate[] = state      = Vector{S}(undef,nstep)                                                                           
    stateᵣ                = State{1,3,3}(copy(initialstate,time=tᵣ))   

    # Prepare assembler
    verbose && @printf("    Preparing assembler\n")
    out,asm,dofgr         = prepare(AssemblyDirect{OX,OU,IA},model,dis;fastresidual,kwargs...)   # model assembler for all arrays   

    verbose && @printf("    Computing matrices\n")
    out.matrices          = true
    assemble!(out,asm,dis,model,stateᵣ,(dbg...,solver=:FreqXU,phase=:matrices))            # assemble all model matrices - in class-blocks
    pattern               = makepattern(out)
    L2                    = Vector{Sparse𝕣2}(undef,5)
    L2[1],L2bigasm,L1bigasm,Ldis  = prepare(pattern)  
    λxu_dofgr             = allΛXUdofs(model,dis)                                        # NB same ordering of dofs in rhs as implied by pattern                                          
    for ider              = 2:5
        L2[ider]          = copy(L2[1])
    end    
    assemblebigmat!(L2,L2bigasm,asm,model,dis,out,(dbg...,solver=:FreqXU))              # assemble all complete model matrices into L2

    verbose && @printf("    Computing rhs\n")
    ndof                  = size(L2[1],1)
    L1𝕔                   = ntuple(ider->𝕔2(undef,nω,ndof)       ,3)
    L1𝕣                   = ntuple(ider->reinterpret(𝕣,L1𝕔[ider]),3)
    out.matrices          = false
    #TODO Multithread
    for (step,timeᵢ)      = enumerate(time)
        L1ᵢ               = ntuple(ider->view(L1𝕣[ider],step,:),3)
        state[step]       = State(timeᵢ,deepcopy(stateᵣ.Λ),deepcopy(stateᵣ.X),deepcopy(stateᵣ.U),stateᵣ.A,nothing,stateᵣ.model,stateᵣ.dis)
        assemblebigvec!(L1ᵢ,L1bigasm,asm,model,dis,out,state[step],dbg)
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


