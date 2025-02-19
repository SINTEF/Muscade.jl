function makepattern(out) 
    L2(α,β) = out.L2[α,β][1,1]
    α       = [2,3,1,2,3,1,2,3]  #   [0 . .]
    β       = [1,1,2,2,2,3,3,3]  #   [. . .]
    return sparse(α,β,L2.(α,β))  # = [. . .]
end

const Mder = (((1,1),           ),
              ((1,2),(2,1)      ),
              ((1,3),(2,2),(3,1)),
              ((3,2),(2,3)      ),
              ((3,3),           ) )

function assemblebigmat!(L2::Vector{Sparse𝕣2},L2bigasm,asm,model,dis,out::AssemblyDirect{OX,OU,0},state,dbg) where{OX,OU}
    scale = 1 ###########
    zero!.(L2)
    out.matrices = true
    assemble!(out,asm,dis,model,state,(dbg...,asm=:assemblebigmat!))
    for  mder ∈ Mder
        for (ider,mderᵢ)∈enumerate(mder)
            sgn = isodd(ider) ? +1 : -1
            for     α ∈ λxu 
                for β ∈ λxu
                    Lαβ = out.L2[α,β]
    
                    @show size(Lαβ) # by design, unrequired derivatives are not stored.
                    addin!(L2bigasm,L2[ider],Lαβ[mderᵢ...],α,β,sgn*scale) 
                end
            end
        end
    end 
end

"""
	FreqXU{OX,OU}

"""
struct FreqXU{OX,OU} <: AbstractSolver end 

function solve(::Type{FreqXU{OX,OU}},pstate,verbose::𝕓,dbg;
    Δt::𝕣, p::𝕫, t₀::𝕣=0., 
    initialstate::State,
    fastresidual::𝔹=false,
    kwargs...) where{OX,OU}

    #  Mostly constants
    local LU
    #nder                  = (1,OX+1,OU+1)
    model,dis             = initialstate.model, initialstate.dis
    nstep                 = 2^p
    time                  = range(;start=t₀,step=Δt,length=nstep)
    IA                    = 0

    # State storage
    S                     = State{1,OX+1,OU+1,Nothing}
    pstate[] = state      = Vector{S}(undef,nstep)                                                                           
    state₁                = State{1,OX+1,OU+1}(copy(initialstate,time=t₀))   

    for (step,timeᵢ)      = enumerate(time)
        state[step]       = step==1 ? state₁ : State(timeᵢ,deepcopy(state₁.Λ),deepcopy(state₁.X),deepcopy(state₁.U),state₁.A,nothing,state₁.model,state₁.dis)
    end
    L2                    = Vector{Sparse𝕣2}(undef,5)

    # Prepare assembler
    verbose && @printf("\n    Preparing assembler\n")
    out,asm,dofgr         = prepare(AssemblyDirect{OX,OU,IA},model,dis;fastresidual,kwargs...)      
    assemble!(out,asm,dis,model,state[1],(dbg...,solver=:DirectXUA,phase=:sparsity))     # assemble all model matrices - in blocks
    pattern               = makepattern(out)
    L2[1],L2bigasm,L1bigasm,L1dis  = prepare(pattern)

    for ider = 2:5
        L2[ider] = copy(L2[1])
    end    
    assemblebigmat!(L2,L2bigasm,asm,model,dis,out,state[1],(dbg...,solver=:FreqXU))

    @show L2[1]
    @show L2[2]
    @show L2[3]
    @show L2[4]
    @show L2[5]


        # try 
        #     if iter==1 LU = lu(Lvv) 
        #     else       lu!(LU ,Lvv)
        #     end 
        # catch 
        #     verbose && @printf("\n")
        #     muscadeerror(@sprintf("Lvv matrix factorization failed at iter=%i",iter));
        # end
        # Δv               = LU\Lv # use ldiv! to save allocation

        # decrementbigmat!(state,Δ²,Lvdis,dofgr,Δv,nder,Δt,nstep)

    return
end


