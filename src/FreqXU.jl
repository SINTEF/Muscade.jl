function makepattern(out) 
    L2(α,β) = out.L2[α,β][1,1]
    α       = [2,3,1,2,3,1,2,3]  #   [0 . .]
    β       = [1,1,2,2,2,3,3,3]  #   [. . .]
    return sparse(α,β,L2.(α,β))  # = [. . .]
end

function assemblebigmat!(L2::Vector{Sparse𝕣2},L2bigasm,asm,model,dis,out::AssemblyDirect{OX,OU,0},dbg) where{OX,OU}
    # does not call assemble!: solve has previously called assemble! to prepare bigasm, so out.L2 is already set,
    zero!.(L2)
    for     α ∈ λxu 
        for β ∈ λxu
            Lαβ = out.L2[α,β]
            for     αder = 1:size(Lαβ,1)
                for βder = 1:size(Lαβ,2)
                    ider =  αder+βder-1   
                    sgn  = isodd(αder) ? +1 : -1 
                    addin!(L2bigasm,L2[ider],Lαβ[αder,βder],α,β,sgn) 
                end
            end
        end
    end
end
function assemblebigvec!(L1::Vector{𝕣1},L1bigasm,asm,model,dis,out::AssemblyDirect{OX,OU,0},state,dbg) where{OX,OU}
    zero!.(L1)
    out.matrices = false
    assemble!(out,asm,dis,model,state,(dbg...,asm=:assemblebigmat!))
    for β ∈ λxu
        Lβ = out.L1[β]
        for βder = 1:size(Lβ,1)
            addin!(L2bigasm,L1[ider],Lβ[βder],β,scale) 
        end
    end
end
"""
	FreqXU{OX,OU}

"""
struct FreqXU{OX,OU} <: AbstractSolver end 

function solve(::Type{FreqXU{OX,OU}},pstate,verbose::𝕓,dbg;
    Δt::𝕣, p::𝕫, t₀::𝕣=0.,tᵣ::𝕣=t₀, 
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
    stateᵣ                = State{1,OX+1,OU+1}(copy(initialstate,time=tᵣ))   

    for (step,timeᵢ)      = enumerate(time)
        state[step]       = State(timeᵢ,deepcopy(stateᵣ.Λ),deepcopy(stateᵣ.X),deepcopy(stateᵣ.U),stateᵣ.A,nothing,stateᵣ.model,stateᵣ.dis)
    end
    L2                    = Vector{Sparse𝕣2}(undef,5)

    # Prepare assembler
    verbose && @printf("\n    Preparing assembler\n")
    out,asm,dofgr         = prepare(AssemblyDirect{OX,OU,IA},model,dis;fastresidual,kwargs...)   # model assembler for all arrays   
    assemble!(out,asm,dis,model,stateᵣ,(dbg...,solver=:FreqXU,phase=:matrices))            # assemble all model matrices - in class-blocks
    pattern               = makepattern(out)
    L2[1],L2bigasm,L1bigasm,L1dis  = prepare(pattern)                                            
    for ider = 2:5
        L2[ider] = copy(L2[1])
    end    
    assemblebigmat!(L2,L2bigasm,asm,model,dis,out,(dbg...,solver=:FreqXU))              # assemble all model matrices, no blocks

    # out.matrices = false
    # for (step,timeᵢ)∈enumerate(time)
    #     assemble!(out,asm,dis,model,state[step],(dbg...,solver=:FreqXU,phase=:matrices))

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


