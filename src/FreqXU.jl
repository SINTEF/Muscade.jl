# TODO: everything!!!

const λxu    = 1:3
const ind    = (Λ=1,X=2,U=3,A=4)
const nclass = length(ind) 

## Assembly of sparse
arrnum(α  )  =          α
arrnum(α,β)  = nclass + β + nclass*(α-1) 
mutable struct AssemblyFreqXU{OX,OU,T1,T2}  <:Assembly
    L1 :: T1   
    L2 :: T2   
    fastresidual :: 𝔹
end  
function prepare(::Type{AssemblyFreqXU{OX,OU}},model,dis;Xwhite=false,XUindep=false,UAindep=false,XAindep=false,fastresidual=false) where{OX,OU}
    dofgr    = (allΛdofs(model,dis),allXdofs(model,dis),allUdofs(model,dis),allAdofs(model,dis))
    ndof     = getndof.(dofgr)
    neletyp  = getneletyp(model)
    asm      = Matrix{𝕫2}(undef,nclass+nclass^2,neletyp)
    nder     = (1,OX+1,OU+1)
    L1 = Vector{Vector{Vector{𝕣}}}(undef,4)
    for α∈λxua
        nα = nder[α]
        av = asmvec!(view(asm,arrnum(α),:),dofgr[α],dis)
        L1[α] = Vector{Vector{𝕣}}(undef,nα)
        for αder=1:nα 
            L1[α][αder] = copy(av)
        end
    end
    L2 = Matrix{Matrix{SparseMatrixCSC{Float64, Int64}}}(undef,4,4)
    for α∈λxua, β∈λxua
        am = asmmat!(view(asm,arrnum(α,β),:),view(asm,arrnum(α),:),view(asm,arrnum(β),:),ndof[α],ndof[β])
        nα,nβ = nder[α], nder[β]
        if            α==β==ind.Λ          nα,nβ=0,0 end   # Lλλ is always zero
        if Xwhite  && α==β==ind.X          nα,nβ=1,1 end   # X-measurement error is white noise process
        if XUindep && α==ind.X && β==ind.U nα,nβ=0,0 end   # X-measurements indep of U
        if XUindep && α==ind.U && β==ind.X nα,nβ=0,0 end   # X-measurements indep of U
        if XAindep && α==ind.X && β==ind.A nα,nβ=0,0 end   # X-measurements indep of A
        if XAindep && α==ind.A && β==ind.X nα,nβ=0,0 end   # X-measurements indep of A
        if UAindep && α==ind.U && β==ind.A nα,nβ=0,0 end   # U-load indep of A
        if UAindep && α==ind.A && β==ind.U nα,nβ=0,0 end   # U-load indep of A
        L2[α,β] = Matrix{SparseMatrixCSC{Float64, Int64}}(undef,nα,nβ)
        for αder=1:nα,βder=1:nβ
            L2[α,β][αder,βder] = copy(am)
        end
    end
    out      = AssemblyFreqXU{OX,OU,typeof(L1),typeof(L2)}(L1,L2,fastresidual)
    return out,asm,dofgr
end
function zero!(out::AssemblyFreqXU)
    for α∈λxua 
        zero!.(out.L1[α])
        for β∈λxua
            zero!.(out.L2[α,β])
        end
    end
end
function addin!(out::AssemblyFreqXU{OX,OU,T1,T2},asm,iele,scale,eleobj::Eleobj,  Λ::NTuple{1  ,SVector{Nx}},
                                                                                         X::NTuple{NDX,SVector{Nx}},
                                                                                         U::NTuple{NDU,SVector{Nu}},
                                                                                         A::           SVector{Na} ,t,SP,dbg) where{OX,OU,NDX,NDU,T1,T2,Nx,Nu,Na,Eleobj} 
    @assert NDX==OX+1 @sprintf("got OX=%i and NDX=%i. Expected OX+1==NDX",OX,NDX)
    @assert NDX==OX+1 @sprintf("got OU=%i and NDU=%i. Expected OU+1==NDU",OU,NDU)
    ndof   = (Nx, Nx, Nu)
    nder   = (1,OX+1,OU+1)
    Np     = Nx + Nx*(OX+1) + Nu*(OU+1)  # number of partials

    Λ∂ =              SVector{Nx}(∂²ℝ{1,Np}(Λ[1   ][idof],                           idof, scale.Λ[idof])   for idof=1:Nx)
    X∂ = ntuple(ider->SVector{Nx}(∂²ℝ{1,Np}(X[ider][idof],Nx+Nx*(ider-1)            +idof, scale.X[idof])   for idof=1:Nx),OX+1)
    U∂ = ntuple(ider->SVector{Nu}(∂²ℝ{1,Np}(U[ider][idof],Nx+Nx*(OX+1)  +Nu*(ider-1)+idof, scale.U[idof])   for idof=1:Nu),OU+1)
    L,FB = getlagrangian(eleobj, Λ∂,X∂,U∂,A ,t,SP,dbg)
    ∇L           = ∂{2,Np}(L)
    pα           = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
    for α∈λxua, i=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
        iα       = pα.+(1:ndof[α])
        pα      += ndof[α]
        Lα = out.L1[α]
        if i≤size(Lα,1)  # ...but only add into existing vectors of L1, for speed
            add_value!(out.L1[α][i] ,asm[arrnum(α)],iele,∇L,iα)
        end
        pβ       = 0
        for β∈λxua, j=1:nder[β]
            iβ   = pβ.+(1:ndof[β])
            pβ  += ndof[β]
            Lαβ = out.L2[α,β]
            if i≤size(Lαβ,1) && j≤size(Lαβ,2) # ...but only add into existing matrices of L2, for better sparsity
                add_∂!{1}(out.L2[α,β][i,j],asm[arrnum(α,β)],iele,∇L,iα,iβ)
            end
        end
    end
end

"""
	FreqXU{OX,OU}

A linear frequency domain solver for optimisation FEM.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
state           = solve(FreqXU{OX,OU};initialstate,time=0:1.:5)
```

The solver does not yet support interior point methods. 

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
Setting the following flags to `true` will improve the sparsity of the system. But setting
a flag to `true` when the condition isn't met causes the Hessian to be wrong, which is detrimental for convergence.                      
- `Xwhite=false`      `true` if response measurement error is a white noise process.
- `XUindep=false`     `true` if response measurement error is independant of `U`

# Output

A vector of length equal to that of `time` containing the state of the optimized model at each of these steps.                       

See also: [`solve`](@ref), [`SweepX`](@ref), [`FreqXU`](@ref)
"""
struct FreqXU{OX,OU} <: AbstractSolver end 
function solve(TS::Type{FreqXU{OX,OU}},pstate,verbose::𝕓,dbg;
    time::AbstractRange{𝕣},
    initialstate::State,
    maxiter::ℤ=50,
    maxΔλ::ℝ=1e-5,maxΔx::ℝ=1e-5,maxΔu::ℝ=1e-5,maxΔa::ℝ=1e-5,
    saveiter::𝔹=false,
    fastresidual:: 𝔹=false,
    kwargs...) where{OX,OU}

    #  Mostly constants
    local LU
    nstep                 = length(time)
    Δt                    = (last(time)-first(time))/(nstep-1)
    γ                     = 0.
    nder                  = (1,OX+1,OU+1)
    model,dis             = initialstate.model, initialstate.dis
    if IA==1  Δ², maxΔ²   = 𝕣1(undef,4), [maxΔλ^2,maxΔx^2,maxΔu^2,maxΔa^2] 
    else      Δ², maxΔ²   = 𝕣1(undef,3), [maxΔλ^2,maxΔx^2,maxΔu^2        ] 
    end

    # State storage
    S                     = State{1,OX+1,OU+1,@NamedTuple{γ::Float64,iter::Int64}}
    state                 = Vector{S}(undef,nstep)
    s                     = State{1,OX+1,OU+1}(copy(initialstate,time=time[1],SP=(γ=0.,iter=1)))   
    for (step,timeᵢ)      = enumerate(time)
        state[step]       = step==1 ? s : State(timeᵢ,deepcopy(s.Λ),deepcopy(s.X),deepcopy(s.U),s.A,s.SP,s.model,s.dis)
    end
    if saveiter
        stateiter         = Vector{Vector{S}}(undef,maxiter) 
        pstate[]          = stateiter
    else
        pstate[]          = state                                                                            
    end    

    # Prepare assembler
    verbose && @printf("\n    Preparing assembler\n")
    out,asm,dofgr         = prepare(AssemblyFreqXU{OX,OU},model,dis;fastresidual,kwargs...)      # mem and assembler for system at any given step
    assemble!(out,asm,dis,model,state[1],(dbg...,solver=:FreqXU,phase=:sparsity))     # create a sample "out" for preparebig
    Lvv,Lv,Lvvasm,Lvasm,Lvdis = preparebig(OX,OU,nstep,out)                             # mem and assembler for big system

    for iter              = 1:maxiter
        verbose && @printf("\n    Iteration %3d\n",iter)

        verbose && @printf("        Assembling")
        SP = (γ=γ,iter=iter)
        assemblebig!(Lvv,Lv,Lvvasm,Lvasm,asm,model,dis,out,state,nstep,Δt,SP,(dbg...,solver=:FreqXU,iter=iter))

        verbose && @printf(", solving")
        try 
            if iter==1 LU = lu(Lvv) 
            else       lu!(LU ,Lvv)
            end 
        catch 
            verbose && @printf("\n")
            muscadeerror(@sprintf("Lvv matrix factorization failed at iter=%i",iter));
        end
        Δv               = LU\Lv # use ldiv! to save allocation

        verbose && @printf(", decrementing.\n")
        decrementbig!(state,Δ²,Lvdis,dofgr,Δv,nder,Δt,nstep)
        
        if saveiter
            stateiter[iter]     = copy.(state) 
        end
        verbose          && @printf(  "        maxₜ(|ΔΛ|)=%7.1e ≤ %7.1e  \n",√(Δ²[ind.Λ]),√(maxΔ²[ind.Λ]))
        verbose          && @printf(  "        maxₜ(|ΔX|)=%7.1e ≤ %7.1e  \n",√(Δ²[ind.X]),√(maxΔ²[ind.X]))
        verbose          && @printf(  "        maxₜ(|ΔU|)=%7.1e ≤ %7.1e  \n",√(Δ²[ind.U]),√(maxΔ²[ind.U]))
        verbose && IA==1 && @printf(  "             |ΔA| =%7.1e ≤ %7.1e  \n",√(Δ²[ind.A]),√(maxΔ²[ind.A]))
        if all(Δ².≤maxΔ²)  
            verbose      && @printf("\n    Converged in %3d iterations.\n",iter)
            verbose      && @printf(  "    nel=%d, nvar=%d, nstep=%d\n",getnele(model),length(Lv),nstep)
            break#out of iter
        end
        iter<maxiter || muscadeerror(@sprintf("no convergence after %3d iterations. \n",iter))
    end # for iter
    return
end


