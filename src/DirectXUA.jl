# dis.dis[ieletyp].index[iele].X|U|A[ieledof]       - disassembling model state into element dofs
# dis.dis[ieletyp].scale.Λ|X|U|A[ieledof]           - scaling each element type 
# dis.scaleΛ|X|U|A[imoddof]                         - scaling the model state
# dis.field  X|U|A[imoddof]                         - field of dofs in model state
# asm[iarray,ieletyp][ieledof|ientry,iele] -> idof|inz
# out.L1[α  ][αder     ][idof] -> gradient     α∈λxua
# out.L2[α,β][αder,βder][inz ] -> Hessian      α∈λxua, β∈λxua
const λxua   = 1:4
const λxu    = 1:3
const xua    = 2:4
const xu     = 2:3
const ind    = (Λ=1,X=2,U=3,A=4)
const nclass = length(ind) 

## Assembly of sparse
arrnum(α  )  =          α
arrnum(α,β)  = nclass + β + nclass*(α-1) 
mutable struct AssemblyDirect{OX,OU,IA}  <:Assembly
    L1 :: Vector{Vector{𝕣1      }}    # L1[α  ][αder     ]  α∈ λ,x,u,a
    L2 :: Matrix{Matrix{Sparse𝕣2}}    # L2[α,β][αder,βder]
    matrices     :: 𝔹
end  
function prepare(::Type{AssemblyDirect{OX,OU,IA}},model,dis;Xwhite=false,XUindep=false,UAindep=false,XAindep=false,matrices=true) where{OX,OU,IA}
    dofgr    = (allΛdofs(model,dis),allXdofs(model,dis),allUdofs(model,dis),allAdofs(model,dis))
    ndof     = getndof.(dofgr)
    neletyp  = getneletyp(model)
    asm      = Matrix{𝕫2}(undef,nclass+nclass^2,neletyp)
    nder     = (1,OX+1,OU+1,IA)
    L1       = Vector{Vector{Vector{𝕣}}}(undef,4)
    for α∈λxua
        nα   = nder[α]
        av   = asmvec!(view(asm,arrnum(α),:),dofgr[α],dis)
        L1[α] = Vector{Vector{𝕣}}(undef,nα)
        for αder = 1:nα 
            L1[α][αder] = copy(av)
        end
    end
    L2    = Matrix{Matrix{Sparse𝕣2}}(undef,4,4)
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
        L2[α,β] = Matrix{Sparse𝕣2}(undef,nα,nβ)
        for αder=1:nα,βder=1:nβ
            L2[α,β][αder,βder] = copy(am)
        end
    end
    out      = AssemblyDirect{OX,OU,IA}(L1,L2,matrices)
    return out,asm,dofgr
end
function zero!(out::AssemblyDirect)
    for L1∈out.L1 
        for ℓ1∈L1
            zero!(ℓ1)
        end
    end

    if out.matrices
        for L2∈out.L2 
            for ℓ2∈L2
                zero!(ℓ2)
            end
        end
    end
end

function addin!(out::AssemblyDirect,asm,iele,scale,eleobj::Eleobj,Λ,X,U,A,t,SP,dbg) where{Eleobj} 
    addin!(out::AssemblyDirect,asm,iele,scale,eleobj,nosecondorder(Eleobj),Λ,X,U,A,t,SP,dbg)
end

function addin!(out::AssemblyDirect{OX,OU,IA},asm,iele,scale,eleobj::Eleobj,fastresidual::Val{true}, 
                                Λ::NTuple{1  ,SVector{Nx}},
                                X::NTuple{NDX,SVector{Nx}},
                                U::NTuple{NDU,SVector{Nu}},
                                A::           SVector{Na} ,t,SP,dbg) where{OX,OU,IA,NDX,NDU,Nx,Nu,Na,Eleobj} 
    @assert NDX==OX+1 @sprintf("got OX=%i and NDX=%i. Expected OX+1==NDX",OX,NDX)
    @assert NDX==OX+1 @sprintf("got OU=%i and NDU=%i. Expected OU+1==NDU",OU,NDU)
    ndof   = (Nx, Nx, Nu, Na)
    nder   = (1,OX+1,OU+1,IA)
    Npfast =      Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials
    Np     = Nx + Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials

    X∂ = ntuple(ider->SVector{Nx}(∂ℝ{1,Npfast}(X[ider][idof],   Nx*(ider-1)            +idof, scale.X[idof])   for idof=1:Nx),OX+1)
    U∂ = ntuple(ider->SVector{Nu}(∂ℝ{1,Npfast}(U[ider][idof],   Nx*(OX+1)  +Nu*(ider-1)+idof, scale.U[idof])   for idof=1:Nu),OU+1)
    if IA == 1
        A∂   =        SVector{Na}(∂ℝ{1,Npfast}(A[      idof],   Nx*(OX+1)  +Nu*(OU+1)  +idof, scale.A[idof])   for idof=1:Na)
        R,FB = residual(eleobj, X∂,U∂,A∂,t,SP,dbg)
    else
        R,FB = residual(eleobj, X∂,U∂,A ,t,SP,dbg)
    end        
    iλ   = 1:ndof[ind.Λ]
    Lλ   = out.L1[ind.Λ]
    add_value!(Lλ[1] ,asm[arrnum(ind.Λ)],iele,R,ia=iλ)
    if out.matrices
        pβ       = 0
        for β∈xua, j=1:nder[β]
            iβ   = pβ.+(1:ndof[β])
            pβ  += ndof[β]
            Lλβ  = out.L2[ind.Λ,β]
            Lβλ  = out.L2[β,ind.Λ]
            if j≤size(Lλβ,2) # ...but only add into existing matrices of L2, for better sparsity
                add_∂!{1           }(Lλβ[1,j],asm[arrnum(ind.Λ,β)],iele,R,ia=iλ,ida=iβ)
                add_∂!{1,:transpose}(Lβλ[j,1],asm[arrnum(β,ind.Λ)],iele,R,ia=iλ,ida=iβ)
            end
        end
    end 
end
function addin!(out::AssemblyDirect{OX,OU,IA},asm,iele,scale,eleobj::Eleobj,fastresidual::Val{false}, 
    Λ::NTuple{1  ,SVector{Nx}},
    X::NTuple{NDX,SVector{Nx}},
    U::NTuple{NDU,SVector{Nu}},
    A::           SVector{Na} ,t,SP,dbg) where{OX,OU,IA,NDX,NDU,Nx,Nu,Na,Eleobj} 

    @assert NDX==OX+1 @sprintf("got OX=%i and NDX=%i. Expected OX+1==NDX",OX,NDX)
    @assert NDX==OX+1 @sprintf("got OU=%i and NDU=%i. Expected OU+1==NDU",OU,NDU)
    ndof   = (Nx, Nx, Nu, Na)
    nder   = (1,OX+1,OU+1,IA)
    Npfast =      Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials
    Np     = Nx + Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials

    Λ∂ =              SVector{Nx}(∂²ℝ{1,Np}(Λ[1   ][idof],                           idof, scale.Λ[idof])   for idof=1:Nx)
    X∂ = ntuple(ider->SVector{Nx}(∂²ℝ{1,Np}(X[ider][idof],Nx+Nx*(ider-1)            +idof, scale.X[idof])   for idof=1:Nx),OX+1)
    U∂ = ntuple(ider->SVector{Nu}(∂²ℝ{1,Np}(U[ider][idof],Nx+Nx*(OX+1)  +Nu*(ider-1)+idof, scale.U[idof])   for idof=1:Nu),OU+1)
    if IA == 1
        A∂   =        SVector{Na}(∂²ℝ{1,Np}(A[      idof],Nx+Nx*(OX+1)  +Nu*(OU+1)  +idof, scale.A[idof])   for idof=1:Na)
        L,FB = getlagrangian(eleobj, Λ∂,X∂,U∂,A∂,t,SP,dbg)
    else
        L,FB = getlagrangian(eleobj, Λ∂,X∂,U∂,A ,t,SP,dbg)
    end
    ∇L           = ∂{2,Np}(L)
    pα           = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
    for α∈λxua, i=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
        iα       = pα.+(1:ndof[α])
        pα      += ndof[α]
        Lα = out.L1[α]
        if i≤size(Lα,1)  # ...but only add into existing vectors of L1, for speed
            add_value!(out.L1[α][i] ,asm[arrnum(α)],iele,∇L,ia=iα)
        end
        if out.matrices
            pβ       = 0
            for β∈λxua, j=1:nder[β]
                iβ   = pβ.+(1:ndof[β])
                pβ  += ndof[β]
                Lαβ = out.L2[α,β]
                if i≤size(Lαβ,1) && j≤size(Lαβ,2) # ...but only add into existing matrices of L2, for better sparsity
                    add_∂!{1}(out.L2[α,β][i,j],asm[arrnum(α,β)],iele,∇L,ia=iα,ida=iβ)
                end
            end
        end
    end
end

## Assembly of bigsparse
function makepattern(OX,OU,IA,nstep,out) 
    # Looking at all steps, class, order of fdiff and Δstep, for rows and columns: which blocks are actualy nz?
    # return a sparse matrix of sparse matrices
    nder     = (1,OX+1,OU+1)
    maxblock = 1 + nstep*90  
    αblk     = 𝕫1(undef,maxblock)
    βblk     = 𝕫1(undef,maxblock)
    nz       = Vector{Sparse𝕣2}(undef,maxblock)
    nblock   = 0
    for step = 1:nstep
        for     α∈λxu 
            for β∈λxu
                Lαβ = out.L2[α,β]
                for     αder = 1:size(Lαβ,1)
                    for βder = 1:size(Lαβ,2)
                        for     iα ∈ finitediff(αder-1,nstep,step)
                            for iβ ∈ finitediff(βder-1,nstep,step)
                                nblock += 1   
                                αblk[nblock]=3*(step+iα.Δs-1)+α
                                βblk[nblock]=3*(step+iβ.Δs-1)+β
                                nz[  nblock]=Lαβ[1,1]  
                            end
                        end
                    end 
                end
            end
        end
    end   

    if IA==1
        Ablk = 3*nstep+1
        nblock +=1
        αblk[nblock] = Ablk                      
        βblk[nblock] = Ablk                    
        nz[  nblock] = out.L2[ind.A,ind.A][1,1]
        for step = 1:nstep
            for     α∈λxu 
                # loop over derivatives and finitediff is optimized out, as time derivatives will only 
                # be added into superbloc already reached by non-derivatives. No, it's not a bug...
                if size(out.L2[ind.A,α],1)>0
                    nblock += 1
                    αblk[nblock] = Ablk                
                    βblk[nblock] = 3*(step-1)+α          
                    nz[  nblock] = out.L2[ind.A,α][1,1]
                    nblock += 1
                    αblk[nblock] = 3*(step-1)+α            
                    βblk[nblock] = Ablk                  
                    nz[  nblock] = out.L2[α,ind.A][1,1]  
                end
            end
        end
    end
    u    = unique(i->(αblk[i],βblk[i]),1:nblock)

    return sparse(αblk[u],βblk[u],nz[u])
end
function preparebig(OX,OU,IA,nstep,out) 
    # create an assembler and allocate for the big linear system
    pattern                  = makepattern(OX,OU,IA,nstep,out)
    Lvv,Lvvasm,Lvasm,Lvdis   = prepare(pattern)
    Lv                       = 𝕣1(undef,size(Lvv,1))
    return Lvv,Lv,Lvvasm,Lvasm,Lvdis
end
function assemblebig!(Lvv,Lv,Lvvasm,Lvasm,asm,model,dis,out::AssemblyDirect{OX,OU,IA},state,nstep,Δt,SP,dbg) where{OX,OU,IA}
    zero!(Lvv)
    zero!(Lv )
    for step = 1:nstep
        state[step].SP   = SP
        
        assemble!(out,asm,dis,model,state[step],(dbg...,asm=:assemblebig!,step=step))

        for β∈λxu
            Lβ = out.L1[β]
            for βder = 1:size(Lβ,1)
                s = Δt^(1-βder)
                for iβ ∈ finitediff(βder-1,nstep,step)  # TODO transpose or not? Potential BUG to be revealed when cost on time derivative of X or U
                    βblk = 3*(step+iβ.Δs-1)+β
                    addin!(Lvasm,Lv ,Lβ[βder],βblk,iβ.w*s) 
                end
            end
        end
        for     α∈λxu 
            for β∈λxu
                Lαβ = out.L2[α,β]
                for     αder = 1:size(Lαβ,1)
                    for βder = 1:size(Lαβ,2)
                        s = Δt^(2-αder-βder)
                        for     iα ∈ finitediff(αder-1,nstep,step) # No transposition here, that's thoroughly checked against decay.
                            for iβ ∈ finitediff(βder-1,nstep,step) # No transposition here, that's thoroughly checked against decay.
                                αblk = 3*(step+iα.Δs-1)+α
                                βblk = 3*(step+iβ.Δs-1)+β
                                addin!(Lvvasm,Lvv,Lαβ[αder,βder],αblk,βblk,iα.w*iβ.w*s) 
                            end
                        end
                    end 
                end
            end
        end
        if IA==1
            Ablk = 3*nstep+1   
            addin!(Lvasm ,Lv ,out.L1[ind.A      ][1  ],Ablk     )
            addin!(Lvvasm,Lvv,out.L2[ind.A,ind.A][1,1],Ablk,Ablk)
            for α∈λxu
                Lαa = out.L2[α    ,ind.A]
                Laα = out.L2[ind.A,α    ]
                for αder = 1:size(Lαa,1)  # size(Lαa,1)==size(Laα,2) because these are 2nd derivatives of L
                    s = Δt^(1-αder)
                    for iα ∈finitediff(αder-1,nstep,step) # TODO transpose or not? BUG to be revealed when cost on time derivative sof X or U
                        αblk = 3*(step+iα.Δs-1)+α
                        addin!(Lvvasm,Lvv,Lαa[αder,1   ],αblk,Ablk,iα.w*s) 
                        addin!(Lvvasm,Lvv,Laα[1   ,αder],Ablk,αblk,iα.w*s) 
                    end
                end
            end
        end
    end   
end
function decrementbig!(state,Δ²,Lvasm,dofgr,Δv,nder,Δt,nstep) 
    Δ²                  .= 0.
    for (step,stateᵢ)    ∈ enumerate(state)
        for β            ∈ λxu
            for βder     = 1:nder[β]
                s        = Δt^(1-βder)
                for iβ   ∈ finitediff(βder-1,nstep,step)
                    βblk = 3*(step+iβ.Δs-1)+β   
                    Δβ   = disblock(Lvasm,Δv,βblk)
                    decrement!(stateᵢ,βder,Δβ.*iβ.w*s,dofgr[β])
                    if βder==1 
                        Δ²[β] = max(Δ²[β],sum(Δβ.^2)) 
                    end
                end
            end
        end
    end    
    if nder[4]==1
        Δa               = disblock(Lvasm,Δv,3*nstep+1)
        Δ²[ind.A]        = sum(Δa.^2)
        decrement!(state[1],1,Δa,dofgr[ind.A]) # all states share same A, so decrement only once
    end
end

"""
	DirectXUA{OX,OU,IA}

A non-linear direct solver for optimisation FEM.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
stateXUA        = solve(DirectXUA{OX,OU,IA};initialstate,time=0:1.:5)
```

The solver does not yet support interior point methods. 

# Parameters
- `OX`                0 for static analysis
                      1 for first order problems in time (viscosity, friction, measurement of velocity)
                      2 for second order problems in time (inertia, measurement of acceleration) 
- `OU`                0 for white noise prior to the unknown load process
                      2 otherwise
- `IA`                0 for XU problems (variables of class A will be unchanged)
                      1 for XUA problems                                                  

# Named arguments
- `dbg=(;)`           a named tuple to trace the call tree (for debugging).
- `verbose=true`      set to false to suppress printed output (for testing).
- `silenterror=false` set to true to suppress print out of error (for testing) .
- `initialstate`      a `State`.
- `time`              an `AbstractRange` of times at which to compute the steps.  Example: 0:0.1:5.                       
- `maxiter=50`        maximum number of Newton-Raphson iterations. 
- `maxΔλ=1e-5`        convergence criteria: a norm of the scaled `Λ` increment.
- `maxΔx=1e-5`        convergence criteria: a norm of the scaled `X` increment. 
- `maxΔu=1e-5`        convergence criteria: a norm of the scaled `U` increment. 
- `maxΔa=1e-5`        convergence criteria: a norm of the scaled `A` increment.
- `saveiter=false`    set to true so that the output `state` is a vector (over the iterations) of 
                      vectors (over the steps) of `State`s of the model (for debugging 
                      non-convergence). 
Setting the following flags to `true` will improve the sparsity of the system. But setting
a flag to `true` when the condition isn't met causes the Hessian to be wrong, which is detrimental for convergence.                      
- `Xwhite=false`      `true` if response measurement error is a white noise process.
- `XUindep=false`     `true` if response measurement error is independant of `U`
- `UAindep=false`     `true` if `U` is independant of `A`
- `XAindep=false`     `true` if response measurement error is independant of `A`

# Output

A vector of length equal to that of `time` containing the state of the optimized model at each of these steps.                       

See also: [`solve`](@ref), [`initialize!`](@ref), [`SweepX`](@ref), [`FreqXU`](@ref)
"""
struct DirectXUA{OX,OU,IA} <: AbstractSolver end 
function solve(::Type{DirectXUA{OX,OU,IA}},pstate,verbose::𝕓,dbg;
    time::AbstractRange{𝕣},
    initialstate::State,
    maxiter::ℤ=50,
    maxΔλ::ℝ=1e-5,maxΔx::ℝ=1e-5,maxΔu::ℝ=1e-5,maxΔa::ℝ=1e-5,
    saveiter::𝔹=false,
    kwargs...) where{OX,OU,IA}

    #  Mostly constants
    local LU
    nstep                 = length(time)
    Δt                    = (last(time)-first(time))/(nstep-1)
    γ                     = 0.
    nder                  = (1,OX+1,OU+1,IA)
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
    out,asm,dofgr         = prepare(AssemblyDirect{OX,OU,IA},model,dis;kwargs...)      # mem and assembler for system at any given step
    assemble!(out,asm,dis,model,state[1],(dbg...,solver=:DirectXUA,phase=:sparsity))     # create a sample "out" for preparebig
    Lvv,Lv,Lvvasm,Lvasm,Lvdis = preparebig(OX,OU,IA,nstep,out)                             # mem and assembler for big system

    for iter              = 1:maxiter
        verbose && @printf("\n    Iteration %3d\n",iter)

        verbose && @printf("        Assembling")
        SP = (γ=γ,iter=iter)
        assemblebig!(Lvv,Lv,Lvvasm,Lvasm,asm,model,dis,out,state,nstep,Δt,SP,(dbg...,solver=:DirectXUA,iter=iter))

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


