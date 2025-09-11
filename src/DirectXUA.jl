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

## Assembly of sparse for a single time step
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


function addin!(out::AssemblyDirect,asm,iele,scale,eleobj::Acost,A::SVector{Na},dbg) where{Na} # addin Atarget element
    A∂  = SVector{Na,∂ℝ{2,Na,∂ℝ{1,Na,𝕣}}}(∂²ℝ{1,Na}(A[idof],idof, scale.A[idof])   for idof=1:Na)
    ø   = nothing
    C,_ = lagrangian(eleobj,ø,ø,ø,A∂,ø,ø ,dbg)
    ∇ₐC = ∂{2,Na}(C)
    add_value!(out.L1[ind.A][1],asm[arrnum(ind.A)],iele,∇ₐC)
    if out.matrices
        add_∂!{1}(out.L2[ind.A,ind.A][1,1],asm[arrnum(ind.A,ind.A)],iele,∇ₐC)
    end
end
addin!(out::AssemblyDirect,asm,iele,scale,eleobj::Eleobj,Λ,X,U,A,t,SP,dbg) where{Eleobj} =
    addin!(out::AssemblyDirect,asm,iele,scale,eleobj,no_second_order(Eleobj),Λ,X,U,A,t,SP,dbg)
function addin!(out::AssemblyDirect{OX,OU,IA},asm,iele,scale,eleobj::Eleobj,fastresidual::Val{true}, 
                                Λ::NTuple{1  ,SVector{Nx}},
                                X::NTuple{NDX,SVector{Nx}},
                                U::NTuple{NDU,SVector{Nu}},
                                A::           SVector{Na} ,t,SP,dbg) where{OX,OU,IA,NDX,NDU,Nx,Nu,Na,Eleobj} 
    @assert NDX==OX+1 @sprintf("got OX=%i and NDX=%i. Expected OX+1==NDX",OX,NDX)
    @assert NDU==OU+1 @sprintf("got OU=%i and NDU=%i. Expected OU+1==NDU",OU,NDU)
    ndof   = (Nx, Nx, Nu, Na)
    nder   = (1,OX+1,OU+1,IA)
    Npfast =      Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials  
    Np     = Nx + Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials  

    T  = ∂ℝ{1,Npfast,𝕣} 
    X∂ = ntuple(ider->SVector{Nx,T}(∂ℝ{1,Npfast}(X[ider][idof],   Nx*(ider-1)            +idof, scale.X[idof])   for idof=1:Nx), Val(NDX))
    U∂ = ntuple(ider->SVector{Nu,T}(∂ℝ{1,Npfast}(U[ider][idof],   Nx*(OX+1)  +Nu*(ider-1)+idof, scale.U[idof])   for idof=1:Nu), Val(NDU))
    if IA == 1
        A∂   =        SVector{Na,T}(∂ℝ{1,Npfast}(A[      idof],   Nx*(OX+1)  +Nu*(OU+1)  +idof, scale.A[idof])   for idof=1:Na)
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
    @assert NDU==OU+1 @sprintf("got OU=%i and NDU=%i. Expected OU+1==NDU",OU,NDU)
    ndof   = (Nx, Nx, Nu, Na)
    nder   = (1,OX+1,OU+1,IA)
    Npfast =      Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials
    Np     = Nx + Nx*(OX+1) + Nu*(OU+1) + Na*IA # number of partials

    T  = ∂ℝ{2,Np,∂ℝ{1,Np,𝕣}}
    Λ∂ =              SVector{Nx,T}(∂²ℝ{1,Np}(Λ[1   ][idof],                           idof, scale.Λ[idof])   for idof=1:Nx)
    X∂ = ntuple(ider->SVector{Nx,T}(∂²ℝ{1,Np}(X[ider][idof],Nx            +Nx*(ider-1)+idof, scale.X[idof])   for idof=1:Nx),Val(NDX))
    U∂ = ntuple(ider->SVector{Nu,T}(∂²ℝ{1,Np}(U[ider][idof],Nx+Nx*(OX+1)  +Nu*(ider-1)+idof, scale.U[idof])   for idof=1:Nu),Val(NDU))
    if IA == 1
        A∂   =        SVector{Na,T}(∂²ℝ{1,Np}(A[      idof],Nx+Nx*(OX+1)  +Nu*(OU+1)  +idof, scale.A[idof])   for idof=1:Na)
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
#            add_value!(out.L1[α][i] ,asm[arrnum(α)],iele,∇L,ia=iα)
            add_value!(Lα[i] ,asm[arrnum(α)],iele,∇L,ia=iα)
        end
        if out.matrices
            pβ       = 0
            for β∈λxua, j=1:nder[β]
                iβ   = pβ.+(1:ndof[β])
                pβ  += ndof[β]
                Lαβ = out.L2[α,β]
                if i≤size(Lαβ,1) && j≤size(Lαβ,2) # ...but only add into existing matrices of L2, for better sparsity
#                    add_∂!{1}(out.L2[α,β][i,j],asm[arrnum(α,β)],iele,∇L,ia=iα,ida=iβ)
                    add_∂!{1}(Lαβ[i,j],asm[arrnum(α,β)],iele,∇L,ia=iα,ida=iβ)
                end
            end
        end
    end
end

## Assembly of bigsparse for all time steps at once
function makepattern(IA,nstep,out) 
    # Looking at all steps, class, order of fdiff and Δstep, for rows and columns: which blocks are actualy nz?
    # return a sparse matrix of sparse matrices
    maxblock = 1 + sum(nstep)*90  
    αblk     = 𝕫1(undef,maxblock)
    βblk     = 𝕫1(undef,maxblock)
    nz       = Vector{Sparse𝕣2}(undef,maxblock)
    nblock   = 0
    cumblk      = 0
    for iexp = 1:length(nstep)
        for istep = 1:nstep[iexp]
            for     α∈λxu 
                for β∈λxu
                    Lαβ = out.L2[α,β]
                    for     αder = 1:size(Lαβ,1)
                        for βder = 1:size(Lαβ,2)
                            for     iα ∈ finitediff(αder-1,nstep[iexp],istep)
                                for iβ ∈ finitediff(βder-1,nstep[iexp],istep)
                                    nblock += 1   
                                    αblk[nblock]=cumblk+3*(istep+iα.Δs-1)+α
                                    βblk[nblock]=cumblk+3*(istep+iβ.Δs-1)+β
                                    nz[  nblock]=Lαβ[1,1]  
                                end
                            end
                        end 
                    end
                end
            end
        end
        cumblk += 3*nstep[iexp]
    end   

    if IA==1
        Ablk = 3*sum(nstep)+1
        nblock +=1
        αblk[nblock] = Ablk                      
        βblk[nblock] = Ablk                    
        nz[  nblock] = out.L2[ind.A,ind.A][1,1]
        cumblk = 0
        for iexp     = 1:length(nstep)
            for istep = 1:nstep[iexp]
                for α∈λxu 
                    # loop over derivatives and finitediff is optimized out, as time derivatives will only 
                    # be added into superbloc already reached by non-derivatives. No, it's not a bug...
                    if size(out.L2[ind.A,α],1)>0
                        nblock += 1
                        αblk[nblock] = Ablk                
                        βblk[nblock] = cumblk+3*(istep-1)+α          
                        nz[  nblock] = out.L2[ind.A,α][1,1]
                        nblock += 1
                        αblk[nblock] = cumblk+3*(istep-1)+α            
                        βblk[nblock] = Ablk                  
                        nz[  nblock] = out.L2[α,ind.A][1,1]  
                    end
                end
            end
            cumblk += 3*nstep[iexp]
        end
    end
    u    = unique(i->(αblk[i],βblk[i]),1:nblock)

    return sparse(αblk[u],βblk[u],nz[u])
end
function preparebig(IA,nstep,out) 
    # create an assembler and allocate for the big linear system
    pattern                  = makepattern(IA,nstep,out)
    # Muscade.spypattern(pattern)
    Lvv,Lvvasm,Lvasm,Lvdis   = prepare(pattern)
    Lv                       = 𝕣1(undef,size(Lvv,1))
    return Lvv,Lv,Lvvasm,Lvasm,Lvdis
end
function assemblebig!(Lvv,Lv,Lvvasm,Lvasm,asm,model,dis,out::AssemblyDirect{OX,OU,IA},state,nstep,Δt,SP,dbg) where{OX,OU,IA}
    zero!(Lvv)
    zero!(Lv )
    cumblk = 0
    if IA==1
        Ablk = 3sum(nstep)+1  
        assembleA!(out,asm,dis,model,state[1][1],(dbg...,asm=:assemblebig!)) 
        addin!(Lvasm ,Lv ,out.L1[ind.A      ][1  ],Ablk     )  
        addin!(Lvvasm,Lvv,out.L2[ind.A,ind.A][1,1],Ablk,Ablk)
    end    
    class = IA==1 ? λxua : λxu
    for iexp = 1:length(nstep)
        for istep = 1:nstep[iexp]
            state[iexp][istep].SP   = SP
            assemble!(out,asm,dis,model,state[iexp][istep],(dbg...,asm=:assemblebig!,step=istep))
            for β∈class
                Lβ = out.L1[β]
                for βder = 1:size(Lβ,1)
                    s    = Δt[iexp]^(1-βder)
                    for iβ ∈ finitediff(βder-1,nstep[iexp],istep)  # TODO transpose or not? Potential BUG to be revealed when cost on time derivative of X or U
                        βblk = β==ind.A ? Ablk : cumblk+3*(istep+iβ.Δs-1)+β
                        addin!(Lvasm,Lv ,Lβ[βder],βblk,iβ.w*s) 
                    end
                end
            end
            for α∈class, β∈class
                Lαβ = out.L2[α,β]
                for αder=1:size(Lαβ,1), βder=1:size(Lαβ,2)
                    s    = Δt[iexp]^(2-αder-βder)
                    for iα∈finitediff(αder-1,nstep[iexp],istep), iβ∈finitediff(βder-1,nstep[iexp],istep) # No transposition here, that's thoroughly checked against decay.
                        αblk = α==ind.A ? Ablk : cumblk+3*(istep+iα.Δs-1)+α
                        βblk = β==ind.A ? Ablk : cumblk+3*(istep+iβ.Δs-1)+β
                        addin!(Lvvasm,Lvv,Lαβ[αder,βder],αblk,βblk,iα.w*iβ.w*s) 
                    end
                end
            end
        end 
        cumblk += 3*nstep[iexp] 
    end 
end
function decrementbig!(state,Δ²,Lvasm,dofgr,Δv,nder,Δt,nstep) 
    Δ²                      .= 0.
    cumblk                   = 0
    for iexp                 = 1:length(nstep)
        for istep            = 1:nstep[iexp]    
            for β            ∈ λxu
                for βder     = 1:nder[β]
                    s        = Δt[iexp]^(1-βder)
                    for iβ   ∈ finitediff(βder-1,nstep[iexp],istep)
                        βblk = cumblk+3*(istep+iβ.Δs-1)+β   
                        Δβ   = disblock(Lvasm,Δv,βblk)
                        decrement!(state[iexp][istep],βder,Δβ.*iβ.w*s,dofgr[β])
                        if βder==1 
                            Δ²[β] = max(Δ²[β],sum(Δβ.^2)) 
                        end
                    end
                end
            end
        end
        cumblk += 3*nstep[iexp]
    end    
    if nder[4]==1 # IA==1
        Δa               = disblock(Lvasm,Δv,3*sum(nstep)+1)
        Δ²[ind.A]        = sum(Δa.^2)
        decrement!(state[1][1],1,Δa,dofgr[ind.A]) # all states share same A, so decrement only once
    end
end



"""
	DirectXUA{OX,OU,IA}

A non-linear direct solver for optimisation FEM.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
stateXUA        = solve(DirectXUA{OX,OU,IA};initialstate=[initialstate],time=[0:1.:5])
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
- `initialstate`      an `AbstractVector` of `State`: one initial state for each experiment
- `time`              an `AbstractVector` (of same length as `initialstate`) of `AbstractRange` 
                      of times at which to compute the steps.  Example: 0:0.1:5.                       
- `maxiter=50`        maximum number of Newton-Raphson iterations. 
- `maxΔλ=1e-5`        convergence criteria: a norm of the scaled `Λ` increment.
- `maxΔx=1e-5`        convergence criteria: a norm of the scaled `X` increment. 
- `maxΔu=1e-5`        convergence criteria: a norm of the scaled `U` increment. 
- `maxΔa=1e-5`        convergence criteria: a norm of the scaled `A` increment.
- `saveiter=false`    set to true so that the output `state` contains the states 
                      at each Newton-Raphson iteration (for debugging 
                      non-convergence). 
Setting the following flags to `true` will improve the sparsity of the system. But setting
a flag to `true` when the condition isn't met causes the Hessian to be wrong, which is detrimental for convergence.                      
- `Xwhite=false`      `true` if response measurement error is a white noise process.
- `XUindep=false`     `true` if response measurement error is independant of `U`
- `UAindep=false`     `true` if `U` is independant of `A`
- `XAindep=false`     `true` if response measurement error is independant of `A`

# Output

- `state`, where `state[iexp][itime]` contains the state of the optimized model at each of these steps, or if `saveiter=true` then `state[iiter][iexp][itime]` is a state.

See also: [`solve`](@ref), [`initialize!`](@ref), [`SweepX`](@ref), [`FreqXU`](@ref)
"""
struct DirectXUA{OX,OU,IA} <: AbstractSolver end 
function solve(::Type{DirectXUA{OX,OU,IA}},pstate,verbose::𝕓,dbg;
    time::AbstractVector{AR},
    initialstate::AbstractVector{STATE},
    maxiter::ℤ=50,
    maxΔλ::ℝ=1e-5,maxΔx::ℝ=1e-5,maxΔu::ℝ=1e-5,maxΔa::ℝ=1e-5,
    saveiter::𝔹=false,
    kwargs...) where{OX,OU,IA,AR<:AbstractRange{𝕣},STATE<:State}

    #  Mostly constants
    local LU
    nexp,nstep,Δt         = length(time),length.(time),step.(time)
    length(initialstate)== nexp || muscadeerror("initialstate and time must be of the same length") 
    γ                     = 0.
    nder                  = (1,OX+1,OU+1,IA)
    model,dis             = initialstate[1].model, initialstate[1].dis
    if IA==1  Δ², maxΔ²   = 𝕣1(undef,4), [maxΔλ^2,maxΔx^2,maxΔu^2,maxΔa^2] 
    else      Δ², maxΔ²   = 𝕣1(undef,3), [maxΔλ^2,maxΔx^2,maxΔu^2        ] 
    end

    # State storage
    S                     = State{1,OX+1,OU+1,@NamedTuple{γ::Float64,iter::Int64}}
    state                 = [Vector{S}(undef,nstep[iexp]) for iexp=1:nexp] # state[iexp][istep]
    # for iexp              = 1:nexp
    #     s                 = State{1,OX+1,OU+1}(copy(initialstate[iexp],time=time[iexp][1],SP=(γ=0.,iter=1)))   
    #     for (istep,timeᵢ) = enumerate(time[iexp])
    #         state[iexp][istep] = istep==1 ? s : State(timeᵢ,deepcopy(s.Λ),deepcopy(s.X),deepcopy(s.U),s.A,s.SP,s.model,s.dis) # all state[iexp][istep].A are === 
    #     end
    # end
    s = initialstate[1]
    for (iexp,initialstateᵢ) ∈ enumerate(initialstate)
        for (istep,timeᵢ) = enumerate(time[iexp])
            state[iexp][istep] = State{1,OX+1,OU+1}(timeᵢ,deepcopy(initialstateᵢ.Λ),deepcopy(initialstateᵢ.X),deepcopy(initialstateᵢ.U),s.A,SP=(γ=0.,iter=1),s.model,s.dis) # all state[iexp][istep].A are === 
        end
    end

    if saveiter
        stateiter         = Vector{Vector{Vector{S}}}(undef,maxiter) # stateiter[iiter][iexp][istep] 
        pstate[]          = stateiter
    else
        pstate[]          = state                                                                            
    end    

    # Prepare assembler
    verbose && @printf("\n    Preparing assembler\n")
    out,asm,dofgr         = prepare(AssemblyDirect{OX,OU,IA},model,dis;kwargs...)          # mem and assembler for system at any given step
    assemble!(out,asm,dis,model,state[1][1],(dbg...,solver=:DirectXUA,phase=:sparsity))    # create a sample "out" for preparebig
    Lvv,Lv,Lvvasm,Lvasm,Lvdis = preparebig(IA,nstep,out)                                   # mem and assembler for big system

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
            stateiter[iter]     = deepcopy.(state) # deep, to avoid common A across iterations
        end
        verbose          && @printf(  "        maxₜ(|ΔΛ|)=%7.1e ≤ %7.1e  \n",√(Δ²[ind.Λ]),√(maxΔ²[ind.Λ]))
        verbose          && @printf(  "        maxₜ(|ΔX|)=%7.1e ≤ %7.1e  \n",√(Δ²[ind.X]),√(maxΔ²[ind.X]))
        verbose          && @printf(  "        maxₜ(|ΔU|)=%7.1e ≤ %7.1e  \n",√(Δ²[ind.U]),√(maxΔ²[ind.U]))
        verbose && IA==1 && @printf(  "             |ΔA| =%7.1e ≤ %7.1e  \n",√(Δ²[ind.A]),√(maxΔ²[ind.A]))
        if all(Δ².≤maxΔ²)  
            verbose      && @printf("\n    Converged in %3d iterations.\n",iter)
            verbose      && @printf(  "    nel=%d, nvar=%d, nstep=%d\n",getnele(model),length(Lv),sum(nstep))
            break#out of iter
        end
        iter<maxiter || muscadeerror(@sprintf("no convergence after %3d iterations. \n",iter))
    end # for iter
    return
end


