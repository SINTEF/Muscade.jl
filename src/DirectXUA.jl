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


mutable struct Wanted{NDΛ,NDX,NDU,NDA}  
    L1 :: Vector{Vector{𝕓}}    # L1[α  ][αder     ]  α∈ λ,x,u,a
    L2 :: Matrix{Matrix{𝕓}}    # L2[α,β][αder,βder]
end  
function wantL1(vectors::Symbol,nder)
    vectors==:all || muscadeerror("Invalid vector spec in Wanted")
    return [ [true           for i=1:nder[α]] for α∈λxua]
end
function wantL1(vectors::Tuple,nder)
    L1   = [ [false          for i=1:nder[α]] for α∈λxua]
    for s ∈ vectors
        L1[s[1]][s[2]] = true
    end
    return L1
end
function wantL2(matrices::Symbol,nder)
    matrices==:all || muscadeerror("Invalid matrix spec Wanted")
    return [ [~(α==β==ind.Λ) for i=1:nder[α],j=1:nder[β]] for α∈λxua, β∈λxua]
end
function wantL2(matrices::Tuple,nder)
    L2   = [ [false          for i=1:nder[α],j=1:nder[β]] for α∈λxua, β∈λxua]
    for s ∈ matrices
        L2[s[1],s[2]][s[3],s[4]] = true
    end
    return L2
end
function Wanted{NDΛ,NDX,NDU,NDA}(vectors::Union{Symbol,Tuple},matrices::Union{Symbol,Tuple}) where{NDΛ,NDX,NDU,NDA} 
    nder = (NDΛ,NDX,NDU,NDA)
    return Wanted{NDΛ,NDX,NDU,NDA}(wantL1(vectors,nder),wantL2(matrices,nder))
end
mutable struct AssemblyDirect{NDΛ,NDX,OU,NDA}  <:Assembly
    L1 :: Vector{Vector{𝕣1      }}    # L1[α  ][αder     ]  α∈ λ,x,u,a
    L2 :: Matrix{Matrix{Sparse𝕣2}}    # L2[α,β][αder,βder]
end  
function prepare(::Type{AssemblyDirect},model,dis,wanted::Wanted{NDΛ,NDX,NDU,NDA}) where{NDΛ,NDX,NDU,NDA}
    dofgr    = (allΛdofs(model,dis),allXdofs(model,dis),allUdofs(model,dis),allAdofs(model,dis))
    ndof     = getndof.(dofgr)
    neletyp  = getneletyp(model)
    asm      = Matrix{𝕫2}(undef,nclass+nclass^2,neletyp) 
    nder     = (NDΛ,NDX,NDU,NDA)
    L1       = Vector{Vector{𝕣1}}(undef,4)
    for α∈λxua
        nα   = nder[α]
        av   = asmvec!(view(asm,arrnum(α),:),dofgr[α],dis)
        L1[α] = Vector{𝕣1}(undef,nα)
        for αder = 1:nα 
            if wanted.L1[α][αder] 
                L1[α][αder] = copy(av)
            end
        end
    end
    L2    = Matrix{Matrix{Sparse𝕣2}}(undef,4,4)
    for α∈λxua, β∈λxua
        am    = asmmat!(view(asm,arrnum(α,β),:),view(asm,arrnum(α),:),view(asm,arrnum(β),:),ndof[α],ndof[β])
        nα,nβ = nder[α], nder[β]
        L2[α,β] = Matrix{Sparse𝕣2}(undef,nα,nβ)
        for αder=1:nα,βder=1:nβ
            if wanted.L2[α,β][αder,βder] 
                L2[α,β][αder,βder] = copy(am)
            end
        end
    end
    out      = AssemblyDirect{NDΛ,NDX,NDU,NDA}(L1,L2)
    return out,asm,dofgr
end
function zero!(out::AssemblyDirect) 
    for L1∈out.L1 
        for i ∈ eachindex(L1)
            isassigned(L1,i) && zero!(L1[i])
        end
    end
    for L2∈out.L2 
        for i ∈ eachindex(L2)
            isassigned(L2,i) && zero!(L2[i])
        end
    end
end

function addin!{mission}(out::AssemblyDirect,asm,iele,scale,eleobj::Acost,A::SVector{Na},dbg) where{Na,mission} # addin Atarget element
    if     mission==:matrices     P=2
    elseif mission==:vectors      P=1
    end
    ∂A  = revariate{P}(A)
    ø   = nothing
    C,_ = lagrangian(eleobj,ø,ø,ø,∂A,ø,ø ,dbg)
    ∇ₐC = ∂{P,Na}(C)
    isassigned(out.L1[ind.A],1) && add_value!(out.L1[ind.A][1],asm[arrnum(ind.A)],iele,∇ₐC)
    if mission==:matrices
        isassigned(out.L2[ind.A,ind.A],1,1) && add_∂!{1}(out.L2[ind.A,ind.A][1,1],asm[arrnum(ind.A,ind.A)],iele,∇ₐC)
    end
end
function addin!{mission}(out::AssemblyDirect,asm,iele,scale,eleobj::Eleobj,Λ,X,U,A,t,Δt,SP,dbg) where{Eleobj,mission} 
    addin!{mission}(out::AssemblyDirect,asm,iele,scale,eleobj,no_second_order(Eleobj),Λ,X,U,A,t,Δt,SP,dbg)
end
function addin!{mission}(out::AssemblyDirect{NDΛ,NDX,NDU,NDA},asm,iele,scale,eleobj::Eleobj,no_second_order::Val{true}, 
                                Λ::NTuple{1   ,SVector{Nx}},
                                X::NTuple{NDX_,SVector{Nx}},
                                U::NTuple{NDU_,SVector{Nu}},
                                A::            SVector{Na} ,t,Δt,SP,dbg) where{mission,NDΛ,NDX,NDU,NDA,NDX_,NDU_,Nx,Nu,Na,Eleobj} 
    ndof   = (Nx, Nx, Nu, Na)
    nder   = (NDΛ,NDX,NDU,NDA)
    if     mission==:matrices     P=1
    elseif mission==:vectors      P=0
    end
    if NDA == 1
        ∂X,∂U,∂A = revariate{1}((;X,U,A),(;X=scale.X,U=scale.U,A=scale.A)) # TODO as many partials as there are XUAs and their derivatives.  Possibly more than `out` wants
        R,FB     = residual(eleobj, ∂X,∂U,∂A,t,SP,dbg) # only elements with "residual" should set no_second_order=Val{true}
    else
        ∂X,∂U    = revariate{1}((;X,U  ),(;X=scale.X,U=scale.U))
        R,FB     = residual(eleobj, ∂X,∂U,  A,t,SP,dbg)
    end        
    iλ   = 1:ndof[ind.Λ]
    Lλ   = out.L1[ind.Λ]
    isassigned(Lλ,1) && add_value!(Lλ[1] ,asm[arrnum(ind.Λ)],iele,R,iλ;Δt)   
    if mission==:matrices
        pβ       = 0
        for β∈xua, βder=1:nder[β]
            iβ   = pβ.+(1:ndof[β])
            pβ  += ndof[β]
            Lλβ  = out.L2[ind.Λ,β]
            Lβλ  = out.L2[β,ind.Λ]
            isassigned(Lλβ,1,βder) && add_∂!{1                 }(Lλβ[1,βder],asm[arrnum(ind.Λ,β)],iele,R,iλ,iβ;Δt)
            isassigned(Lλβ,βder,1) && add_∂!{1,:plus,:transpose}(Lβλ[βder,1],asm[arrnum(β,ind.Λ)],iele,R,iλ,iβ;Δt)
        end
    end 
end
struct   DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA} end
function DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA}(out,asm,L,iele,Δt) where{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA}
    if     mission==:matrices     P=2
    elseif mission==:vectors      P=1
    end
    ndof         = (Nx, Nx, Nu, Na)
    nder         = (NDΛ,NDX,NDU,NDA)
    Np           = Nx + Nx*(NDX) + Nu*(NDU) + Na*NDA # number of partials
    λxua         = 1:4
    ∇L           = ∂{P,Np}(L)
    pα           = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
    for α∈λxua, αder=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
        iα       = pα.+(1:ndof[α])
        pα      += ndof[α]
        Lα       = out.L1[α]
        isassigned(Lα,αder) && add_value!(Lα[αder] ,asm[arrnum(α)],iele,∇L,iα;Δt)
        if mission==:matrices
            pβ       = 0
            for β∈λxua, βder=1:nder[β]
                iβ   = pβ.+(1:ndof[β])
                pβ  += ndof[β]
                Lαβ = out.L2[α,β]
                isassigned(Lαβ,αder,βder) && add_∂!{1}(Lαβ[αder,βder],asm[arrnum(α,β)],iele,∇L,iα,iβ;Δt)
            end
        end
    end
end    
function addin!{mission}(out::AssemblyDirect{NDΛ,NDX,NDU,NDA},asm,iele,scale,eleobj::Eleobj,no_second_order::Val{false}, 
                           Λ::NTuple{1   ,SVector{Nx}},
                           X::NTuple{NDX_,SVector{Nx}},
                           U::NTuple{NDU_,SVector{Nu}},
                           A::            SVector{Na} ,t,Δt,SP,dbg) where{mission,NDΛ,NDX,NDU,NDA,NDX_,NDU_,Nx,Nu,Na,Eleobj} 

    if     mission==:matrices     P=2
    elseif mission==:vectors      P=1
    end
    if NDA == 1   
        ∂Λ,∂X,∂U,∂A = revariate{P}((;Λ=Λ[1],X,U,A),(;Λ=scale.Λ,X=scale.X,U=scale.U,A=scale.A))
        L,FB        = getlagrangian(eleobj, ∂Λ,∂X,∂U,∂A,t,SP,dbg)
    else
        ∂Λ,∂X,∂U    = revariate{P}((;Λ=Λ[1],X,U),(;Λ=scale.Λ,X=scale.X,U=scale.U))
        L,FB        = getlagrangian(eleobj, ∂Λ,∂X,∂U,A  ,t,SP,dbg)
    end
    DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA}(out,asm,L,iele,Δt)
end
# Specialised to accelerate ElementCost and ElementConstraint
function addin!{mission}(out::AssemblyDirect{NDΛ,NDX,NDU,NDA},asm,iele,scale,eleobj::ElementCost,no_second_order::Val{false}, # TODO can we compact this?
                           Λ::NTuple{1   ,SVector{Nx}},
                           X::NTuple{NDX_,SVector{Nx}},
                           U::NTuple{NDU_,SVector{Nu}},
                           A::            SVector{Na} ,t,Δt,SP,dbg) where{mission,NDΛ,NDX,NDU,NDA,NDX_,NDU_,Nx,Nu,Na} 
         addin!{mission}(out,asm,iele,scale,eleobj,Val(true),Λ,X,U,A,t,Δt,SP,dbg) 
end
function addin!{mission}(out::AssemblyDirect{NDΛ,NDX,NDU,NDA},asm,iele,scale,eleobj::ElementCost,no_second_order::Val{true}, 
                           Λ::NTuple{1   ,SVector{Nx}},
                           X::NTuple{NDX_,SVector{Nx}},
                           U::NTuple{NDU_,SVector{Nu}},
                           A::            SVector{Na} ,t,Δt,SP,dbg) where{mission,NDΛ,NDX,NDU,NDA,NDX_,NDU_,Nx,Nu,Na} 
    if     mission==:matrices     P=2
    elseif mission==:vectors      P=1
    end
    if     NDA == 1  # NB: compile-time condition
        ∂X,∂U,∂A    = revariate{P-1}((;X,U,A),(;X=scale.X,U=scale.U,A=scale.A))
        R,FB,eleres = residual(eleobj.eleobj, ∂X,∂U,∂A,t,SP,(dbg...,via=:ElementCostAccelerator),eleobj.req.eleres)  
    elseif NDA == 0
        ∂X,∂U       = revariate{P-1}((;X,U ),(;X=scale.X,U=scale.U))
        R,FB,eleres = residual(eleobj.eleobj, ∂X,∂U,  A,t,SP,(dbg...,via=:ElementCostAccelerator),eleobj.req.eleres)  
    end
    Releres         = revariate{P}(eleres)
    
    Rcost           = eleobj.cost(Releres,t,eleobj.costargs...)
    cost            = chainrule(Rcost,to_order{P}(eleres))
    L               = Λ[1] ∘₁ R + cost
    DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA}(out,asm,L,iele,Δt)
end
function addin!{mission}(out::AssemblyDirect{NDΛ,NDX,NDU,NDA},asm,iele,scale,eleobj::ElementConstraint,no_second_order::Val{false}, 
                         Λ::NTuple{1   ,SVector{Nx}},
                         X::NTuple{NDX_,SVector{Nx}},
                         U::NTuple{NDU_,SVector{Nu}},
                         A::            SVector{Na} ,t,Δt,SP,dbg) where{mission,NDΛ,NDX,NDU,NDA,NDX_,NDU_,Nx,Nu,Na} 
         addin!{mission}(out,asm,iele,scale,eleobj,Val(true),Λ,X,U,A,t,Δt,SP,dbg) 
end
function addin!{mission}(out::AssemblyDirect{NDΛ,NDX,NDU,NDA},asm,iele,scale,eleobj::ElementConstraint,no_second_order::Val{true}, 
                           Λ::NTuple{1   ,SVector{Nx}},
                           X::NTuple{NDX_,SVector{Nx}},
                           U::NTuple{NDU_,SVector{Nu}},
                           A::            SVector{Na} ,t,Δt,SP,dbg) where{mission,NDΛ,NDX,NDU,NDA,NDX_,NDU_,Nx,Nu,Na} 
# TODO Specialised code to accelerate constraints in DirectXUA, but... it does not set FB, and DIrectXUA/solve has no line search...                                
    if     mission==:matrices     P=2
    elseif mission==:vectors      P=1
    end
    u               = getsomedofs(U,SVector{Nu}(1:Nu-1))
    λ               = ∂0(U)[Nu]
    γ               = default{:γ}(SP,0.)
    m               = eleobj.mode(t)
    if     NDA == 1  # NB: compile-time condition
        ∂X,∂U,∂A    = revariate{P-1}((X=X,U=U,A=A),(;X=scale.X,U=scale.U,A=scale.A))
        R,FB,eleres = residual(eleobj.eleobj, ∂X,∂U,∂A,t,SP,(dbg...,via=:ElementCoonstraintAccelerator),eleobj.req)  
    elseif NDA == 0
        ∂X,∂U       = revariate{P-1}((X=X,U=U    ),(;X=scale.X,U=scale.U))
        R,FB,eleres = residual(eleobj.eleobj, ∂X,∂U,  A,t,SP,(dbg...,via=:ElementConstraintAccelerator),eleobj.req)  
    end
    Releres         = revariate{P}(eleres)
    Rgap            = eleobj.gap(eleres,t,eleobj.gargs...)
    gap             = chainrule(Rgap,to_order{P}(eleres))
    L               = Λ[1] ∘₁ R +   if      m==:equal;    -gap*λ   
                                    elseif  m==:positive; -KKT(λ,gap,γ) 
                                    elseif  m==:off;      -0.5λ^2 
                                    end
    DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA}(out,asm,L,iele,Δt)
end


## Assembly of bigsparse for all time steps at once
function makepattern(NDA,nstep,out) 
    # Looking at all steps, class, order of fdiff and Δstep, for rows and columns: which blocks are actualy nz?
    # return a sparse matrix of sparse matrices
    maxblock = 1 + sum(nstep)*90  
    αblk     = 𝕫1(undef,maxblock)
    βblk     = 𝕫1(undef,maxblock)
    nz       = Vector{Sparse𝕣2}(undef,maxblock)
    nblock   = 0
    cumblk   = 0
    for iexp = 1:length(nstep)
        for istep = 1:nstep[iexp]
            for     α∈λxu 
                for β∈λxu
                    Lαβ = out.L2[α,β]
                    for     αder = 1:size(Lαβ,1)
                        for βder = 1:size(Lαβ,2)
                            if isassigned(Lαβ,αder,βder)
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
        end
        cumblk += 3*nstep[iexp]
    end   

    if NDA==1
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
                    Laα = out.L2[ind.A,α]
                    if size(Laα,1)>0 && isassigned(Laα,1)
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
function preparebig(NDA,nstep,out) 
    # create an assembler and allocate for the big linear system
    pattern                  = makepattern(NDA,nstep,out)
    # Muscade.spypattern(pattern)
    Lvv,Lvvasm,Lvasm,Lvdis   = prepare(pattern)
    Lv                       = 𝕣1(undef,size(Lvv,1))
    return Lvv,Lv,Lvvasm,Lvasm,Lvdis
end
struct assemblebig!{mission} end
function assemblebig!{mission}(Lvv,Lv,Lvvasm,Lvasm,asm,model,dis,out::AssemblyDirect{NDΛ,NDX,NDU,NDA},state,nstep,Δt,SP,dbg) where{mission,NDΛ,NDX,NDU,NDA}
    zero!(Lvv)
    zero!(Lv )
    cumblk = 0
    if NDA==1
        Ablk = 3sum(nstep)+1  
        assembleA!{mission}(out,asm,dis,model,state[1][1],(dbg...,asm=:assemblebig!)) 
        La, Laa = out.L1[ind.A], out.L2[ind.A,ind.A]
        isassigned(La ,1  ) && addin!(Lvasm ,Lv ,out.L1[ind.A      ][1  ],Ablk     )  
        isassigned(Laa,1,1) && addin!(Lvvasm,Lvv,out.L2[ind.A,ind.A][1,1],Ablk,Ablk)
    end    
    class = NDA==1 ? λxua : λxu
    for iexp = 1:length(nstep)
        for istep = 1:nstep[iexp]
            state[iexp][istep].SP   = SP
            assemble!{mission}(out,asm,dis,model,state[iexp][istep],idmult,(dbg...,asm=:assemblebig!,step=istep))
            for β∈class
                Lβ = out.L1[β]
                for βder = 1:size(Lβ,1)
                    if isassigned(Lβ,βder)
                        s    = Δt[iexp]^(1-βder)
                        for iβ ∈ finitediff(βder-1,nstep[iexp],istep)  # TODO transpose or not? Potential BUG to be revealed when cost on time derivative of X or U
                            βblk = β==ind.A ? Ablk : cumblk+3*(istep+iβ.Δs-1)+β
                            addin!(Lvasm,Lv ,Lβ[βder],βblk,iβ.w*s) 
                        end
                    end
                end
            end
            for α∈class, β∈class
                Lαβ = out.L2[α,β]
                for αder=1:size(Lαβ,1), βder=1:size(Lαβ,2)
                    if isassigned(Lαβ,αder,βder)
                        s    = Δt[iexp]^(2-αder-βder)
                        for iα∈finitediff(αder-1,nstep[iexp],istep), iβ∈finitediff(βder-1,nstep[iexp],istep) # No transposition here, that's thoroughly checked against decay.
                            αblk = α==ind.A ? Ablk : cumblk+3*(istep+iα.Δs-1)+α
                            βblk = β==ind.A ? Ablk : cumblk+3*(istep+iβ.Δs-1)+β
                            addin!(Lvvasm,Lvv,Lαβ[αder,βder],αblk,βblk,iα.w*iβ.w*s) 
                        end
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
    if nder[4]==1 # NDA==1
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
- `initialstate`      an `AbstractVector` of `State`: one initial state for each experiment.
                      `initialstate` must be with zero time derivatives.  It does not
                      provide initial conditions for the problem, but an initial guess
                      for the iterative solver.        
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
    wanted::Wanted=Wanted{1,OX+1,OU+1,IA}(:all,:all)) where{OX,OU,IA,AR<:AbstractRange{𝕣},STATE<:State}

    #  Mostly constants
    local LU
    nexp,nstep,Δt         = length(time),length.(time),step.(time)
    length(initialstate)== nexp || muscadeerror("initialstate and time must be of the same length=number of experiments") 
    γ                     = 0.
    nder                  = (1,OX+1,OU+1,IA)
    model,dis             = initialstate[1].model, initialstate[1].dis
    if IA==1  Δ², maxΔ²   = 𝕣1(undef,4), [maxΔλ^2,maxΔx^2,maxΔu^2,maxΔa^2] 
    else      Δ², maxΔ²   = 𝕣1(undef,3), [maxΔλ^2,maxΔx^2,maxΔu^2        ] 
    end

    # State storage
    S                     = State{1,OX+1,OU+1,@NamedTuple{γ::Float64,iter::Int64}}
    state                 = [Vector{S}(undef,nstep[iexp]) for iexp=1:nexp] # state[iexp][istep]
    s                     = initialstate[1]
    for (iexp,initialstateᵢ) ∈ enumerate(initialstate)
        for (istep,timeᵢ) = enumerate(time[iexp])
            state[iexp][istep] = State{1,OX+1,OU+1}(timeᵢ,deepcopy(initialstateᵢ.Λ),deepcopy(initialstateᵢ.X),deepcopy(initialstateᵢ.U),s.A,(γ=0.,iter=1),s.model,s.dis) # all state[iexp][istep].A are === 
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
    out,asm,dofgr         = prepare(AssemblyDirect,model,dis,wanted)          # mem and assembler for system at any given step
    assemble!{:matrices}(out,asm,dis,model,state[1][1],idmult,(dbg...,solver=:DirectXUA,phase=:sparsity))    # create a sample "out" for preparebig
    Lvv,Lv,Lvvasm,Lvasm,Lvdis = preparebig(IA,nstep,out)                                   # mem and assembler for big system
    cLvv                  = copy(Lvv)
    for iter              = 1:maxiter
        verbose && @printf("\n    Iteration %3d\n",iter)

        verbose && @printf("        Assembling")
        SP = (γ=γ,iter=iter)
        assemblebig!{:matrices}(Lvv,Lv,Lvvasm,Lvasm,asm,model,dis,out,state,nstep,Δt,SP,(dbg...,solver=:DirectXUA,iter=iter))
        sparser!(cLvv,Lvv,1e-20) # TODO user defined parameter
        verbose && @printf(", solving")
        try 
            LU = lu(cLvv) 
        catch 
            verbose && @printf("\n")
            muscadeerror("Lvv matrix factorization failed");
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


