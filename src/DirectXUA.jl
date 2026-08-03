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
    _,_,∂A  = variate{P}(A)
#    ø   = nothing
#    C,_ = getlagrangian(eleobj,ø,ø,ø,∂A,ø,ø ,dbg)
    C,_ = getlagrangian(eleobj,∂A,nothing,dbg)
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
    if     mission==:matrices     P=1
    elseif mission==:vectors      P=0
    end
    if NDA == 1
        _,_,(∂X,∂U,∂A) = variate{1}((X,U,A),scale=(AllElements(scale.X),AllElements(scale.U),scale.A)) 
        R,FB     = getresidual(eleobj, ∂X,∂U,∂A,t,SP,dbg) 
    else
        _,_,(∂X,∂U)   = variate{1}((X,U ),scale=(AllElements(scale.X),AllElements(scale.U)))
        R,FB     = getresidual(eleobj, ∂X,∂U,  A,t,SP,dbg)
    end        
    ndof   = (Nx, Nx, Nu, Na)
    nder   = (NDΛ,NDX,NDU,NDA)
    iλ     = 1:ndof[ind.Λ]  ### TODO CHECK
    Lλ     = out.L1[ind.Λ]
    isassigned(Lλ,1) && add_value!(Lλ[1] ,asm[arrnum(ind.Λ)],iele,R,iλ;Δt)   
    if mission==:matrices
        pβ       = 0
        for β∈xua, βder=1:nder[β]
            iβ   = pβ.+(1:ndof[β])
            pβ  += ndof[β]
            Lλβ  = out.L2[ind.Λ,β]
            Lβλ  = out.L2[β,ind.Λ]
            isassigned(Lλβ,1,βder) && add_∂!{1                 }(Lλβ[1,βder],asm[arrnum(ind.Λ,β)],iele,R,iλ,iβ;Δt)
            isassigned(Lβλ,βder,1) && add_∂!{1,:plus,:transpose}(Lβλ[βder,1],asm[arrnum(β,ind.Λ)],iele,R,iλ,iβ;Δt)
        end
    end 
end
struct   DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA} end
function DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA}(out,asm,L,iele,Δt,class=λxua) where{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA}
    if     mission==:matrices     P=2
    elseif mission==:vectors      P=1
    end
    ndof         = (Nx, Nx, Nu, Na)
    nder         = (NDΛ,NDX,NDU,NDA)
    Np           = npartial(L)
    @assert Np   == Nx*NDΛ + Nx*NDX + Nu*NDU + Na*NDA # number of partials
    ∇L           = ∂{P,Np}(L)
    pα           = 0   # points into the partials, 1 entry before the start of relevant partial derivative in α,ider-loop
    for α∈class, αder=1:nder[α]   # we must loop over all time derivatives to correctly point into the adiff-partials...
        iα       = pα.+(1:ndof[α])
        pα      += ndof[α]
        Lα       = out.L1[α]
#        @show α,αder,iα,isassigned(Lα,αder)
        isassigned(Lα,αder) && add_value!(Lα[αder] ,asm[arrnum(α)],iele,∇L,iα;Δt)
        if mission==:matrices
            pβ       = 0
            for β∈class, βder=1:nder[β]
                iβ   = pβ.+(1:ndof[β])
                pβ  += ndof[β]
                Lαβ = out.L2[α,β]
#                @show α,β,αder,βder,iα,iβ,isassigned(Lαβ,αder,βder)
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
        _,_,(∂Λ,∂X,∂U,∂A) = variate{P}((Λ[1],X,U,A),scale=(scale.Λ,AllElements(scale.X),AllElements(scale.U),scale.A))
        L,FB              = getlagrangian(eleobj, ∂Λ,∂X,∂U,∂A,t,SP,dbg)
    else
        _,_,(∂Λ,∂X,∂U)  = variate{P}((Λ[1],X,U),scale=(scale.Λ,AllElements(scale.X),AllElements(scale.U)))
        L,FB            = getlagrangian(eleobj, ∂Λ,∂X,∂U,A  ,t,SP,dbg)
    end
    DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,NDΛ,NDX,NDU,NDA}(out,asm,L,iele,Δt)
end

function addin!{mission}(out::AssemblyDirect{NDΛ,NDX,NDU,NDA},asm,iele,scale,o::ElementCostAndConstraint{NSO,TargetElement,λxinod,λuinod},no_second_order::Val{true}, 
                           Λ::NTuple{1   ,SVector{Nx}},
                           X::NTuple{NDX_,SVector{Nx}},
                           U::NTuple{NDU_,SVector{Nu}},
                           A::            SVector{Na} ,t,Δt,SP,dbg) where{mission,NDΛ,NDX,NDU,NDA,NSO,TargetElement,λxinod,λuinod,NDX_,NDU_,Nx,Nu,Na} 

# Possible optimisation: adiff wrt eX,eU,A only, and create functionality to expand partials (for KKT operations with λX, λU)
# This affect the add_∂! calls
    local L
    iλX,ieX,iλU,ieU,iλA,ieA = dofpartition(typeof(o))
    ieΛ,Nλ                  = ieX,Nx                  
    ncstr                   = length(iλX)+length(iλU) 
    if     mission==:matrices     P=2
    elseif mission==:vectors      P=1
    end
    if     NDA == 1  
        _,N∂,(∂X,∂U,∂A)  = variate{P-1}((X,U,A),scale=(AllElements(scale.X),AllElements(scale.U),scale.A))
        ∂eX              = getsomedofs(∂X,ieX) 
        ∂eU              = getsomedofs(∂U,ieU) 
        eR,FB,eleres     = getresidual(o.eleobj, ∂eX,∂eU,∂A,t,SP,(dbg...,via=:ElementDofAndCoonstraintAccelerator),o.req.eleres)  
        class            = xua
    elseif NDA == 0
        _,N∂,(∂X,∂U)     = variate{P-1}((X,U  ),scale=(AllElements(scale.X),AllElements(scale.U)        ))
        ∂eX              = getsomedofs(∂X,ieX) 
        ∂eU              = getsomedofs(∂U,ieU) 
        eR,FB,eleres     = getresidual(o.eleobj, ∂eX,∂eU, A,t,SP,(dbg...,via=:ElementDofAndConstraintAccelerator),o.req.eleres)  
        class            = xu
    end
    Lλ                   = out.L1[ind.Λ]
    # out[asm[iasm,iele]] += R[ia]
    isassigned(Lλ,1) && add_value!(Lλ[1] ,asm[arrnum(ind.Λ)],iele,eR;iasm=ieΛ,Δt) #  I have a vector R from o.eleobj, and I assemble it at the end of o's vector
    if mission==:matrices # K,C,M
        ndof             = (Nλ , Nx, Nu, Na)
        iedof            = (ieΛ,ieX,ieU,ieA) # where in target's 2nd dim to put
        nder             = (NDΛ,NDX,NDU,NDA)
        iα               = SVector{length(ieΛ)}(1:length(ieΛ)) # and NOT ieα: where in eR (not R) to pick
        pβ               = 0
        for β∈class, βder=1:nder[β]
            ieβ          = pβ.+iedof[β]   # which of eR's partials to pick
            pβ          += ndof[β]
            Lλβ          = out.L2[ind.Λ,β]
            Lβλ          = out.L2[β,ind.Λ]
            # ∀i,j out[asm[ieΛᵢ+Nλ*(iedof[β]-1),iele]] += R[iαᵢ].dx[ieβⱼ]
            isassigned(Lλβ,1,βder) && add_∂!{1,:plus,:notranspose}(Lλβ[1,βder],asm[arrnum(ind.Λ,β)],iele,eR,iα,ieβ;nasm=Nλ,ndasm=ndof[β],iasm=iedof[ind.Λ],idasm=iedof[β],Δt)  
            isassigned(Lβλ,βder,1) && add_∂!{1,:plus,  :transpose}(Lβλ[βder,1],asm[arrnum(β,ind.Λ)],iele,eR,iα,ieβ;nasm=Nλ,ndasm=ndof[β],iasm=iedof[ind.Λ],idasm=iedof[β],Δt)  
        end
    end
    Releres            = revariate{P}(eleres) 
    if typeof(o.cost)  ≠ Nothing  
        Rcost          = o.cost(Releres,t,o.costargs...)
        cost           = chainrule(Rcost,to_order{P,N∂}(eleres))
        DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,0,NDX,NDU,NDA}(out,asm,cost,iele,Δt)
    end
    if typeof(o.gap)   ≠ Nothing   
        γ              = default{:γ}(SP,0.)
        λ_             = SVector(∂0(∂X)[iλX]...,∂0(∂U)[iλU]...) 
        λ              = to_order{P,N∂}(λ_) 
        Rgap           = o.gap(Releres,t,o.gapargs...)           #             2nd order derivative of gap wrt eleres
        gap            = chainrule(Rgap,to_order{P,N∂}(eleres))  # approximate 2nd order derivative of gap wrt ∂X,∂U,∂A
        mode           = o.mode(t)
        gap  isa SVector{ncstr} || error( "Functor `gap` must return a SVector{length(λxinod)+length(λuinod)}")
        mode isa SVector{ncstr} || error("Functor `mode` must return a SVector{length(λxinod)+length(λuinod)}")
        for iλ ∈ eachindex(λ)
            λᵢ,mᵢ,gapᵢ = λ[iλ],mode[iλ],gap[iλ]
            ΔL_=if   mᵢ==:equal;    gapᵢ * λᵢ     
            elseif   mᵢ==:positive; KKT(λᵢ,gapᵢ,γ) 
            elseif   mᵢ==:off;      0.5*(λᵢ * λᵢ) 
            end
            if iλ==1 L  = -ΔL_
            else     L -=  ΔL_  
            end
        end
        DirectXUA_lagrangian_addition!{mission,Nx,Nu,Na,0,NDX,NDU,NDA}(out,asm,L,iele,Δt)
    end   
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
primerstate    = initialize!(model)
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
- `silenterror=false` set to true to suppress print out of error (for testing).
- `primerstate`       an `AbstractVector` of `State` (1) or an `AbstractVector` of `Vector{State}` (2). `primerstate` must be with zero time derivatives.
                      It does not provide initial conditions for the problem, but an initial guess for the iterative solver.
                      As DirectXUA solves multiple experiments at once, the i-th element of `primerstate` is an initial guess for the i-th experiment. 
                      If provided with (1), the i-th initial guess is constructed internally with replicating the i-th input `State` for each timestep of the i-th timerange (see `time`).
                      If provided with (2), the i-th initial guess is the input `Vector{State}`, if the length corresponds to the i-th timerange (see `time`).
                      By convention the Adofs for building this primer are taken from the first experiment, first state.
- `time`              an `AbstractVector` (of same length as `primerstate`) of `AbstractRange` 
                      of times at which to compute the steps, one for each experiment.  Example: 0:0.1:5.                       
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
    primerstate::Union{Nothing,AbstractVector},
    maxiter::ℤ=50,
    maxΔλ::ℝ=1e-5,maxΔx::ℝ=1e-5,maxΔu::ℝ=1e-5,maxΔa::ℝ=1e-5,
    saveiter::𝔹=false,
    wanted::Wanted=Wanted{1,OX+1,OU+1,IA}(:all,:all)) where{OX,OU,IA,AR<:AbstractRange{𝕣},STATE<:State}

    #  Mostly constants
    local LU
    nexp,nstep,Δt         = length(time),length.(time),step.(time)
    γ                     = 0.
    nder                  = (1,OX+1,OU+1,IA)
    if IA==1  Δ², maxΔ²   = 𝕣1(undef,4), [maxΔλ^2,maxΔx^2,maxΔu^2,maxΔa^2] 
    else      Δ², maxΔ²   = 𝕣1(undef,3), [maxΔλ^2,maxΔx^2,maxΔu^2        ] 
    end
    
    # State storage
    S                     = State{1,OX+1,OU+1,@NamedTuple{γ::Float64,iter::Int64}}
    state                 = [Vector{S}(undef,nstep[iexp]) for iexp=1:nexp] # state[iexp][istep]
    if eltype(primerstate) <: Muscade.State
        model,dis,Ainit = primerstate[1].model, primerstate[1].dis, primerstate[1].A 
    else
        model,dis,Ainit = primerstate[1][1].model, primerstate[1][1].dis, primerstate[1][1].A
    end

    # Test input dimension compatibility
    length(primerstate)== nexp || muscadeerror("primerstate length must be equal to the number of experiments")

    for (iexp,primerstateᵢ) ∈ enumerate(primerstate)
        for (istep,timeᵢ) = enumerate(time[iexp])
            if eltype(primerstate) <: Muscade.State
                # Case where primerstate is an AbstractVector of State. A primer trajectory of state is created from a repetition of this state.
                primer = primerstateᵢ
            else
                # Case where primerstate is already an AbstractVector of Vector{State}.
                length(primerstateᵢ)== nstep[iexp] || muscadeerror("trajectories in `primerstate` vector must match their associated time grid in `time` vector")
                primer = primerstateᵢ[istep]
            end
            state[iexp][istep] = State{1,OX+1,OU+1}(timeᵢ,deepcopy(primer.Λ),deepcopy(primer.X),deepcopy(primer.U),Ainit,(γ=0.,iter=1),model,dis)
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