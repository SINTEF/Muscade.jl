# mder: the order+1 of time derivatives of Λ,X,U,A with respect to which the Lagrangian is to be differentiated in `out` 
# nder: the order+1 of time derivatives of Λ,X,U,A given to addin! 
#
# We make a distinction between nΛder==nAder==1, nXder=length(X), nUder=length(U) on the one hand, and m?der ≤ n?der.  This allows
# 1) to freeze A for XU algo (or any class)
# 2) not to compute cost on U′ or U′′ if these costs are known to be zero (same with X)                                      

# dis.dis[ieletyp].index[iele].X|U|A[ieledof]       - disassembling model state into element dofs
# dis.dis[ieletyp].scale.Λ|X|U|A[ieledof]           - scaling each element type 
# dis.scaleΛ|X|U|A[imoddof]                         - scaling the model state
# dis.field  X|U|A[imoddof]                         - field of dofs in model state
# asm1[iarray,ieletyp][ieledof|ientry,iele] -> idof|inz
# out1.L1[α  ][αder     ][idof] -> gradient     α∈λxua
# out1.L2[α,β][αder,βder][inz ] -> Hessian      α∈λxua, β∈λxua
const λxua   = 1:4
const λxu    = 1:3
const ind    = (Λ=1,X=2,U=3,A=4)
const nder   = 3
const nclass = 4 
const nvec   = nclass
const nmat   = nclass^2  # we leave undef subarrays in asm1 for unwanted derivatives.
arrnum(α  )  =        α
arrnum(α,β)  = nvec + β + nclass*(α-1) 
mutable struct AssemblyDirect{Mder,T1,T2}  <:Assembly
    L1 :: T1   
    L2 :: T2   
end  
function prepare(::Type{AssemblyDirect},model,dis,mder) 
    dofgr    = (allΛdofs(model,dis),allXdofs(model,dis),allUdofs(model,dis),allAdofs(model,dis))
    ndof     = getndof.(dofgr)
    neletyp  = getneletyp(model)
    asm      = Matrix{𝕫2}(undef,nvec+nmat,neletyp)
    L1       = [[asmvec!(view(asm,arrnum(α  ),:),dofgr[α],dis   )                                              for αder=1:mder[α]               ] for α∈λxua        ] # recomputes asm three mder  times
    L2       = [[asmmat!(view(asm,arrnum(α,β),:),view(asm,arrnum(α),:),view(asm,arrnum(β),:),ndof[α],ndof[β])  for αder=1:mder[α],βder=1:mder[β]] for α∈λxua, β∈λxua] # recomputes asm three mder² times
    out      = AssemblyDirect{mder,typeof(L1),typeof(L2)}(L1,L2)
    return out,asm
end
function zero!(out::AssemblyDirect)
    for α∈λxua 
        zero!.(out.L1[α])
        for β∈λxua
            zero!.(out.L2[α,β])
        end
    end
end
function addin!(out::AssemblyDirect{Mder,T1,T2},asm,iele,scale,eleobj,Λ::NTuple{nΛder,SVector{Nx}},
                                                                X::NTuple{nXder,SVector{Nx}},
                                                                U::NTuple{nUder,SVector{Nu}},
                                                                A::             SVector{Na}   ,t,SP,dbg) where{Mder,T1,T2,nΛder,nXder,nUder,Nx,Nu,Na} 
# asm[iarray         ][ieledof|ientry,iele] -> idof|inz
# mder: the derivatives wanted in out 
# nder: the derivatives given to addin! 
# αder ≤ nder & ider ≤ mder : take the time derivative and variate it
# nder < ider ≤ mder        : variate 0.  So a dynamic analysis from a static state will return zero inertial force, but non-zero mass matrix
# mder < ider ≤ nder        : do not pass to element.  So a static analysis starting from a dynamic state will return neither inertial forces nor mass matrix
    ndof  = (Nx  ,Nx   ,Nu   ,Na  )
    Nz    = Nx+Mder[2]*Nx+Mder[3]*Nu+Na

    Λ∂ =               SVector{Nx}(  ∂²ℝ{1,Nz}(Λ[1   ][idof],  idof)   for idof=1:Nx)
    X∂ = ntuple(Mder[2]) do ider 
        qx = Nx+(ider-1)*Nx
        if ider≤nXder SVector{Nx}(  ∂²ℝ{1,Nz}(X[ider][idof],qx+idof)   for idof=1:Nx)
        else          SVector{Nx}(  ∂²ℝ{1,Nz}(0.           ,qx+idof)   for idof=1:Nx)
        end
    end
    U∂ = ntuple(Mder[3]) do ider 
        qu = Nx+Mder[2]*Nx+(ider-1)*Nu
        if ider≤nUder SVector{Nu}(  ∂²ℝ{1,Nz}(U[ider][idof],qu+idof)   for idof=1:Nu)
        else          SVector{Nu}(  ∂²ℝ{1,Nz}(0.           ,qu+idof)   for idof=1:Nu)
        end
    end
    qa = Nx+Mder[2]*Nx+Mder[3]*Nu
    A∂ =              SVector{Na}(  ∂²ℝ{1,Nz}(A[      idof],qa+idof)   for idof=1:Na)

    L,FB         = getlagrangian(eleobj, Λ∂,X∂,U∂,A∂,t,SP,dbg)
 
    ∇L           = ∂{2,Nz}(L)
    pα           = 0   # point 1 under the start of relevant partial derivative in α,ider-loop
    for α∈λxua, i=1:Mder[α]
        iα       = pα.+(1:ndof[α])
        pα      += ndof[α]
        add_value!(out.L1[α][i] ,asm[arrnum(α)],iele,∇L,iα)
        pβ       = 0
        for β∈λxua, j=1:Mder[β]
            iβ   = pβ.+(1:ndof[β])
            pβ  += ndof[β]
            add_∂!{1}( out.L2[α,β][i,j],asm[arrnum(α,β)],iele,∇L,iα,iβ)
        end
    end
end

######################

mutable struct AssemblyDirectLine  <:Assembly
    ming  :: 𝕣
    minλ  :: 𝕣
    Σλg   :: 𝕣
    npos  :: 𝕫
end  
struct AssemblerDirectLine end
function prepare(::Type{AssemblyDirectLine},model) 
    neletyp  = getneletyp(model)
    asm      = Matrix{𝕫2}(undef,0,neletyp)
    out      = AssemblyDirectLine(∞,∞,0.,0)
    return out,asm
end
function zero!(out::AssemblyDirectLine)
    out.ming = ∞    
    out.minλ = ∞
    out.Σλg  = 0.
    out.npos = 0    
end
function addin!(out::AssemblyDirectLine,asm,iele,scale,eleobj,Λ,X,U,A,t,SP,dbg) 
    L,FB    = getlagrangian(eleobj,∂0(Λ),X,U,A,t,SP,dbg)
    if hasfield(typeof(FB),:mode) && FB.mode==:positive
        out.ming   = min(out.ming,FB.g)
        out.minλ   = min(out.minλ,FB.λ)
        out.Σλg   += FB.g*FB.λ
        out.npos  += 1
    end
end

######################

function preparebig(ND,NA,nstep,out)
    istep,jstep              = FDsparsity(ND-1,nstep)
    ndiff                    = length(istep)                      # number of 3*3 superblocks in the XU part of pattern 
    nrow = ncol              = length(λxu)*nstep + NA             # number of rows and cols in pattern
    nΛXUblock                = ndiff*length(λxu)^2                  
    nAblock                  = NA*(2*(length(λxu)*nstep) + 1) 
    nblock                   = nΛXUblock + nAblock                 # nnz of pattern

    irow                     = 𝕫1(undef,nblock) 
    icol                     = 𝕫1(undef,nblock) 
    nz                       = Vector{SparseMatrixCSC{𝕣,𝕫}}(undef,nblock)
    iblock                   = 0
    for i ∈ 1:ndiff
        for α∈λxu, β∈λxu
            iblock          += 1
            irow[iblock]     = (istep[i]-1)*length(λxu) + α
            icol[iblock]     = (jstep[i]-1)*length(λxu) + β
            nz[  iblock]     = out.L2[α,β][1,1]
        end 
    end 
    if NA == 1
        for istep            = 1:nstep   
            for α∈λxu
                istep,α
                iblock      += 1
                irow[iblock] = (istep-1)*length(λxu) + α
                icol[iblock] = ncol
                nz[  iblock] = out.L2[α      ,ind[:A]][1,1]
                iblock += 1
                irow[iblock] = nrow
                icol[iblock] = (istep-1)*length(λxu) + α
                nz[  iblock] = out.L2[ind[:A],α     ][1,1]
            end
        end        
        iblock              += 1
        irow[iblock]         = nrow
        icol[iblock]         = ncol
        nz[  iblock]         = out.L2[ind[:A],ind[:A]][1,1]
    end
    pattern                  = sparse(irow,icol,nz)
    Lvv,bigasm               = prepare(pattern)
    Lv                       = 𝕣1(undef,size(Lvv,1))

    return Lv,Lvv,bigasm
end

function assemblebig!(Lvv,Lv,bigasm,asm,model,dis,out,state,nstep,Δt,NA,γ,dbg)
    zero!(Lvv)
    zero!(Lv )
    for step = 1:nstep
        state[step].SP   = (γ=γ ,)
        
        assemble!(out,asm,dis,model,state[step],(dbg...,asm=:assemblebig!,step=step))

        for β∈λxu
            for βder = 1:size(out.L1[β],1)
                for iβ ∈ finitediff(βder-1,nstep,step;transposed=true)
                    βblk = 3*(step+iβ.Δs-1)+β
                    addin!(bigasm,Lv,out.L1[β][βder],βblk,iβ.w/Δt)
                end
            end
        end
        for     α∈λxu 
            for β∈λxu
                for     αder = 1:size(out.L2[α,β],1)
                    for βder = 1:size(out.L2[α,β],2)
                        for     iα ∈ finitediff(αder-1,nstep,step;transposed=true)
                            for iβ ∈ finitediff(βder-1,nstep,step;transposed=true)
                                αblk = 3*(step+iα.Δs-1)+α
                                βblk = 3*(step+iβ.Δs-1)+β
                                if αblk==11 && βblk==2
                                    @show step,α,β,αder,βder,iα,iβ
                                    @show βder-1,nstep,step
                                end
                                addin!(bigasm,Lvv,out.L2[α,β][αder,βder],αblk,βblk,iα.w*iβ.w/Δt^2)
                            end
                        end
                    end 
                end
            end
        end
        if NA==1
            Ablk = length(λxu)*nstep+1
            addin!(bigasm,Lv ,out.L1[ind.A      ][1  ],Ablk     )
            addin!(bigasm,Lvv,out.L2[ind.A,ind.A][1,1],Ablk,Ablk)
            for α∈λxu
                for αder = 1:size(out.L2[α,ind.A],1)
                    for iα ∈finitediff(αder-1,nstep,step;transposed=true)
                        αblk = 3*(step+iα.Δs-1)+α
                        addin!(bigasm,Lvv,out.L2[α,ind.A][αder,1],αblk,Ablk,iα.w/Δt) # TTODO audit weight iα.w/Δt
                    end
                end
            end
            for β∈λxu
                for βder = 1:size(out.L2[ind.A,β],2)
                    for iβ ∈finitediff(βder-1,nstep,step;transposed=true)
                        βblk = 3*(step+iβ.Δs-1)+β
                        addin!(bigasm,Lvv,out.L2[ind.A,β][1,βder],Ablk,βblk,iβ.w/Δt)
                    end
                end
            end
        end
    end   
end


"""
	DirectXUA

A non-linear direct solver for optimisation FEM.

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
stateX          = solve(SweepX{0}  ;initialstate=initialstate,time=[0.,1.])
stateXUA        = solve(DirectXUA;initialstate=stateX)
```

The interior point algorithm requires a starting point that is
strictly primal feasible (at all steps, all inequality constraints must have 
positive gaps) and strictly dual feasible (at all steps, all associated Lagrange 
multipliers must be strictly positive). Note the use of `setdof!` in the example
above to ensure dual feasibility.

# Named arguments
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a vector of `state`s, one for each load case in the optimization problem, 
                      obtained from one or several previous `SweepX` analyses
- `maxiter=50`        maximum number of Newton-Raphson iterations 
- `maxΔa=1e-5`        "outer" convergence criteria: a norm on the scaled `A` increment 
- `maxΔy=1e-5`        "inner" convergence criteria: a norm on the scaled `Y=[ΛXU]` increment 
- `saveiter=false`    set to true so that the output `state` is a vector (over the Aiter) of 
                      vectors (over the steps) of `State`s of the model (for debugging 
                      non-convergence). 
- `maxLineIter=50`    maximum number of iterations in the linear search that ensure interior points   
- `β=0.5`             `β∈]0,1[`. In the line search, if conditions are not met, then a new line-iteration is done
                      with `s *= β` where  `β→0` is a hasty backtracking, while `β→1` stands its ground.            
- `γfac=0.5`          `γfac∈[0,1[`. At each iteration, the barrier parameter γ is taken as `γ = (∑ⁿᵢ₌₁ λᵢ gᵢ)/n*γfac` where
                      `(∑ⁿᵢ₌₁ λᵢ gᵢ)/n` is the complementary slackness, and `n` the number of inequality constraints.
- `γbot=1e-8`         `γ` will not be reduced to under the original complementary slackness divided by `γbot`,
                      to avoid conditioning problems.                                               

# Output

A vector of length equal to that of `initialstate` containing the state of the optimized model at each of these steps.                       

See also: [`solve`](@ref), [`SweepX`](@ref), [`setdof!`](@ref) 
"""
struct DirectXUA{NA,ND} <: AbstractSolver end 
function solve(::Type{DirectXUA{NA,ND}},pstate,verbose::𝕓,dbg;
    time::AbstractVector{𝕣},
    initialstate::State,
    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxΔa::ℝ=1e-5,
    saveiter::𝔹=false,
    maxLineIter::ℤ=50,β::𝕣=.5,γfac::𝕣=.5,γbot::𝕣=1e-8) where{NA,ND}

    model,dis             = initialstate.model, initialstate.dis
    out1,asm1             = prepare(AssemblyDirect    ,model,dis,(1,ND,ND,NA))
    out2,asm2             = prepare(AssemblyDirectLine,model)
    nstep                 = length(time)
    assemble!(out1,asm1,dis,model,initialstate,(dbg...,solver=:DirectXUA,phase=:sparsity))
    Lv,Lvv,bigasm         = preparebig(ND,NA,nstep,out1)

    state                 = [State{1,ND,ND,@NamedTuple{γ::Float64}}(copy(initialstate)) for timeᵢ ∈ time]    # TODO set state[step].time
    pstate[]              = state                                                                            # TODO pstate typestable???
    if saveiter
        stateiter         = Vector{typeof(state)}(undef,maxiter) 
        pstate[]          = stateiter
    end    
    assemble!(out2,asm2,dis,model,initialstate,(dbg...,solver=:DirectXUA,phase=:preliminary))
    out2.ming ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly primal-feasible"))
    out2.minλ ≤ 0 && muscadeerror(@sprintf("Initial point is not strictly dual-feasible"  ))
    γ = γ₀ = out2.Σλg/max(1,out2.npos)*γfac

    Δy²                   = Vector{𝕣 }(undef,nstep)

    Δt = time[2]-time[1]                                                                                      # TODO reformat input for cst step size

    local LU
    for iter              = 1                                                                                 # TODO 1:maxiter
        verbose && @printf("    iteration %3d, γ=%g\n",iter,γ)

        assemblebig!(Lvv,Lv,bigasm,asm1,model,dis,out1,state,nstep,Δt,NA,γ,(dbg...,solver=:DirectXUA,iter=iter))



#         try if iter==1 LU = lu(Lvv) 
#         else           lu!(LU ,Lvv)
#         end catch; muscadeerror(@sprintf("Lvv matrix factorization failed at iter=%i",iter));end
#         Δv               = LU\Lv 

#         Δa               = getblock(Δv,bigasm,nblock)
#         Δa²              = sum(Δa.^2)
#         for (step,state)   ∈ enumerate(state)
#             Δy           = getblock(Δv,bigasm,step  )
#             Δy²[step]    = sum(Δy.^2)
#             decrement!(state,0,Δy,Ydofgr)
#             decrement!(state,0,Δa,Adofgr)
#         end    
        
#         s  = 1.  
#         local  Σλg,npos 
#         for iline = 1:maxLineIter
#             ΣLa              .= 0   
#             minλ,ming         = ∞,∞
#             Σλg,npos          = 0.,0
#             for (step,state)  ∈ enumerate(state)
#                 assemble!(out2,asm2,dis,model,state,(dbg...,solver=:DirectXUAstepwise,phase=:linesearch,iter=iter,iline=iline,step=step))
#                 ΣLa         .+= out2.La 
#                 minλ          = min(minλ,out2.minλ)
#                 ming          = min(ming,out2.ming)
#                 Σλg          += out2.Σλg
#                 npos         += out2.npos
#             end
#             if minλ>0 && ming>0 
#                 verbose && @printf("    %3d line-iterations\n",iline)
#                 break#out of line search
#             end
#             iline==maxLineIter && muscadeerror(@sprintf("Line search failed at iter=%3d, iline=%3d, s=%7.1e",iter,iline,s))
#             Δs                = s*(β-1)
#             s                += Δs
#             for (step,state)  ∈ enumerate(state)
#                 decrement!(state,0,Δs*getblock(Δv,bigasm,step),Ydofgr)
#                 decrement!(state,0,Δs*Δa                      ,Adofgr)
#             end
#         end
#         γ                     = max(Σλg/max(1,npos)*γfac, γ₀*γbot)

#         if saveiter
#             stateiter[iter]     = copy.(state) 
#         end

#         if all(Δy².≤cΔy²)  && Δa²≤cΔa²  
#             verbose && @printf("\n    DirectXUA converged in %3d iterations.\n",iter)
#             verbose && @printf(  "    maxₜ(|ΔY|)=%7.1e  |ΔA|=%7.1e  \n",√(maximum(Δy²)),√(Δa²) )
#             verbose && @printf(  "    nel=%d, nvariables=%d, nstep=%d, niter=%d\n",getnele(model),nV,nstep,iter)
#             break#out of iter
#         end
#         iter<maxiter || muscadeerror(@sprintf("no convergence after %3d iterations. |ΔY|=%7.1e  |ΔA|=%7.1e \n",iter,√(maximum(Δy²)),√(Δa²)))
       end
#     return
end


