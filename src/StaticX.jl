
###--------------------- ASMstaticX: for good old static FEM

mutable struct OUTstaticX{Tλ,Tλx} 
    Lλ    :: Tλ
    Lλx   :: Tλx 
    α     :: 𝕣
end   
function prepare(::Type{OUTstaticX},model,dis) 
    dofgr              = allXdofs(model,dis)
    ndof               = getndof(dofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lλ                 = asmvec!(view(asm,1,:),dofgr,dis) 
    Lλx                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
    out                = OUTstaticX(Lλ,Lλx,∞)
    return out,asm,dofgr
end
function zero!(out::OUTstaticX)
    zero!(out.Lλ)
    zero!(out.Lλx)
    out.α = ∞    
end
function addin!(out::OUTstaticX,asm,iele,scale,eleobj::E,Λ,X::NTuple{Nxdir,<:SVector{Nx}},U,A, t,γ,dbg) where{E,Nxdir,Nx}
    if Nx==0; return end # don't waste time on Acost elements...  
    ΔX         = δ{1,Nx,𝕣}(scale.X)                 # NB: precedence==1, input must not be Adiff 
    Lλ,α       = getresidual(implemented(eleobj)...,eleobj,(∂0(X)+ΔX,),U,A, t,γ,dbg)
    Lλ         = Lλ .* scale.X
    add_value!(out.Lλ ,asm[1],iele,Lλ)
    add_∂!{1}( out.Lλx,asm[2],iele,Lλ)
    out.α      = min(out.α,α)
end


# function scaledresidual(scale,eleobj::AbstractElement, Xs::NTuple{Nxder},Us::NTuple{Nuder},As, t,γ,dbg) where{Nxder,Nuder} 
#     X     = NTuple{Nxder}(xs.*scale.X for xs∈Xs)  
#     U     = NTuple{Nuder}(us.*scale.U for us∈Us)
#     A     =       As.*scale.A
#     R,α   = getresidual(implemented(eleobj)...,eleobj, X,U,A, t,γ,dbg) 
#     hasnan(R) && muscadeerror(dbg,"NaN in a residual or its partial derivatives")
#     return R.*scale.Λ ,α
# end
###---------------------
struct StaticX end
getnder(::Type{StaticX}) = (nXder=1,nUder=1)
function solve(::Type{StaticX},pstate,verbose,dbg;time::AbstractVector{𝕣},
                    initialstate::State,
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxresidual::ℝ=∞,
                    saveiter::𝔹=false,γ0::𝕣=1.,γfac1::𝕣=.5,γfac2::𝕣=100.)
    # important: this code assumes that there is no χ in state.
    model,dis        = initialstate.model,initialstate.dis
    out,asm,dofgr    = prepare(OUTstaticX,model,dis)
    asmt,solt,citer  = 0.,0.,0
    cΔx²,cLλ²        = maxΔx^2,maxresidual^2
    state            = allocate(pstate,Vector{State{1,1}}(undef,saveiter ? maxiter : length(time))) # state is not a return argument so that data is not lost in case of exception
    s                = State{1,1}(initialstate) 
    local facLλx 
    for (step,t)     ∈ enumerate(time)
        s            = settime(s,t)
        γ            = γ0
        for iiter    = 1:maxiter
            citer   += 1
            asmt+=@elapsed assemble!(out,asm,dis,model,s, γ,(dbg...,solver=:StaticX,step=step,iiter=iiter))
            solt+=@elapsed try if step==1 && iiter==1
                facLλx = lu(out.Lλx) 
            else
                lu!(facLλx,out.Lλx) 
            end catch; muscadeerror(@sprintf("matrix factorization failed at step=%i, iiter=%i",step,iiter)) end
#            @show cond(Array(out.Lλx))
            solt+=@elapsed Δx  = facLλx\out.Lλ
            Δx²,Lλ²  = sum(Δx.^2),sum(out.Lλ.^2)
            solt+=@elapsed decrement!(s,0,Δx,dofgr)
            γ       *= γfac1*exp(-(out.α/γfac2)^2)
            verbose && saveiter && @printf("        iteration %3d, γ= %7.1e\n",iiter,γ)
            saveiter && (state[iiter]=State(s.Λ,deepcopy(s.X),s.U,s.A,s.time,γ,model,dis))
            if Δx²≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iiter √(Δx²) √(Lλ²)
                ~saveiter && (state[step]=State(s.Λ,deepcopy(s.X),s.U,s.A,s.time,γ,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf("no convergence in step %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",step,iiter,√(Δx²),maxΔx,√(Lλ²)^2,maxresidual))
        end
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(dofgr) length(time) citer citer/length(time)
    verbose && @printf "    Build  time = %s, (per iteration: %s, per iteration and element: %s)\n" showtime(asmt)  showtime(asmt/citer)  showtime(asmt/citer/getnele(model))
    verbose && @printf "    Solve  time = %s, (per iteration: %s, per iteration and dof:     %s)\n" showtime(solt)  showtime(solt/citer)  showtime(solt/citer/getndof(dofgr))
    return
end
