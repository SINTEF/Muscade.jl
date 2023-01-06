
###--------------------- ASMstaticX: for good old static FEM

struct OUTstaticX  
    Lλ    :: 𝕣1
    Lλx   :: SparseMatrixCSC{𝕣,𝕫} 
end   
function prepare(::Type{OUTstaticX},model,dis) 
    dofgr              = allXdofs(model,dis)
    ndof               = getndof(dofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{𝕫2}(undef,narray,neletyp)  
    Lλ                 = asmvec!(view(asm,1,:),dofgr,dis) 
    Lλx                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
    out                = OUTstaticX(Lλ,Lλx)
    return out,asm,dofgr
end
function zero!(out::OUTstaticX)
    out.Lλ        .= 0
    out.Lλx.nzval .= 0
end
function addin!(out::OUTstaticX,asm,iele,scale,eleobj,Λ,X,U,A, t,ε,dbg) 
    Nx                       = length(Λ)                   
    ΔX                       = δ{1,Nx,𝕣}()                 # NB: precedence==1, input must not be Adiff 
    Lλ                       = scaledresidual(scale,eleobj, (∂0(X)+ΔX,),U,A, t,ε,dbg)
    addin!(out.Lλ       ,asm[1],iele,value{1}(Lλ) )
    addin!(out.Lλx.nzval,asm[2],iele,∂{1,Nx}(Lλ)  )
end

###---------------------

function staticX(pstate,dbg;model::Model,time::AbstractVector{𝕣},
                    initial::State=State(model,Disassembler(model)),
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxresidual::ℝ=∞,
                    verbose::𝕓=true,saveiter::𝔹=false)
    # important: this code assumes that there is no χ in state.
    verb             = verbose
    verb && @printf "    staticX solver\n\n"
    dis              = initial.dis
    out,asm,dofgr    = prepare(OUTstaticX,model,dis)
    asmt,solt,citer  = 0.,0.,0
    cΔx²,cLλ²        = maxΔx^2,maxresidual^2
    state            = allocate(pstate,Vector{State}(undef,saveiter ? maxiter : length(time))) # state is not a return argument so that data is not lost in case of exception
    s                = initial 
    for (step,t)     ∈ enumerate(time)
        verb && @printf "    step %3d" step
        s            = settime(s,t)
        for iiter    = 1:maxiter
            citer   += 1
            asmt+=@elapsed assemble!(out,asm,dis,model,s, 0.,(dbg...,solver=:StaticX,step=step,iiter=iiter))
            solt+=@elapsed Δx = try out.Lλx\out.Lλ catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            Δx²,Lλ²  = sum(Δx.^2),sum(out.Lλ.^2)
            decrement!(s,0,Δx,dofgr)
            saveiter && (state[iiter]=State(s.Λ,deepcopy(s.X),s.U,s.A,s.time,0.,model,dis))
            if Δx²≤cΔx² && Lλ²≤cLλ² 
                verb && @printf " converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" iiter √(Δx²) √(Lλ²)
                ~saveiter && (state[step]=State(s.Λ,deepcopy(s.X),s.U,s.A,s.time,0.,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf(" no convergence after %3d iterations |Δx|:%g / %g, |Lλ|:%g / %g",iiter,√(Δx²),maxΔx,√(Lλ²)^2,maxresidual))
        end
    end
    verb && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(dofgr) length(time) citer citer/length(time)
    verb && @printf "    Build  time = %s, (per iteration: %s, per iteration and element: %s)\n" showtime(asmt)  showtime(asmt/citer)  showtime(asmt/citer/getnele(model))
    verb && @printf "    Solve  time = %s, (per iteration: %s, per iteration and dof:     %s)\n" showtime(solt)  showtime(solt/citer)  showtime(solt/citer/getndof(dofgr))
    return
end
