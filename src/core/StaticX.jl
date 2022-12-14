
###--------------------- ASMstaticX: for good old static FEM

struct OUTstaticX  
    L位    :: 饾暎1
    L位x   :: SparseMatrixCSC{饾暎,饾暙} 
end   
function prepare(::Type{OUTstaticX},model,dis) 
    dofgr              = allXdofs(model,dis)
    ndof               = getndof(dofgr)
    narray,neletyp     = 2,getneletyp(model)
    asm                = Matrix{饾暙2}(undef,narray,neletyp)  
    L位                 = asmvec!(view(asm,1,:),dofgr,dis) 
    L位x                = asmmat!(view(asm,2,:),view(asm,1,:),view(asm,1,:),ndof,ndof) 
    out                = OUTstaticX(L位,L位x)
    return out,asm,dofgr
end
function zero!(out::OUTstaticX)
    out.L位        .= 0
    out.L位x.nzval .= 0
end
function addin!(out::OUTstaticX,asm,iele,scale,eleobj,螞,X,U,A, t,蔚,dbg) 
    Nx                       = length(螞)                   
    螖X                       = 未{1,Nx,饾暎}()                 # NB: precedence==1, input must not be Adiff 
    L位                       = scaledresidual(scale,eleobj, (鈭?0(X)+螖X,),U,A, t,蔚,dbg)
    addin!(out.L位       ,asm[1],iele,value{1}(L位) )
    addin!(out.L位x.nzval,asm[2],iele,鈭倇1,Nx}(L位)  )
end

###---------------------

function staticX(pstate,dbg;model::Model,time::AbstractVector{饾暎},
                    initial::State=State(model,Disassembler(model)),
                    maxiter::鈩?=50,max螖x::鈩?=1e-5,maxresidual::鈩?=鈭?,
                    verbose::饾晸=true,saveiter::饾敼=false)
    # important: this code assumes that there is no 蠂 in state.
    verb             = verbose
    verb && @printf "    staticX solver\n\n"
    dis              = initial.dis
    out,asm,dofgr    = prepare(OUTstaticX,model,dis)
    asmt,solt,citer  = 0.,0.,0
    c螖x虏,cL位虏        = max螖x^2,maxresidual^2
    state            = allocate(pstate,Vector{State}(undef,saveiter ? maxiter : length(time))) # state is not a return argument so that data is not lost in case of exception
    s                = initial 
    for (step,t)     鈭? enumerate(time)
        verb && @printf "    step %3d" step
        s            = settime(s,t)
        for iiter    = 1:maxiter
            citer   += 1
            asmt+=@elapsed assemble!(out,asm,dis,model,s, 0.,(dbg...,solver=:StaticX,step=step,iiter=iiter))
            solt+=@elapsed 螖x = try out.L位x\out.L位 catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            螖x虏,L位虏  = sum(螖x.^2),sum(out.L位.^2)
            decrement!(s,0,螖x,dofgr)
            saveiter && (state[iiter]=State(s.螞,deepcopy(s.X),s.U,s.A,s.time,0.,model,dis))
            if 螖x虏鈮螖x虏 && L位虏鈮L位虏 
                verb && @printf " converged in %3d iterations. |螖x|=%7.1e |L位|=%7.1e\n" iiter 鈭?(螖x虏) 鈭?(L位虏)
                ~saveiter && (state[step]=State(s.螞,deepcopy(s.X),s.U,s.A,s.time,0.,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf(" no convergence after %3d iterations |螖x|:%g / %g, |L位|:%g / %g",iiter,鈭?(螖x虏),max螖x,鈭?(L位虏)^2,maxresidual))
        end
    end
    verb && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(dofgr) length(time) citer citer/length(time)
    verb && @printf "    Build  time = %s, (per iteration: %s, per iteration and element: %s)\n" showtime(asmt)  showtime(asmt/citer)  showtime(asmt/citer/getnele(model))
    verb && @printf "    Solve  time = %s, (per iteration: %s, per iteration and dof:     %s)\n" showtime(solt)  showtime(solt/citer)  showtime(solt/citer/getndof(dofgr))
    return
end
