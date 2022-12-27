###### DofGroup AllXdofs
abstract type DofGroup end
struct AllXdofs <: DofGroup 
    scale :: 𝕣1
end
function AllXdofs(model::Model,dis)
    scale  = Vector{𝕣}(undef,getndof(model,:X))
    for di ∈ dis
        for i ∈ di.index
            scale[i.X] = di.scale.X
        end
    end
    return AllXdofs(scale)
end
function decrement!(s::State,x::𝕣1,gr::AllXdofs) 
    s.X[1] .-= x.*gr.scale
end
Base.getindex(s::State,gr::AllXdofs) = s.X[1] # get values of a dofgroup (not used by solver)
getndof(gr::AllXdofs) = length(gr.scale)

##### ASMstaticX: for good old static FEM
struct ASMstaticX <: Assembler 
    dis   :: Vector{Any}          # naïve version! 
    Lλ    :: 𝕣1
    Lλx   :: SparseMatrixCSC{𝕣,𝕫} 
end #  
function ASMstaticX(model::Model,dis) 
    nX       = getndof(model,:X)
    return ASMstaticX(dis,zeros(nX),sparse(Int64[],Int64[],Float64[],nX,nX))
end
function zero!(asm::ASMstaticX)
    asm.Lλ  .= 0
    asm.Lλx .= 0
end
function addin!(asm::ASMstaticX,index,scale,eleobj,Λ,X,U,A, t,ε,dbg) 
    Nx            = length(Λ)                   
    ΔX            = δ{1,Nx,𝕣}()                 # NB: precedence==1, input must not be Adiff 
    Lλ            = scaledresidual(scale,eleobj, (∂0(X)+ΔX,),U,A, t,ε,dbg)
    i             = Vector(index.X)    
    asm.Lλ[ i  ] += value{1}(Lλ)            
    asm.Lλx[i,i] += ∂{1,Nx}(Lλ)                     
end
function staticX(pstate,dbg;model::Model,time::AbstractVector{𝕣},
                    initial::State=State(model,Disassembler(model)),
                    maxiter::ℤ=50,maxΔx::ℝ=1e-5,maxresidual::ℝ=∞,
                    verbose::𝕓=true,saveiter::𝔹=false)
    # important: this code assumes that there is no χ in state.
    verb             = verbose
    verb && @printf "    staticX solver\n\n"
    dis              = initial.dis
    asm              = ASMstaticX(model,dis)
    dofgr            = AllXdofs(model,dis)
    asmt,solt,citer  = 0.,0.,0
    cΔx²,cLλ²        = maxΔx^2,maxresidual^2
    state            = allocate(pstate,Vector{State}(undef,saveiter ? maxiter : length(time))) # state is not a return argument so that data is not lost in case of exception
    s                = initial 
    for (step,t)     ∈ enumerate(time)
        verb && @printf "    step %3d" step
        s            = settime(s,t)
        for iiter    = 1:maxiter
            citer   += 1
            asmt+=@elapsed assemble!(asm,dis,model,s, 0.,(dbg...,solver=:StaticX,step=step,iiter=iiter))
            solt+=@elapsed Δx = try asm.Lλx\asm.Lλ catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            Δx²,Lλ²  = sum(Δx.^2),sum(asm.Lλ.^2)
            decrement!(s,Δx,dofgr)
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
