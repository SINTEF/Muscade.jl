###### DofGroups
abstract type DofGroup end
struct AllXdofs <: DofGroup 
    scale :: 𝕣1
end
function AllXdofs(model::Model,dis)
    scale  = Vector{𝕣}(undef,getndof(model,:X))
    for di ∈ dis
        for d ∈ di
            scale[d.index.X] = d.scale.X
        end
    end
    return AllXdofs(scale)
end
function Base.setindex!(s::State,x::𝕣1,gr::AllXdofs) 
    s.X[1] .= x.*gr.scale
end
Base.getindex(s::State,gr::AllXdofs) = s.X[1]./gr.scale
getndof(gr::AllXdofs) = length(gr.scale)

##### Solvers and their Addin
# NB: A solver may require several Assemblers.  Assemblers are object, solvers are functions.

# ASMstaticX: for good old static FEM
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
function addin!(asm::ASMstaticX,scale,ieletyp,iele,eleobj::E,Λ,X,U,A, t,ε,dbg)  where{E<:AbstractElement}
    Nx            = length(Λ)                   
    ΔX            = δ{1,Nx,𝕣}()                 # NB: precedence==1, input must not be Adiff 
    Lλ            = scaledresidual(scale,eleobj, (∂0(X)+ΔX,),U,A, t,ε,dbg)
    i             = Vector(asm.dis[ieletyp][iele].index.X)    
    asm.Lλ[ i  ] += value{1}(Lλ)            
    asm.Lλx[i,i] += ∂{1,Nx}(Lλ)                     
end
function StaticX(pstate,dbg;model::Model,time::AbstractVector{𝕣},
                    initial::State=State(model,Disassembler(model)),
                    maxiter::ℤ=50,maxΔy::ℝ=1e-5,maxresidual::ℝ=∞,
                    verbose::𝕓=true,saveiter::𝔹=false)
    # important: this code assumes that there is no χ in state.
    verb             = verbose
    verb && @printf "    StaticX solver\n\n"
    dis              = initial.dis
    asm              = ASMstaticX(model,dis)
    dofgr            = AllXdofs(model,dis)
    asmt,solt,citer  = 0.,0.,0
    cΔy²,cLλ²        = maxΔy^2,maxresidual^2
    state            = allocate(pstate,Vector{State}(undef,saveiter ? maxiter : length(time))) # state is not a return argument so that data is not lost in case of exception
    for (step,t)     ∈ enumerate(time)
        verb && @printf "    step %3d" step
        old          = step==1 ? initial : state[step-1]
        s            = State(old.Λ,deepcopy(old.X),old.U,old.A,t,0.,model,dis)
        y            = s[dofgr] # includes scaling
        for iiter    = 1:maxiter
            citer   += 1
            asmt+=@elapsed assemble!(asm,dis,model,s, 0.,(dbg...,solver=:StaticX,step=step,iiter=iiter))
            solt+=@elapsed Δy = try asm.Lλx\-asm.Lλ catch; muscadeerror(@sprintf("Incremental solution failed at step=%i, iiter=%i",step,iiter)) end
            Δy²,Lλ²  = sum(Δy.^2),sum(asm.Lλ.^2)
            y      .+= Δy
            s[dofgr] = y  # includes descaling

            saveiter && (state[iiter]=deepcopy(s))
            if Δy²≤cΔy² && Lλ²≤cLλ² 
                verb && @printf " converged in %3d iterations. |Δy|=%7.1e |Lλ|=%7.1e\n" iiter √(Δy²) √(Lλ²)
                ~saveiter && (state[step]=s)
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf(" no convergence after %3d iterations |Δy|:%g / %g, |R|:%g / %g",iiter,√(Δy²),maxΔy,√(Lλ²)^2,maxresidual))
        end
    end
    verb && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(dofgr) length(time) citer citer/length(time)
    verb && @printf "    Build  time = %s, (per iteration: %s, per iteration and element: %s)\n" showtime(asmt)  showtime(asmt/citer)  showtime(asmt/citer/getnele(model))
    verb && @printf "    Solve  time = %s, (per iteration: %s, per iteration and dof:     %s)\n" showtime(solt)  showtime(solt/citer)  showtime(solt/citer/getndof(dofgr))
end
