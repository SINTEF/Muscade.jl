### Assembler



"""
	EigX{ORDER}

A non-linear, time domain solver, that solves the problem time-step by time-step.
Only the `X`-dofs of the model are solved for, while `U`-dofs and `A`-dofs are unchanged.

- `SweepX{0}` is Newton-Raphson, with feasibility line-search, to handle inequality constraints. 
- `SweepX{1}` is implicit Euler, with feasibility line-search. 
- `SweepX{2}` is Newmark-β, with Newton-Raphson iterations and feasibility line search

IMPORTANT NOTE: Muscade does not allow elements to have state variables, for example, plastic strain,
or shear-free position for dry friction.  Where the element implements such physics, this 
is implemented by introducing the state as a degree of freedom of the element, and solving
for its evolution, *even in a static problem*, requires the use of `ORDER≥1`

An analysis is carried out by a call with the following syntax:

```
initialstate    = initialize!(model)
setdof!(initialstate,1.;class=:U,field=:λcsr)
states           = solve(SweepX{2};initialstate=initialstate,time=0:10)
```
# Named arguments to `solve`:
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `initialstate`      a `State`, obtain from `ìnitialize!` or `SweepX`.
- `time`              maximum number of Newton-Raphson iterations 
- `β=1/4`,`γ=1/2`     parameters to the Newmark-β algorithm. Dummy if `ORDER<2`
- `maxiter=50`        maximum number of equilibrium iterations at each step.
- `maxΔx=1e-5`        convergence criteria: norm of `X`. D
- `maxLλ=∞`           convergence criteria: norm of the residual. 
- `saveiter=false`    set to true so that output `states` contains the state
                      at the iteration of the last step analysed.  Useful to study
                      a step that fails to converge. 
- `maxLineIter=50`    Maximum number of iteration in the feasibility line search.
                      set to 0 to skip the line search (not recommended for models
                      with inequality constraints).
- `sfac=0.5`          Parameter in the line search for a feasible point. If a 
                      tentative result is not feasible, backtrack by a factor `sfac`.
                      If still not feasible, backtrack what is left by a factor `sfac`,
                      and so forth, up to `maxLineIter` times.
- `γfac=0.5`          Parameter for feasibility. For an inequality constraint `g(X)`
                      with reaction force `λ`, require `g(X)*λ==γ`, and multiply
                      `γ *= γfac` at each iteration.                            

# Output

A vector of length equal to that of the named input argument `time` containing the states at the time steps.                       

See also: [`solve`](@ref), [`initialize!`](@ref), [`findlastassigned`](@ref), [`studysingular`](@ref), [`DirectXUA`](@ref), [`FreqXU`](@ref)  
"""
struct        EigX{TYPE} <: AbstractSolver end
function solve(::Type{EigX{ℝ}},pstate,verbose,dbg;
                    initialstate::State,
                    nmode::𝕣=10) 
    model,dis        = initialstate.model,initialstate.dis
    out,asm,Xdofgr   = prepare(AssemblySweepX{2},model,dis)  
    ndof             = getndof(Xdofgr)

    state            = State{1,3,1}(copy(initialstate,SP=(γ=0.,))) 


    states           = allocate(pstate,Vector{typeof(state)}(undef,saveiter ? maxiter : length(time))) # states is not a return argument of this function.  Hence it is not lost in case of exception
    local facLλx 
        assemble!(out,asm,dis,model,state,(dbg...,solver=:SweepX,phase=:preliminary,step=step))
        for iiter    = 1:maxiter
            citer   += 1
            out.firstiter = firstiter = iiter==1
            out.line = false
            assemble!(out,asm,dis,model,state,(dbg...,solver=:SweepX,step=step,iiter=iiter))
            try if step==1 && firstiter  facLλx = lu(out.Lλx) 
            else                         lu!(facLλx, out.Lλx) 
            end catch; muscadeerror(@sprintf("matrix factorization failed at step=%i, iiter=%i",step,iiter)) end
            Δx       = facLλx\out.Lλ
            Δx²,Lλ²  = sum(Δx.^2),sum(out.Lλ.^2)
            if     ORDER==0  decr0!(state,Δx ,Xdofgr                      )
            elseif ORDER==1  decr1!(state,Δx ,Xdofgr,out.c                )
            elseif ORDER==2  decr2!(state,Δx ,Xdofgr,out.c,firstiter,x′,x″)
            end


            verbose && saveiter && @printf("        iteration %3d, γ= %7.1e\n",iiter,γ)
            saveiter && (states[iiter]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
            if Δx²≤cΔx² && Lλ²≤cLλ² 
                verbose && @printf "    step %3d converged in %3d iterations. |Δx|=%7.1e |Lλ|=%7.1e\n" step iiter √(Δx²) √(Lλ²)
                ~saveiter && (states[step]=State(state.time,state.Λ,deepcopy(state.X),state.U,state.A,state.SP,model,dis))
                break#out of the iiter loop
            end
            iiter==maxiter && muscadeerror(@sprintf("no convergence of step %3d after %3d iterations |Δx|=%g / %g, |Lλ|=%g / %g",step,iiter,√(Δx²),maxΔx,√(Lλ²)^2,maxLλ))
            state.SP     = (γ=state.SP.γ*γfac,)
        end
    end
    verbose && @printf "\n    nel=%d, ndof=%d, nstep=%d, niter=%d, niter/nstep=%5.2f\n" getnele(model) getndof(Xdofgr) length(time) citer citer/length(time)
    return
end
function decr2!(state,Δx ,Xdofgr,c,firstiter,x′,x″)
    a₁,a₂,a₃,b₁,b₂,b₃ = c.a₁,c.a₂,c.a₃,c.b₁,c.b₂,c.b₃
    if firstiter
        getdof!(state,1,x′,Xdofgr) 
        getdof!(state,2,x″,Xdofgr) 
        a        = a₂*x′+a₃*x″
        b        = b₂*x′+b₃*x″
        Δx′      = a₁*Δx + a
        Δx″      = b₁*Δx + b
    else
        Δx′      = a₁*Δx 
        Δx″      = b₁*Δx 
    end
    decrement!(state,1,Δx ,Xdofgr)
    decrement!(state,2,Δx′,Xdofgr)
    decrement!(state,3,Δx″,Xdofgr)
end
function decr1!(state,Δx ,Xdofgr,c)
    Δx′      = c.a₁*Δx 
    decrement!(state,1,Δx ,Xdofgr)
    decrement!(state,2,Δx′,Xdofgr)
end
function decr0!(state,Δx ,Xdofgr)
    decrement!(state,1,Δx ,Xdofgr)
end
