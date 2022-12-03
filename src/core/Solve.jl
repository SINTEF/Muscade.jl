######## state and initstate
# at each step, contains the complete, unscaled state of the system
struct State{Nxder,Nuder,D}
    Λ     :: 𝕣1
    X     :: NTuple{Nxder,𝕣1}
    U     :: NTuple{Nuder,𝕣1}
    A     :: 𝕣1
    time  :: 𝕣
    ε     :: 𝕣
    model :: Model
    dis   :: D
end
# a constructor that provides an initial state
State(model::Model,dis;t=-∞) = State(zeros(getndof(model,:X)),(zeros(getndof(model,:X)),),(zeros(getndof(model,:U)),),zeros(getndof(model,:A)),t,0.,model,dis)
settime(s,t) = State(s.Λ,s.X,s.U,s.A,t,0.,s.model,s.dis)  


## find the last assigned array-element in a vector 
lastassigned(state) = state
function lastassigned(v::Vector)
    i = findlast([isassigned(v,i) for i=1:length(v)])
    return isnothing(i) ? nothing : lastassigned(v[i])
end

######### error management for solver
function solve(solver!::Function;verbose::𝕓=true,kwargs...) # e.g. solve(SOLstaticX,model,time=1:10)
    verbose && printstyled("\n\n\nMuscade\n",bold=true,color=:cyan)
    pstate = Ref{Any}()
    dbg    = ()
    try
        solver!(pstate,dbg;verbose=verbose,kwargs...) # 
    catch exn
        verbose && report(exn)
        verbose && printstyled("\nAborting the analysis.",color=:red)
        verbose && println(" Function `solve` should still be returning results obtained so far.")
    end
    verbose && printstyled("\nMuscade done.\n\n\n",bold=true,color=:cyan)
    return pstate[]
end

