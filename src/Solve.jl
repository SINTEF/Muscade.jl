using Printf
######### error management for solver
function solve(solver!;dbg=NamedTuple(),verbose::𝕓=true,silenterror::𝕓=false,kwargs...) # e.g. solve(SOLstaticX,model,time=1:10)
    verbose && printstyled("\n\n\nMuscade:",bold=true,color=:cyan)
    verbose && printstyled(@sprintf(" %s solver\n\n",nameof(solver!)),color=:cyan)

    pstate = Ref{Vector{State}}() # state is not a return argument of solver!, so that partial results are not lost on error
    try
        solver!(pstate,(dbg...,solver=nameof(solver!));verbose=verbose,kwargs...) # 
    catch exn
        silenterror || report(exn)
        silenterror || printstyled("\nAborting the analysis.",color=:red)
        silenterror || println(" Function `solve` should still be returning results obtained so far.")
    end
    verbose && printstyled("\nMuscade done.\n\n\n",bold=true,color=:cyan)
    return pstate[]
end

