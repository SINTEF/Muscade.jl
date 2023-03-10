using Printf
######### error management for solver
function solve(solver;dbg=NamedTuple(),verbose::𝕓=true,silenterror::𝕓=false,kwargs...) 
    verbose && printstyled("\n\n\nMuscade:",bold=true,color=:cyan)
    verbose && printstyled(@sprintf(" %s solver\n\n",Symbol(solver)),color=:cyan)

    nXdir,nUdir = getnder(solver)
    pstate = Base.RefValue{Vector{State{nXdir,nUdir}}}() # state is not a return argument of the solver, so that partial results are not lost on error
    try
        t = @elapsed solve(solver,pstate,verbose,(dbg...,solver=Symbol(solver));kwargs...)  
        verbose && @printf("    %s time: %s\n",Symbol(solver),showtime(t))
    catch exn
        silenterror || report(exn)
        silenterror || printstyled("\nAborting the analysis.",color=:red)
        silenterror || println(" Function 'solve' should still be returning results obtained so far.")
    end
    verbose && printstyled("Muscade done.\n\n\n",bold=true,color=:cyan)
    return pstate[]
end

