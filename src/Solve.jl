using Printf
abstract type AbstractSolver  end
"""
    solve(Solver;dbg=(;),verbose=true,silenterror=false,kwargs...)

Execute an analysis using `Solver`, and safeguard partial results in the
case of error. 

# Named arguments
- `dbg=(;)`           a named tuple to trace the call tree (for debugging)
- `verbose=true`      set to false to suppress printed output (for testing)
- `silenterror=false` set to true to suppress print out of error (for testing) 
- `kwargs...`         Further arguments passed on to the method `solve` provided by the solver

This will call the method `solve` provided by the solver with
`solve(Solver,pstate,verbose,(dbg...,solver=Symbol(Solver));kwargs...)`

See also: [`SweepX`](@ref), [`DirectXUA`](@ref), [`initialize!`](@ref) 
"""
function solve(Solver::Type{<:AbstractSolver};dbg=NamedTuple(),verbose::𝕓=true,silenterror::𝕓=false,catcherror::𝕓=true,kwargs...) 
    verbose && printstyled("\n\n\nMuscade:",bold=true,color=:cyan)
    verbose && printstyled(@sprintf(" %s solver\n\n",Symbol(Solver)),color=:cyan)
    pstate = Ref{Any}() # state is not a return argument of the solver, so that partial results are not lost on error
    if catcherror
        try
            t = @elapsed solve(Solver,pstate,verbose,(dbg...,solver=Symbol(Solver));kwargs...)  
            verbose && @printf("    %s time: %s\n",Symbol(Solver),showtime(t))
        catch exn
            silenterror || report(exn)
            silenterror || printstyled("\nAborting the analysis.",color=:red)
            silenterror || println(" Function 'solve' still returns any results obtained before the exception.")
        end
    else
        t = @elapsed solve(Solver,pstate,verbose,(dbg...,solver=Symbol(Solver));kwargs...)  
        verbose && @printf("    %s time: %s\n",Symbol(Solver),showtime(t))
    end
    verbose && printstyled("Muscade done.\n\n\n",bold=true,color=:cyan)
    return pstate[]
end

