
######### error management for solver
function solve(solver!::Function;dbg=(),verbose::𝕓=true,kwargs...) # e.g. solve(SOLstaticX,model,time=1:10)
    verbose && printstyled("\n\n\nMuscade\n\n",bold=true,color=:cyan)
    pstate = Ref{Any}() # state is not a return argument of solver!, so that partial results are not lost on error
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

