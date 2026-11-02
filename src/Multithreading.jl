## Schedulers

"""
    scheduler = SingleDependencyScheduler(dependencies)

`scheduler` is intended as an input to [`team`](@ref), and constrains how a collection
of jobs can be sequenced over multiple parallel threads.

`dependencies[i]::Int64` contains the index of the job that job `i` depends on.
Set `dependencies[i]` to `0` if job `i` does not require any other job to be completed.

See also: [team`](@ref), [`SequentialScheduler`](@ref)

"""
struct SingleDependencyScheduler 
    dependson  :: 𝕫1
    dispatched :: 𝕓1 #Vector{Bool }
    done       :: 𝕓1 #Vector{Bool }
end
SingleDependencyScheduler(dependson:: Vector{Int64}) = SingleDependencyScheduler(dependson,falses(length(dependson)),falses(length(dependson)))
function next!(o::SingleDependencyScheduler)
    for ijob ∈ eachindex(o.dependson)
        if ~o.dispatched[ijob] 
            if o.dependson[ijob]==0 || o.done[o.dependson[ijob]]  
                o.dispatched[ijob]=true
                return ijob,false,false # job, waiting, done
            end
        end
    end  
    done    = all(o.done)
    waiting = ~done
    waiting && (o.dispatched==o.done) && muscadeerror("Circular dependencies") # we have not completed all jobs (waiting), but all dispatched jobs are done, so no new jobs will become doable
    return 0, waiting, done 
end
function done!(o::SingleDependencyScheduler,ijob)
    o.done[ijob]=true
end

"""
    scheduler = SequentialScheduler(list)

`scheduler` is intended as an input to [`team`](@ref), and constrains how a collection
of jobs can be sequenced over multiple parallel threads.

`list[i]::Int64` contains the index of the job to be run at the `i`-th position.

    scheduler = SequentialScheduler(n)

where `n` is and `Int64`, is equivalent to `list=1:n`.    

See also: [team`](@ref), [`SingleDependencyScheduler`](@ref)

"""
struct SequentialScheduler
    jobs :: 𝕫1
    SequentialScheduler(order:: AbstractVector{𝕫}) = isperm(order) ?  new(collect(order)) : error("Input must be a permutation")
end
SequentialScheduler(n::𝕫) = SequentialScheduler(1:n)
# next! and done! not implemented, because SequentialScheduler will dispatch to a version of team that does not use these



## team


function worker!(execute!,iworker,jobs,todo_ch,done_ch,error_ch,verbose)
    verbose && println("Worker $iworker: clock in")
    for ijob ∈ todo_ch 
        verbose && println("Worker $iworker: start job ",ijob)
        try
            execute!(jobs[ijob],error_ch)
        catch ex
            put!(error_ch,(ex,catch_backtrace(),ijob,iworker))
            verbose && println("Worker $iworker: clock out due to error in job $ijob")
            return
        else
            verbose && println("Worker $iworker: close job ",ijob)
            put!(done_ch,ijob)
        end
    end
    verbose && println("Worker ",threadid(),": clock out")
    return
end

"""
    team(execute!,jobs,scheduler;nworker=nthreads(),verbose=true) 

Calls `execute!(jobs[i],error_ch)` for all `i`, where `execute!` and `jobs` are user-defined, 
accounting for dependencies between jobs defined by `scheduler`.

`scheduler` must implement `next(scheduler)` which returns the index of a job ready for execution
and `done!(scheduler,i)`, to tell the scheduler that job `i` is completed.

See also: [`SequentialScheduler`](@ref), [`SingleDependencyScheduler`](@ref). `SequentialScheduler` does
not use multithreading and is intended for debugging.

`execute` must execute `isready(error_ch) && return` regularly, so that all threads can be stopped if
one job fails.

"""
function team(execute!,jobs::AbstractVector,scheduler;nworker::Int=Threads.nthreads(),verbose=true) 
    local done
    todo_ch      = Channel{𝕫  }(length(jobs))
    done_ch      = Channel{𝕫  }(length(jobs))
    error_ch     = Channel{Any}(1           )

    for iworker ∈ 1:nworker
        @spawn worker!(execute!,iworker,jobs,todo_ch,done_ch,error_ch,verbose) 
    end

    while true # still work to do
        while true # still work doable
            ijob,waiting,done = next!(scheduler)
            (waiting || done) && break  
            verbose && println("scheduler: dispatching ijob $ijob")
            put!(todo_ch,ijob)  
        end    
        done && break 
        verbose && println("scheduler: go to sleep")   
        while ~isready(done_ch)
            sleep(.01)
            if isready(error_ch)
                verbose && println("scheduler: error in a thread, close all channels")
                close.((todo_ch,done_ch,error_ch))
                sleep(.01)
                throw(ErrorInMultithread(fetch(error_ch)...))
            end
        end
        ijob = take!(done_ch)
        done!(scheduler,ijob) 
        verbose && println("scheduler: wake up (job $ijob completed)")
    end   
    verbose && println("scheduler: close all channels") 
    close.((todo_ch,done_ch,error_ch))
    sleep(.01)
    verbose && println("All done!") 
    return 
end

function team(execute!::Function,jobs::AbstractVector,scheduler::SequentialScheduler;nworker::Int=Threads.nthreads(),verbose=true)
    error_ch     = Channel{Exception}(1)
    for i=1:length(scheduler.jobs)
        ijob = scheduler.jobs[i]
        verbose && println("Starting $ijob")
        try
            execute!(jobs[ijob],error_ch)
        catch ex
            verbose && println("Interrupt sequential work due to error in job $ijob")
            close(error_ch)
            rethrow()
        end
    end
    close(error_ch)
    println("Done")
end




