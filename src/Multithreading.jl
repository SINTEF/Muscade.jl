## Schedulers

"""
    scheduler = SingleDependencyScheduler(dependencies)

`dependencies[i]::Int64` contains the index of the job that job `i` depends on.
Set `dependencies[i]` to `0` if job `i` does not require any other job to be completed.

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

`list[i]::Int64` contains the index of the job to be run at the `i`-th position.

    scheduler = SequentialScheduler(n)

where `n` is and `Int64`, is equivalent to `list=1:n`.    

"""
struct SequentialScheduler
    jobs :: 𝕫1
    SequentialScheduler(order:: AbstractVector{𝕫}) = isperm(order) ?  new(collect(order)) : error("Input must be a permutation")
end
SequentialScheduler(n::𝕫) = SequentialScheduler(1:n)
# next! and done! not implemented, because SequentialScheduler will dispatch to a version of weave that does not use these



## weave


function thread!(execute!,jobs,todo_ch,done_ch,error_ch,verbose)
    ithread = threadid()
    verbose && println("Thread $ithread: start thread")
    for ijob ∈ todo_ch 
        verbose && println("Thread $ithread: start job ",ijob)
        try
            execute!(jobs[ijob],error_ch)
        catch ex
            put!(error_ch,(ex,catch_backtrace(),ijob,ithread))
            verbose && println("Thread $ithread: close thread due to error in job $ijob")
            return
        else
            verbose && println("Thread $ithread: close job ",ijob)
            put!(done_ch,ijob)
        end
    end
    verbose && println("Thread ",threadid(),": close thread")
    return
end

"""
    weave(execute!,jobs,scheduler;nthread=nthreads(),verbose=true) 

Calls `execute!(jobs[i],error_ch)` for all `i`, where `execute!` and `jobs` are user-defined, 
accounting for dependencies between jobs defined by `scheduler`.

`scheduler` must implement `next(scheduler)` which returns the index of a job ready for execution
and `done!(scheduler,i)`, to tell the scheduler that job `i` is completed.

See also: [`SequentialScheduler`](@ref), [`SingleDependencyScheduler`](@ref). `SequentialScheduler` does
not use multithreading and is intended for debugging.

`execute` must execute `isready(error_ch) && return` regularly, so that all threads can be stopped if
one job fails.

"""
function weave(execute!,jobs::AbstractVector,scheduler;nthread::Int=Threads.nthreads(),verbose=true) 
    local done
    todo_ch      = Channel{𝕫  }(length(jobs))
    done_ch      = Channel{𝕫  }(length(jobs))
    error_ch     = Channel{Any}(1           )

    for ithread ∈ 1:nthread
        @spawn thread!(execute!,jobs,todo_ch,done_ch,error_ch,verbose) 
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

function weave(execute!::Function,jobs::AbstractVector,scheduler::SequentialScheduler;nthread::Int=Threads.nthreads(),verbose=true)
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




