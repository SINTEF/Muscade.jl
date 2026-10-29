## Schedulers

struct SingleDependencyScheduler 
    dependson  :: Vector{Int64}
    dispatched :: Vector{Bool}
    done       :: Vector{Bool}
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
    waiting && (o.dispatched==o.done) && error("Circular dependencies") # we have not completed all jobs (waiting), but all dispatched jobs are done, so no new jobs will become doable
    return 0, waiting, done 
end
function done!(o::SingleDependencyScheduler,ijob)
    o.done[ijob]=true
end

struct SequentialScheduler
    jobs :: Vector{Int64}
    SequentialScheduler(order:: AbstractVector{Int64}) = isperm(order) ?  new(collect(order)) : error("Input must be a permutation")
end
# next! and done! not implemented, because it will dispatch to a version of weave that does not use these



## weave


function thread!(execute!,jobs,todo_ch,done_ch,error_ch,verbose)
    ithread = threadid()
    verbose && println("Thread $ithread: start thread")
    for ijob ∈ todo_ch 
        verbose && println("Thread $ithread: start job ",ijob)
        try
            execute!(jobs[ijob],error_ch,verbose)
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


function weave(execute!,jobs::AbstractVector,scheduler;nthread::Int=Threads.nthreads(),verbose=true) 
    local done
    todo_ch      = Channel{Int64}(length(jobs))
    done_ch      = Channel{Int64}(length(jobs))
    error_ch     = Channel{Any  }(1)

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
            sleep(0.01)
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
            execute!(jobs[ijob],error_ch,verbose)
        catch ex
            verbose && println("Interrupt sequential work due to error in job $ijob")
            close(error_ch)
            rethrow()
        end
    end
    close(error_ch)
    println("Done")
end




