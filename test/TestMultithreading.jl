module TestMultithreading

using Test
using Muscade
using Muscade: run_all_jobs!, SingleDependencyScheduler

struct Job
    ijob :: 𝕫
    fail :: 𝕓
    time :: 𝕣
    done :: BitVector
    verbose :: 𝕓
end
function execute!(job::Job,error_ch,dbg;verbose=false,kwargs...) 
    job.done[job.ijob] && muscadeerror("Job required is already completed")
    tend = round(Int64,1000_000_000*job.time) + time_ns()
    while time_ns() < tend 
        if isready(error_ch)
            verbose && println("Job interrupted")
            return
        end
        if job.fail && time_ns()>tend/2
            muscadeerror("1202") # Now see where that got you!
        end
    end 
    job.done[job.ijob] = true
    return 
end
struct Accumulator
    sum :: Ref{𝕣}
end
function reduce!(accum,job,dbg)
    accum.sum[] += 1
end


## Test

# # Sequential scheduler

verbose = false
dbg = (;location=:testing)
@testset "SingleDependencyScheduler sequential" begin
    done     = falses(7) 
    jobs     = [Job(i,false,0.001+0.0001i,done,verbose) for i=1:7] 
    accum    = Accumulator(Ref(0.))
    sched    = SingleDependencyScheduler([4,3,0,0,1,2,5])  
    run_all_jobs!{:sequential}(accum,execute!,reduce!,jobs,sched,verbose,dbg) 
    @test all(done)
    @test accum.sum[] == 7
end

# Sequential scheduler with job error

@testset "SingleDependencyScheduler sequential with job error" begin
    done     = falses(7) 
    jobs     = [Job(i,i==7,0.001+0.0001i,done,verbose) for i=1:7] 
    accum    = Accumulator(Ref(0.))
    sched    = SingleDependencyScheduler([4,3,0,0,1,2,5]) 
    @test_throws Muscade.MuscadeException run_all_jobs!{:sequential}(accum,execute!,reduce!,jobs,sched,verbose,dbg)
    @test done == [1,1,1,1,1,1,0]
    @test accum.sum[] == 6
end

# SingleDependencyScheduler

@testset "SingleDependencyScheduler parallel" begin
    done     = falses(7) 
    sched    = SingleDependencyScheduler([4,3,0,0,1,2,5]) 
    jobs     = [Job(i,false,0.001+0.0001i,done,verbose) for i=1:7] 
    accum    = Accumulator(Ref(0.))
    run_all_jobs!{:parallel}(accum,execute!,reduce!,jobs,sched,verbose,dbg;nworker=3)
    @test all(done)
    @test accum.sum[] == 7
end

# SingleDependencyScheduler with job error

@testset "SingleDependencyScheduler parallel with job error" begin
    done     = falses(7) 
    sched    = SingleDependencyScheduler([4,3,0,0,1,2,5]) 
    jobs     = [Job(i,i==4,0.001+0.0001i,done,verbose) for i=1:7] 
    accum    = Accumulator(Ref(0.))
    @test_throws Muscade.ErrorInMultithread   run_all_jobs!{:parallel}(accum,execute!,reduce!,jobs,sched,verbose,dbg;nworker=3)
    @test all(done .≤ [0,1,1,1,0,1,0]) # what exactly will be done is a race
    @test accum.sum[] ≤ 4
end

end

