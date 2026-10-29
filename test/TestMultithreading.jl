module TestMultithreading

using Test
using Muscade
using Muscade: weave, SequentialScheduler, SingleDependencyScheduler

struct Job
    ijob :: Int64
    fail :: Bool
    time :: Float64
    done :: BitVector
end
function execute!(job::Job,error_ch,verbose) 
    tend = round(Int64,1000_000_000*job.time) + time_ns()
    while time_ns() < tend 
        if isready(error_ch)
            verbose && println("Job interrupted")
            break
        end
        if job.fail && time_ns()>tend/2
            muscadeerror("1202") # Now see where that got you!
        end
    end 
    job.done[job.ijob] && muscadeerror("run twice?")
    job.done[job.ijob] = true
    return 
end


## Test

# Sequential scheduler


@testset "SequentialScheduler" begin
    done     = falses(7) 
    jobs     = [Job(i,false,0.01+0.001i,done) for i=1:7] 
    sched    = SequentialScheduler([4,1,2,3,6,7,5]) 
    weave(execute!,jobs,sched,verbose=false) 
    @test all(done)
end

# Sequential scheduler with job error

@testset "SequentialScheduler with job error" begin
    done     = falses(7) 
    jobs     = [Job(i,i==7,0.01+0.001i,done) for i=1:7] 
    sched    = SequentialScheduler([4,1,2,3,6,7,5]) 
    @test_throws Muscade.MuscadeException weave(execute!,jobs,sched,verbose=false)
    @test done == [1,1,1,1,0,1,0]
end

# SingleDependencyScheduler

@testset "SingleDependencyScheduler" begin
    done     = falses(7) 
    sched    = SingleDependencyScheduler([4,3,0,0,1,2,5]) 
    jobs     = [Job(i,false,0.01+0.001i,done) for i=1:7] 
    weave(execute!,jobs,sched,nthread=3,verbose=false)
    @test all(done)
end

# SingleDependencyScheduler with job error

@testset "SingleDependencyScheduler with job error" begin
    done     = falses(7) 
    sched    = SingleDependencyScheduler([4,3,0,0,1,2,5]) 
    jobs     = [Job(i,i==4,0.01+0.001i,done) for i=1:7] 
    @test_throws Muscade.ErrorInMultithread   weave(execute!,jobs,sched,nthread=3,verbose=false)
    @test all(done .≤ [0,1,1,1,0,1,0]) # what exactly will be done is a race
end

end

