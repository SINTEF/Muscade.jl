function newEntryInBenchmarkHistory(modelID,runDuration)
    commit_hash         = readchomp(`git -C $(@__DIR__) rev-parse HEAD`)
    username            = ENV["USERNAME"]
    pc_id               = ENV["COMPUTERNAME"]
    benchmark_datetime  = Dates.now()
    julia_version       = VERSION
    open(joinpath(@__DIR__, "BenchmarkHistory.txt"), "a") do file
        println(file, "$modelID $commit_hash $julia_version $username $pc_id $benchmark_datetime $runDuration")
    end
end