function newEntryInBenchmarkHistory(modelID,runDuration)
    commit_hash         = readchomp(`git -C $(@__DIR__) rev-parse HEAD`)
    username            = get(ENV, "USER", get(ENV, "USERNAME", "unknown"))
    pc_id               = gethostname()
    machine             = Sys.MACHINE
    benchmark_datetime  = Dates.now()
    julia_version       = VERSION
    open(joinpath(@__DIR__, "BenchmarkHistory.txt"), "a") do file
        println(file, "$modelID $commit_hash $julia_version $username $pc_id $machine $benchmark_datetime $runDuration")
    end
end