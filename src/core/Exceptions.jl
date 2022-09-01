## NOT OPERATIONAL

struct MuscadeException <:Exception
    msg::String
    dbg::NamedTuple
end
"""
    muscadeerror([[dbg,]msg])

Throw a `MuscadeException`, where
- `dbg` is a `NamedTuple` that contains "location information"
(for example: step, increment, iteration, element, quadrature point) that will Basedisplayed with the error message.
- `msg` is a `String` describing the problem.
"""
muscadeerror(dbg::NamedTuple,msg)      = throw(MuscadeException(msg,dbg))
muscadeerror(msg)                      = throw(MuscadeException(msg,(;)))
muscadeerror()                         = throw(MuscadeException("" ,(;)))
relativebacktrace()                    = setdiff(catch_backtrace(),backtrace())[2:end-1]
function showerr(exn::MuscadeException)
    printstyled("Muscade error: ",color=:red)
    print(exn.msg)
    exn.dbg==(;) || print("\n",string(exn.dbg))
end
function showerr(exn::Exception)
    printstyled("Error: ",color=:red)
    showerror(stdout, exn)
end
report(::Exception)    = rethrow()
function report(::MuscadeException)
    println("")
    cs = Base.catch_stack()
    nex = length(cs)
    for iex = 1:nex-1
        showerr(cs[iex][1])
        Base.show_backtrace(stdout,setdiff(cs[iex][2],cs[iex+1][2])[2:end-1])
        print("\n\nthen caused ")
    end
    showerr(cs[nex][1])
    Base.show_backtrace(stdout,setdiff(cs[nex][2],backtrace())[2:end-1])
    print("\n\n")
end
