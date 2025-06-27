
struct MuscadeException <:Exception
    msg::String
    dbg::NamedTuple
end
"""
    muscadeerror([[dbg,]msg])

Throw a `MuscadeException`, where
- `dbg` is a `NamedTuple` that contains "location information"
(for example: solver, step, iteration, element, quadrature point) that will be displayed with the error message.
- `msg` is a `String` describing the problem.
"""
muscadeerror(dbg::NamedTuple,msg)      = throw(MuscadeException(msg,dbg))
muscadeerror(msg)                      = throw(MuscadeException(msg,(;)))
muscadeerror()                         = throw(MuscadeException("" ,(;)))
# relativebacktrace()                    = setdiff(catch_backtrace(),backtrace())[2:end-1]
function Base.showerror(io::IO, e::MuscadeException)
    printstyled(io,"MuscadeException: ",color=:red,bold=true)
    println(io, e.msg)
    if e.dbg≠(;)
        println(io,e.dbg)
    end
end
report(::Exception)    = rethrow()
function report(::MuscadeException)
        print("\n\n")
    cs = Base.catch_stack()
    nex = length(cs)
    for iex = 1:nex-1
        showerror(stderr,cs[iex][1])
        Base.show_backtrace(stdout,setdiff(cs[iex][2],cs[iex+1][2])[2:end-1])
        print("\n\nthen caused ")
    end
    showerror(stderr,cs[nex][1])
    Base.show_backtrace(stdout,setdiff(cs[nex][2],backtrace())[2:end-1])
    print("\n\n")
end
function muscadewarning(str,indent=0)
    for i = 1:indent
        print(" ")
    end
    printstyled("Warning: ",bold=true,color=:red)
    printstyled(str,color=:red)
    print("\n")
end