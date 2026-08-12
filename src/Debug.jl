"""
    Muscade.stackstring()

Produces a concise `Vector{Sring}` representation of the stacktrace.
"""
function stackstring(cut)
    stk = stacktrace()
    str = Vector{String}(undef,0)
    for i∈cut+1:length(stk)
        if stk[i].file == Symbol(".\\boot.jl") || stk[i].func == Symbol("top-level scope")
            break
        end
        push!(str,@sprintf("%s",stk[i].func))
        if str[end]=="macro expansion"
            str[end]="macro"
        end
    end
    str
end
"""
    Muscade.@dbg var1,...

Similar to @show, but prints code location.    

"""
macro dbg(ex)
    path,file = splitdir(String(__source__.file))
    line = __source__.line
    stk     = gensym(:stk)
    stklen  = gensym(:stklen)
    stkᵢ    = gensym(:stkᵢ)
    sho     = gensym(:sho)
    esc(quote
        printstyled("@dbg ";color=:yellow,bold=true)
        printstyled(@sprintf("%-30s ",@sprintf("%s:%i",$file,$line));color=:light_black)
        $stk = Muscade.stackstring(1)
        $stklen = sum(length.($stk).+1)
        for $stkᵢ ∈ $stk  
#            printstyled(">";bold=true)
            print("►")
            printstyled($stkᵢ;color=:light_black)
        end
        $sho = @sprintf("%-15s = %s",$(Meta.quot(ex)),$(ex))
        if $stklen<50 && length($sho)<100  
            print(" "^(50-$stklen))
            println($sho)
        else                                    
            println()
            println($(Meta.quot(ex))," = ",$(ex))
        end
    end)
end

"""
    Muscade.scoop(A,[B...])

Gather variables into a `Scoop` exception, and throw the exception.  This provides a 
mechanism to "scoop up" arbitrary variables from inside a code - at the cost of interrupting it.

The scooped variables can be accessed in the REPL as
`A,B = Muscade.scoop(err)`

or by catching the `error`:

try
    fishy()
catch ex
    @show ex.scooped
end
"""
struct Scoop<: Exception scooped end
Base.showerror(io::IO, ex::Scoop) = print(io,"Scoop, data = Muscade.scoop(err) to retrieve scooped data") 
scoop(args...) = throw(Scoop(args))
scoop(err::Base.ExceptionStack) = err.stack[1].exception.error.scooped
