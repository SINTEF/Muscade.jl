"""
    Muscade.stackstring()

Produces a concise `Vector{String}` representation of the stacktrace, with
inner functions before outer functions
"""
function stackstring(cut)
    stk = stacktrace()
    str = Vector{String}(undef,0)
    for i∈cut+1:length(stk)
        if stk[i].file == Symbol(".\\boot.jl") || stk[i].func == Symbol("top-level scope")
            break
        end
        strᵢ = @sprintf("%s",stk[i].func)
        if !(length(strᵢ)≥3 && strᵢ[1:3]=="#_#") # @inline macros add a slice to the stack
            push!(str,strᵢ)
            if str[end]=="macro expansion"
                str[end]="macro"
            end
        end
    end
    str
end
"""
    Muscade.@dbg var1,...

Similar to @show, but prints code location.    

"""
macro dbg(ex...)
    fullname = String(__source__.file)
    _,file   = splitdir(fullname)
    line     = __source__.line
    stk      = gensym(:stk)
    stkᵢ     = gensym(:stkᵢ)
    code = quote
        #printstyled("● ";color=:light_green) # \mdlgblkcircle
        printstyled(@sprintf("● %s:%i\n",$file,$line);color=:light_green,bold=true)  \mdlgblkcircle
      #  printstyled(@sprintf("  %s\n",$fullname);color=:green)
        $stk = Muscade.stackstring(2)
 #       print("  ")
           printstyled(@sprintf("  %25s = ","call stack");color=:green)
         for $stkᵢ ∈ reverse($stk)  
#            printstyled("▸ ",$stkᵢ," ";color=:green)
            print("▸ ",$stkᵢ," ";color=:green)
        end
        println()
    end
    for exᵢ ∈ ex    
        code = quote 
            $(code)
            printstyled(@sprintf("  %25s = ",$(string(exᵢ)));color=:green)
            length($(string(exᵢ)))>80 && println()
            @printf("%s\n",$(exᵢ))
        end                     
    end
    code = quote 
        $(code)
        println()
    end                     
    return esc(code)
end

macro dbg1(ex)
    path,file = splitdir(String(__source__.file))
    line = __source__.line
    stk     = gensym(:stk)
    stklen  = gensym(:stklen)
    stkᵢ    = gensym(:stkᵢ)
    sho     = gensym(:sho)
    esc(quote
        printstyled("● ";color=:light_green) # \mdlgblkcircle
        printstyled(@sprintf("%-30s ",@sprintf("%s:%i",$file,$line));color=:light_black)
        $stk = Muscade.stackstring(1)
        $stklen = sum(length.($stk).+1)
        for $stkᵢ ∈ reverse($stk)  
            printstyled("▸";color=:yellow) # \smallblacktriangleleft
            printstyled($stkᵢ;color=:light_black)
        end
        $sho = @sprintf("%-25s = %s",$(Meta.quot(ex)),$(ex))
        if $stklen<50 && length($sho)<100  # oneliner
            print(" "^(50-$stklen))
#            println($sho)
            @printf("%-25s ",$(Meta.quot(ex)))
            printstyled(" = ";color=:yellow,bold=true)
            println($(ex))
        else       # multiline                             
            println()
            println("  ",$(Meta.quot(ex))," = ",$(ex))
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
