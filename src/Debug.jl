function showstack()
    s = stacktrace()
    for sᵢ∈s
        @printf(">>>%-40s %s\n", sᵢ.func, sᵢ.file) 
    end
end

"""
    Muscade.stackstring()

Produces a concise string representation of the stacktrace.
"""
function stackstring(cut)
    stk = stacktrace()
    str = ""
    for i∈cut+1:length(stk)
        if stk[i].file == Symbol(".\\boot.jl") || stk[i].func == Symbol("top-level scope")
            break
        end
        str = @sprintf(">%s%s",stk[i].func,str)
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
    stk = gensym()
    sho = gensym()
    esc(quote
        printstyled("@dbg ";color=:yellow)
        printstyled(@sprintf("%-30s ",@sprintf("%s:%i",$file,$line));color=:light_black)
        $stk = Muscade.stackstring(1)
        $sho = @sprintf("%-15s = %s",$(Meta.quot(ex)),$(ex))
        if length($stk)<50 && length($sho)<100    
            printstyled(@sprintf("%-50s "   ,$stk);color=:blue)
            println($sho)
        else                                    
            printstyled(@sprintf("%s\n     ",$stk);color=:blue)
            println($(Meta.quot(ex))," = ",$(ex))
        end
    end)
end