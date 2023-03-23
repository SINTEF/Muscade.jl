using Printf,MacroTools

# TODO
# rewrite the @espy macro, test
# rewrite function codes (should not be necessary as they do not yet contain do loops)
# rewrite all extractor code (Ouput.jl, BasicElements.jl, all tests)    

######################## Helper functions
newsym(name)     = Symbol(name,"_",string(gensym())[3:end])
☼tag(s::Symbol)  = string(s)[1]=='☼' # \sun  
☼tag(s)          = false
♢tag(s::Symbol)  = string(s)[1]=='♢' # \diamond  
♢tag(s)          = false
☼tail(s::Symbol) = Symbol(string(s)[4:end])
tail(s::Symbol)  = Symbol(string(s)[2:end])
☼untag(s::Symbol)= ☼tag(s) ? ☼tail(s) : s
code_tuple(e...) = Expr(:tuple,e...)

######################## Request definition

"""

    req = @request expr

Create a request of internal results wanted from a function. Considering the function
presented as example for [`@espy`](@ref), examples of possible syntax include

    req       = @request gp(s,z,material(a,b))
    req       = @request gp(s)
    req       = @request gp(material(a))

The first expression can be read as follows: "In the function, there is a `do` loop over variable `igp`, and within this
loop a call to a function `material`.  Results `s` and `z` are wanted from within the loop, and results `a` and `b`
from within `material`.

The corresponding datastructure containing the results for each element is a nesting of `NTuples` and `NamedTuples`, 
and can be accessed as 

    out.gp[igp].material.a 

and so forth.        

See also: [`@espy`](@ref), [`@espydbg`](@ref)
"""
macro request(ex)
    if     @capture(ex,(args__,));      prettify(code_tuple(request_.(args)...))
    elseif @capture(ex,name_(args__));  prettify(code_tuple(request_(ex)))
    elseif @capture(ex,symbol_);        prettify(code_tuple(request_(ex)))
    else;                               error("Not a valid request")
    end
end
request_(ex::Symbol) = :($ex=nothing)
request_(ex        ) = @capture(ex,name_(args__)) ? :( $name=$(code_tuple(request_.(args)...))) : error("Not a valid request")


######################## Generate new function code
# ✓check,✔Check,∎QED,⋆star,♢diamond,☼sun,☐Box

## Clean code
code_clean_function(ex::Expr  ) = Expr(ex.head,[code_clean_function(a) for a ∈ ex.args]...)
code_clean_function(ex::Symbol) = ☼tag(ex) || ♢tag(ex) ? ☼tail(ex) : ex
code_clean_function(ex        ) = ex

## Extractor code
# should such a function itself generate the new "out" variable?
function code_write_to_out(out,req,var) 
    newout = newsym(out)
    code   = :($newout = haskey($req,$(QuoteNode(var))) ? ($out..., $var=$var) : $out) 
    return code,newout
end
function code_call(left,foo,args,out) # left = foo(args)
    out_new = newsym(out)
    out_foo = newsym(out)
    code   = quote
        if haskey($req,$(QuoteNode(foo)))
            $left,$out_foo = $foo($req.$foo,$(args...))  
            $out_new       = ($out...,$foo=$out_foo) 
        else
            $left         = $foo($(args...))
            $out_new       = $out
        end 
        $left   
    end 
    return code,out_new
end
function code_doloop(left,foo,igp,ngp,body,out) # foo(args) do var body end
    gp              = tail(igp)   
    req_gp          = Symbol(req,"_",gp)   
    out_gp          = Symbol(out,"_",gp)
    out_new         = newsym(out)
    gpsym           = QuoteNode(gp)
    @capture(body[end],(res__,)) || muscadeerror("body of do-loop must end with a Tuple containing the output of the iteration")
    out_gp_new      = out_gp
    for i=1:length(body)-1
        body[i],out_gp_new = code_recursion(body[i],out_gp,req_gp,trace)
    end
    code            = quote
        $req_gp = haskey($req,$gpsym) ? $req.$gp : nothing
        $out_gp = (;)
        $left   = $foo($(args...)) do $var
            $(body...)
            ($(res...),out=$out_gp_new)
        end
        $out_new = haskey($req,$gpsym) ? ($out..., $gp=NTuple{$ngp}($left[$igp].out for $igp=1:$ngp)) : $out
    end    
    return code,out_new
end

function code_assigment_rhs(left,right,out,req,trace=false) # work with the rhs of an assigment, return the whole assigment
    if @capture(right,  foo_(args__)   ) &&  ☼tag(foo)           
        trace && println("left = ☼foo(args)")
        code,out_new = code_call(left,☼tail(foo),args,out)
    elseif @capture(right, mod_.foo_(args__)  )   &&  ☼tag(foo)  
        trace && println("left = mod.☼foo(args)")
        code,out_new = code_call(left,:($mod.$(☼tail(foo))),args,out)
    elseif @capture(right, ntuple(len_) do igp_ body__ end)      
        trace && println("left = ntuple(ngp) do igp body end")
        code,out_new = code_doloop(left,foo,igp,ngp,body,out)
    else
        trace && println("left = right")
        code = :($left = $right)
        out_new = out
    end
    return code,out_new
end
function code_recursion(ex::Expr,out,req,trace=false)
    return if @capture(ex,  left_ = right_   )
        trace && println("assign")
        if ☼tag(left)                       # ☼a = ...
            trace && println("☼a = ...")
            left              = ☼tail(left)   
            assigment,out_new = code_assigment_rhs(left,right,out,req,trace) 
            write,out_new     = code_write_to_out(out_new,req,left)                                          
            code              = quote
                $assigment
                $write
                $left
            end
        elseif @capture(left,   (args__,)    )                                               # (a,☼b) = ...
            trace && println("(a,☼b) = ...")
            write  = Vector{Expr}(undef,0)                                          # will contain the macros to insert
            left   = ☼untag(left)
            assignement,out_new = code_assigment_rhs(left,right,out,req,trace)
            for arg ∈ args
                if ☼tag(arg) 
                    w,out_new = code_write_to_out(out_new,req,arg)
                    push!(write,w)                                        # @espy_record out req b
                end
            end
            code = quote
                $assignement
                $(write...)
                $left
            end
        elseif ♢tag(left)                       # ♢a = ...
            trace && println("♢a = ...")
            left               = ☼tail(left) 
            varsym            = QuoteNode(left)   
            assigment,out_tmp = code_assigment_rhs(left,right,out,req,trace) 
            out_new           = newsym(out_tmp)                                         
            code              = quote
                if haskey($req,$varsym)                                       # if haskey(req,:x)
                    $assigment
                    $out_new  = ($(out...),$left = $left)
                end
            end
        else
            trace && println("a = ...")
            code,out_new = code_assigment_rhs(left,right,out,req,trace)
        end
    else
        trace && println("recursion")
        code = Expr(ex.head,[code_recursion(a,out,req,trace) for a ∈ ex.args]...)
        out_new = out
    end
    return code,out_new
end
function code_recursion(ex,out,req,trace=false) # to handle line-number nodes, Symbols... etc.
    trace && println("default")
    return ex,out
end
function code_espying_function(ex::Expr,out,req,trace=false) # ex must be a function declaration
    trace && println("function")
    dict            = splitdef(ex)
    dict[:args]     = vcat([out,req],dict[:args])
    a,b = code_recursion(dict[:body],out,req,trace)
    dict[:body],out = code_recursion(dict[:body],out,req,trace)
    return combinedef(dict),out
end        

macro espy(ex)
    esc(quote
        $(code_clean_function(ex))
        $(code_espying_function(ex,newsym(:out),newsym(:req),false))
    end)
end

"""
    @espydbg function ...
    end

Generate the same code as [`@espy`](@ref) and print it (for debug purposes).

See also: [`@espy`](@ref), [`@request`](@ref)
"""
macro espydbg(ex)
    bold = true
    printstyled("@espydbg: ";bold,color=:cyan)
    printstyled("Clean code\n";bold,color=:green)
    println(prettify(code_clean_function(ex)))
    printstyled("\n@espydbg: ";bold,color=:cyan)
    printstyled("Extractor code\n";bold,color=:magenta)
    println(prettify(code_espying_function(ex,:out,:req,true)))
    println("\n")
    esc(quote
        $(code_clean_function(ex))
        $(code_espying_function(ex,newsym(:out),newsym(:req),false))
    end)
end

