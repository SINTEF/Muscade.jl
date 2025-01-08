using Printf,MacroTools

######################## Helper functions
☼tag(s::Symbol)  = string(s)[1]=='☼' # \sun  
☼tag(s)          = false
♢tag(s::Symbol)  = string(s)[1]=='♢' # \diamond  
♢tag(s)          = false
☼tail(s::Symbol) = Symbol(string(s)[4:end])
tail(s::Symbol)  = Symbol(string(s)[2:end])
code_tuple(e...) = Expr(:tuple,e...)

const digits = Set("0123456789")
function newsym(name) 
    s = string(name)
    if s[end]∈digits && s[end-1]∈digits && s[end-2]∈digits && s[end-3]=='_'
        i = parse(Int64,s[end-2:end])
        s = s[1:end-4]
        return Symbol(s,"_",@sprintf("%03i",i+1))
    else
        return Symbol(s,"_001")
    end
end


######################## Request definition

"""

    req = @request expr

Create a request of internal results wanted from a function. Considering the function
presented as example for [`@espy`](@ref), examples of possible syntax include

```
req       = @request gp(s,z,material(a,b))
req       = @request gp(s)
req       = @request gp(material(a))
```

The first expression can be read as follows: "In the function, there is a `do` loop over variable `igp`, and within this
loop a call to a function `material`.  Results `s` and `z` are wanted from within the loop, and results `a` and `b`
from within `material`.

The corresponding datastructure containing the results for each element is a nesting of `NTuples` and `NamedTuples`, 
and can be accessed as `out.gp[igp].material.a` and so forth.        

See also: [`@espy`](@ref), [`@espydbg`](@ref)
"""
macro request(ex)
    if     @capture(ex,(args__,));      prettify(code_tuple(request_.(args)...))
    elseif @capture(ex,name_(args__));  prettify(code_tuple(request_(ex)))
    elseif @capture(ex,symbol_);        prettify(code_tuple(request_(ex)))
    else;                               muscadeerror("Not a valid request")
    end
end
request_(ex::Symbol) = quote $ex=nothing end
request_(ex        ) = @capture(ex,name_(args__)) ? quote $name=$(code_tuple(request_.(args)...)) end : muscadeerror("Not a valid request")

instantiate(::Type{Nothing}        )            = nothing
instantiate(::Type{NamedTuple{k,F}}) where{k,F} = NamedTuple{k}(instantiate.(fieldtypes(F)))
"""

    req = Muscade.merge(o.req)

Elements like [`ElementCost`](@ref) and [`ElementConstraint`](@ref) use requests to apply a cost or a constraint to
requestables from another "target" element. These cost/constraint elements must be coded carefully so that `getresult` can be used to extracted
both requestable internal results from the cost/constraint and from the target element.

`merge` (which is *not exported by Muscade*) is used to merge the requests for the request needed to enforce a cost or constraint, and the user's 
request for element to be obtained from the analysis.  The call to `merge`, to be inserted in the code of `lagrange` for the cost/constraint element
will be modified by `@espy`.  

See the code of `ElementCost`'s constructor and `lagrange` method for an example.

See also: [`ElementCost`](@ref),[`@request`](@ref),[`getresult`](@ref)
"""
@inline merge(x)            = x
merge(x        , ::Nothing) = x
merge(::Nothing, x        ) = x
merge(::Nothing, ::Nothing) = nothing
@generated function merge(NT1::NamedTuple, NT2::NamedTuple) 
    nt1     = instantiate(NT1)
    nt2     = instantiate(NT2)
    kout    = union(keys(nt1), keys(nt2))
    t       = ntuple(length(kout)) do ikey
        key = kout[ikey]
        v1  = get(nt1, key, nothing)
        v2  = get(nt2, key, nothing)
        v   = merge(v1, v2)
        :($key = $v) 
    end
    return code_tuple(t...)
    return c
end


######################## Generate new function code
# ✓check,✔Check,∎QED,⋆star,♢diamond,☼sun,☐Box
spaces = "|  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  "
function printtrace(trace::𝕫,s::String)
    if trace < 0;return end
    println(spaces[1:3trace],s)
end

## Clean code
function code_clean_function(ex::Expr  ) 
    if @capture(ex,@named(res__,))
        Expr(:tuple,(:($s=$s) for s∈res)...)
    elseif @capture(ex,left_=right_) && ♢tag(left)
        :(nothing) 
    else
        Expr(ex.head,[code_clean_function(a) for a ∈ ex.args]...)
    end
end
code_clean_function(ex::Symbol) = ☼tag(ex) ? ☼tail(ex) : ex
code_clean_function(ex        ) = ex

## Extractor code
# should such a function itself generate the new "out" variable?
function code_write_to_out(out,req,var,trace) 
    printtrace(trace,"code_write_to_out")
    newout = newsym(out)
    code   = quote
        $newout = haskey($req,$(QuoteNode(var))) ? ($out..., $var=$var) : $out 
    end 
    printtrace(trace,"done")
    return code,newout
end
function code_merge_req(left,oreq,req,out,trace=-999999) # req = merge(o.req)
    printtrace(trace,"code_merge_req")
    code   = quote
        $left  = merge($req,$oreq)
    end 
    printtrace(trace,"done")
    return code,out
end



function code_call(left,foo,args,req,out,trace=-999999) # left,... = foo(args)
    printtrace(trace,"code_call")
    out_new = newsym(out)
    out_foo = newsym(out)
    if left isa Symbol
        leftout = code_tuple(left,out_foo)
    elseif left isa Expr && left.head==:tuple
        leftout = code_tuple(left.args...,out_foo)
    else
        muscadeerror("Espy breakdown")
    end
    code   = quote
        if haskey($req,$(QuoteNode(foo)))
            $leftout            = $foo($(args...),$req.$foo,)  
            $out_new            = ($out...,$foo=$out_foo) 
        else
            $left               = $foo($(args...))
            $out_new            = $out
        end 
    end 
    printtrace(trace,"done")
    return code,out_new
end
function code_ntupledoloop(left,igp,ngp,body,req,out,trace=-999999) # foo(args) do var body end
    printtrace(trace,"code_ntupledoloop")
    gp              = tail(igp)   
    req_gp          = Symbol(req,"_",gp)   
    out_gp          = Symbol(out,"_",gp)
    out_new         = newsym(out)
    gpsym           = QuoteNode(gp)
    out_gp_new      = out_gp
    for i=1:length(body)-1
        body[i],out_gp_new = code_recursion(body[i],out_gp_new,req_gp,trace+1)
    end
    if     @capture(body[end],@named(res__,)) 
        tup = Expr(:tuple,(:($s=$s) for s∈res)...)  #@named(a,b) → (a=a,b=b)
        push!(tup.args,:(out=$out_gp_new))           #            → (a=a,b=b,out=out_gp_new)
        code            = quote
            $req_gp = haskey($req,$gpsym) ? $req.$gp : nothing
            $left   = ntuple($ngp) do $igp
                $out_gp = (;)
                $(body[1:end-1]...)
                $tup
            end
            $out_new = haskey($req,$gpsym) ? ($out..., $gp=NTuple{$ngp}($left[$igp].out for $igp=1:$ngp)) : $out
        end    
    else   muscadeerror("body of ntuple-do-loop must end with a NamedTuple containing the outputs of the iteration (;a=a) or (a=a,b=b)")
    end
    printtrace(trace,"done")
    return code,out_new
end

function code_assigment_rhs(left,right,out,req,trace=-999999) # work with the rhs of an assigment, return the whole assigment
    printtrace(trace,"code_assigment_rhs")
    trace += 1
    if @capture(right,  merge(oreq_)   )            
        printtrace(trace,"req = merge(o.req)")
        code,out_new = code_merge_req(left,oreq,req,out,trace+1)
        printtrace(trace,"done")
    elseif @capture(right,  foo_(args__)   ) &&  ☼tag(foo)           
        printtrace(trace,"left = ☼foo(args)")
        code,out_new = code_call(left,☼tail(foo),args,req,out,trace+1)
        printtrace(trace,"done")
    elseif @capture(right, mod_.foo_(args__)  )   &&  ☼tag(foo)  
        printtrace(trace,"left = mod.☼foo(args)")
        code,out_new = code_call(left,quote $mod.$(☼tail(foo)) end,args,req,out,trace+1)
        printtrace(trace,"done")
    elseif @capture(right, ntuple(ngp_) do igp_ body__ end)      
        printtrace(trace,"left = ntuple(ngp) do igp body end")
        code,out_new = code_ntupledoloop(left,igp,ngp,body,req,out,trace+1)
        printtrace(trace,"done")
    else
        printtrace(trace,"left = right")
        code = quote 
            $left = $right
        end
        out_new = out
        printtrace(trace,"done")
    end
    trace -= 1
    printtrace(trace,"done")
    return code,out_new
end
function code_recursion(ex::Expr,out,req,trace=-999999)
    printtrace(trace,"code_recursion")
    trace += 1
    if @capture(ex, for i_=it_ body_ end) || # syntaxes which contain an assigment must not recurse the assigment because
       @capture(ex, [a_ for i_=it_])
        out_new = out                     # 1) anotation is forbidden in this context
        code = ex                         # 2) processing an equality introduces a LineNumber, at an illegal place within the syntax
    elseif @capture(ex,  left_ = right_   )
        printtrace(trace,"assign")
        trace += 1
        if ☼tag(left)                       # ☼a = ...
            printtrace(trace,"☼a = ...")
            left              = ☼tail(left)   
            assigment,out_new = code_assigment_rhs(left,right,out,req,trace+1) 
            write,out_new     = code_write_to_out(out_new,req,left,trace+1)                                          
            code              = quote
                $assigment
                $write
                $left
            end
            printtrace(trace,"done")
        elseif @capture(left,   (args__,)    )                                               # (a,☼b) = ...
            printtrace(trace,"(a,☼b) = ...")
            write  = Vector{Expr}(undef,0)                                          # will contain the macros to insert
            @capture(left,   (args__,)    )
            left   = code_clean_function(left)
            assignement,out_new = code_assigment_rhs(left,right,out,req,trace+1)
            for arg ∈ args
                if ☼tag(arg) 
                    arg = ☼tail(arg)
                    w,out_new = code_write_to_out(out_new,req,arg,trace+1)
                    push!(write,w)                                        # @espy_record out req b
                end
            end
            if length(write)==0
                code = assignement
            else
                code = quote
                    $assignement
                    $(write...)
                    $left
                end
            end
            printtrace(trace,"done")
        elseif ♢tag(left)                       # ♢a = ...
            printtrace(trace,"♢a = ...")
            left              = ☼tail(left) 
            varsym            = QuoteNode(left)   
            assigment,out_tmp = code_assigment_rhs(left,right,out,req,trace+1) 
            out_new           = newsym(out_tmp)                                         
            code              = quote
                if haskey($req,$varsym)                                       # if haskey(req,:x)
                    $assigment
                    $out_new  = ($out...,$left = $left)
                else
                    $out_new  = $out
                end
            end
            printtrace(trace,"done")
        else
            printtrace(trace,"a = ...")
            code,out_new = code_assigment_rhs(left,right,out,req,trace+1)
            printtrace(trace,"done")
        end
        trace -= 1
        printtrace(trace,"done assign")
    elseif @capture(ex, return (args__,)) 
        code = quote
            return $(args...),$out    
        end
        out_new = out
    else
        printtrace(trace,"Expr recursion")
        args = ()
        for a ∈ ex.args
            codea,out = code_recursion(a,out,req,trace+1)
            args = (args...,codea)
        end
        code = Expr(ex.head,args...)
        out_new = out
        printtrace(trace,"done Expr recursion")
    end
    trace -= 1
    printtrace(trace,"done code_recursion")
    return code,out_new
end
function code_recursion(ex,out,req,trace=-999999) # to handle line-number nodes, Symbols... etc.
    printtrace(trace,"code_recursion (default)")
    printtrace(trace,"done")
    return ex,out
end
function code_espying_function(ex::Expr,out,req,trace=-999999) # ex must be a function declaration
    printtrace(trace,"code_espying_function")
    if @capture(ex, foo_(args__) = (body__,)) # oneliner functions can not have ☼ or ♢
        return quote
            $foo($(args...),req) = $(body...),nothing
        end
    elseif @capture(ex, foo_(args__) where T_= (body__,))
        return quote
            $foo($(args...),req) where $T = $(body...),nothing
        end
    else
        dict            = splitdef(ex)
        dict[:body]     = quote
            $out = (;)
            $(dict[:body]) 
        end
        push!(dict[:args],req)
        dict[:body],out = code_recursion(dict[:body],out,req,trace+1)
        printtrace(trace,"done")
        return combinedef(dict)
    end
end        

"""
    @espy function ... end

From an anotated function code, generate
    - "clean" code, in which the anotations have been deleted, and with 
      the call syntax `argout... = foo(argin...)`
    - "espying" code, with added input and ouput arguments
      `argout...,res = foo(argin...,req)` where `req` has been
      generated using `@request` and `res` is a nested structure
      of `NamedTuple`s and `NTuple`s containing the requested data.

The macro is not general: it is designed for `residual` and `lagrangian`,
which for performance have to be programmed in "immutable" style: they must
never mutate variables (this implies in particular, no adding into
an array in a loop over Gauss points). So `@espy` only supports the specific
programming constructs needed in this style.

The following is an example of anotated code:

```
@espy function residual(x::Vector{R},y) where{R<:Real}
    ngp=2
    accum = ntuple(ngp) do igp
        ☼z = x[igp]+y[igp]
        ☼s,☼t  = ☼material(z)
        ♢square = s^2
        @named(s) 
    end
    r = sum(i->accum[i].s,ngp)
    return r,nothing,nothing
end
```
- The keyword `function` is preceded by the macro-call `@espy`.  
- The name of requestable variables is preceded by `☼` (`\\sun`). 
  Such anotation must always appear on the left of an assigment. 
- If the name of a variable is preceded by `♢` (`\\diamond`), then
  the variable is evaluated only if requested. Such a notation can only
  be used if there is only one variable left of the assignement.
- The name of a function being called must be preceded by `☼` if the
  definition of the function is itself preceeded by the macro-call `@espy`.
- `for`-loops are not supported. `do`-loops must be used: to be efficient, 
  `residual` and `lagrangian` must not allocate and thus use immutables.
- The keyword `return` must be explicitly used, and if must be followed
  the a comma separated list of output variables. Syntaxes like
  `return if...` are not supported.  
                

See also: [`@request`](@ref), [`@espydbg`](@ref), [`getresult`](@ref)
"""
macro espy(ex)
    cntr = 0
    esc(quote
        $(code_clean_function(ex))
        $(code_espying_function(ex,newsym(:out),newsym(:req)))
    end)
end

"""
    @espydbg function ... end

Generate the same code as [`@espy`](@ref) and print it (for debug purposes).

See also: [`@espy`](@ref), [`@request`](@ref)
"""
macro espydbg(ex)
    cntr = 0
    bold = true
    printstyled("@espydbg: ";bold,color=:cyan)
    printstyled("Clean code\n";bold,color=:green)
    clean = code_clean_function(ex)
    println(prettify(clean))
    printstyled("\n@espydbg: ";bold,color=:cyan)
    printstyled("Espying code\n";bold,color=:yellow)
    espying = code_espying_function(ex,newsym(:out),newsym(:req))#,0) #trace
    println(prettify(espying)) 
    printstyled("\n@espydbg: ";bold,color=:cyan)
    printstyled("Done\n\n";bold,color=:cyan)
    esc(quote
            $clean
            $espying
        end)
end

