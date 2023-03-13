# manipulation of data-structures representing Julia expressions (::Expr)

# REPRISE

using Printf,MacroTools

## Generate unique symbol
newsym(name) = Symbol(name,"_",string(gensym())[3:end])

## Tests
issymbol(e::Symbol)         = true
issymbol(e)                 = false
isquote(e::QuoteNode)       = true
isquote(e)                  = false
isnumber(e::Number)         = true
isnumber(e)                 = false
isline(e::LineNumberNode)   = true
isline(e)                   = false
isexpr(e::Expr)             = true
isexpr(e,f...)              = false
isexpr(s,e::Expr)           = e.head == s

istuple(e)                  = isexpr(e) && e.head == :tuple  # (a,b)
isref(e)                    = isexpr(e) && e.head == :ref    # a[b]
isdot(e)                    = isexpr(e) && e.head == :(.)    # a.b
issub(e)                    = isexpr(e) && e.head == :(.)  && isquote(getright(e))  # a.b
isdeclare(e)                = isexpr(e) && e.head == :(::)   # a::b
isassign(e)                 = isexpr(e) && e.head == :(=)
iscall(e)                   = isexpr(e) && e.head == :call
isfunc(e)                   = isexpr(e) && e.head == :function
isfor(e)                    = isexpr(e) && e.head == :for
isescsymbol(e)              = isexpr(e) && e.head == (:escape) && issymbol(getleft(e))
issymbolish(e)              = issymbol(e)||isescsymbol(e)

## Analysis
getright(e)                 = e.args[2]
getleft(e)                  = e.args[1]
getfor(e)                   = (var=e.args[1].args[1],from=e.args[1].args[2].args[2],to=e.args[1].args[2].args[3],body=e.args[2]) #(var,from,to,body)

## Synthesis
maketuple(e...)             = Expr(:tuple,          e...)
makeref(e...)               = Expr(:ref,            e...)
#makecall(e...)              = Expr(:call,           e...)
makeblock(e...)             = Expr(:block,          e...)
#makewhere(e...)             = Expr(:where,          e...)
#makeassign(e...)            = Expr(:(=),            e...)
#makefunction(f,e...)        = Expr(:function,f,makeblock(e...))
#makefor(e...)               = Expr(:for,e...)
#makemacro(m,e...)           = Expr(:macrocall,      Symbol("@",m),:LineNumberHere,e...)
#makestructure(e...)         = Expr(:struct,false,   e...)
#makemutable(e...)           = Expr(:struct,true ,   e...)
makedeclare(e...)           = Expr(:(::),e[1],e[2])
makedot(e...)               = Expr(:(.) ,e[1],e[2])
makesub(e...)               = Expr(:(.) ,e[1],QuoteNode(e[2]))




# ## For human readable code, use println(ex) or println(deline(ex))
# # To understand the tree structure of an Expression, use pretty(ex)
# function dent(i)
#     for _ = 1:i
#         @printf("   ")
#     end
# end
# pretty(e::Union{Expr,Symbol,LineNumberNode,Symbol,Nothing})   = pretty(e,0)
# function pretty(e::Expr,ind)
#     dent(ind)
#     @printf("%s\n",e.head)
#     for iarg = 1:length( e.args)
#         pretty(e.args[iarg],ind+1)
#     end
# end
# function pretty(e::Symbol,ind)
#     dent(ind)
#     @printf("%s\n",e)
# end
# function pretty(::Nothing,ind)
#     dent(ind)
#     @printf("Nothing()\n")
# end
# function pretty(e::Number,ind)
#     dent(ind)
#     @printf("%g\n",e)
# end
# function pretty(e::LineNumberNode,ind)
#     dent(ind)
#     @printf("LineNumberNode\n")
# end
# function pretty(x,ind)
#     dent(ind)
#     display(x)
# end
# function pretty(m::String)
#     @printf("%s\n",m)
# end
# macro pretty(ex)
#     pretty(ex)
#     return ex
# end

####################### Definition by element of what can be requested
"""

    forloop

Component to build the `requestable` input to [`makekey`](@ref)
See also: [`makekey`](@ref), [`scalar`](@ref)
"""
struct forloop
    range :: Int
    body
end
"""

    scalar

Component to build the `requestable` input to [`makekey`](@ref)
See also: [`makekey`](@ref), [`forloop`](@ref)
"""
const scalar = ()#Int64[]
######################## Request definition

# Julia parses a.(b::Tb,c.d).(e::Te,f::Tf) as (a.(b::Tb,c.d)).(e::Te,f::Tf)
# transform to a.((b::Tb,c.d).(e::Te,f::Tf))
function leftparse(e)
    if     issub(e)     && issub(getleft(e))    leftparse(makedot(getleft(getleft(e)),makedot(    getright(getleft(e)),getright(e))))
    elseif isref(e)     && issub(getleft(e))    leftparse(makedot(getleft(getleft(e)),makeref(    getright(getleft(e))            )))
    elseif isdeclare(e) && issub(getleft(e))    leftparse(makedot(getleft(getleft(e)),makedeclare(getright(getleft(e)),getright(e))))
    elseif isexpr(e)                            Expr(e.head,leftparse.(e.args)...)
    elseif isquote(e)                           e.value
    else                                        e
    end
end

# Given a request, recursively generate the "key" to access big output arrays
# The input is an Expr (or Symbol), the output a tree structure of arrays (over GP) and
# NamedTuple, with leaves that are Int64 or AbstractArray of Int64 - containing indices into the output array
namedtuple(sym,val)         = (;zip(sym,val)...)
namedtuple(sym::Symbol,val) = (;zip((sym,),(val,))...)
isloopreq(ele) = isdot(ele) && isref(getleft(ele))
iscallreq(ele) = isdot(ele) && ~isref(getleft(ele))
function makekey_tuple(cnt,ex,reqabl) # ex=(a,gp[].(...),foo.(...))   reqabl=(a=[...],gp=forloop(ngp,(...)),foo=(...))
    if ex==() || ex==(;) return NamedTuple(),0 end
    if  isquote(ex) ex = maketuple(ex.value) end # allow user to type a.b for (a.(b,),)
    if ~istuple(ex) ex = maketuple(ex)       end
    len   = length(ex.args)
    elkey = Vector{Any   }(undef,len) # "Any" is OK here, will be used to build a tuple
    sym   = Vector{Symbol}(undef,len)
    for (i,ele) = enumerate(ex.args)
        if issymbol(ele)    # :ele
            sym[i]       = ele
            elkey[i],cnt = makekey_symbol(cnt,reqabl[ele])
        elseif isloopreq(ele) # sym[].(body)
            sym[i]       = getleft(getleft(ele))
            body         = getright(ele)
            elkey[i],cnt = makekey_loop(cnt,body,reqabl[sym[i]])
        elseif iscallreq(ele) # sym.(body)
            sym[i]       = getleft(ele)
            body         = getright(ele)
            elkey[i],cnt = makekey_tuple(cnt,body,reqabl[sym[i]])
        else
            error("Illegal request expression")
        end
    end
    return namedtuple(sym,elkey),cnt
end
function makekey_symbol(cnt,siz) # ex=:a  reqabl=[...]
    len     = prod(siz)
    ndim    = length(siz)
    key     = ndim>0 ? collect(reshape(cnt+1:cnt+len,siz)) : cnt+1
    cnt    += len
    return key,cnt
end
function makekey_loop(cnt,body,reqabl) # loop[].(body) reqabl=forloop(range,(rbody))
    range   = reqabl.range
    rbody   = reqabl.body
    tmp,cnt = makekey_tuple(cnt,body,rbody)
    key     = Vector{typeof(tmp)}(undef,range)
    key[1]  = tmp
    for ia  = 2:range
        key[ia],cnt = makekey_tuple(cnt,body,rbody)
    end
    return key,cnt
end

"""
    key = makekey(requested,requestable)

Create a "key" i.e. a data structure of indices into an array `out` of internal
results, returned by the code generated by `@espy`. 

Inputs are
- `requested` a data structure defining a request. This input is provided
   by the user of the code to specify what results are to be extracted.
- `requestable` a named tuple defining the names and sizes of intermediate results
   that can be requested from a given function: this input is provided

# Example

    requestable  = (gp=forloop(2, (z=scalar,s=scalar, material=(a=scalar,b=scalar))),)
    requested    = @request gp[].(s,z,material.(a,b))
    key,nkey     = makekey(requested,requestable)

returns `key` such that

    key.gp[1] == (s=1, z=2, material = (a=3, b=4))
    key.gp[2] == (s=5, z=6, material = (a=7, b=8))
    key.gp[2].material.a == 7
    nkey      == 8

See also: [`@espy`](@ref), [`@espydbg`](@ref), [`@request`](@ref), [`forloop`](@ref), [`scalar`](@ref)
"""
makekey(requested,requestable) = makekey_tuple(0,leftparse(requested),requestable)
"""

    req = @request expr

Create a request of internal results wanted from a function. Considering the function
presented as example for [`@espy`](@ref), examples of possible syntax include

    req       = @request gp[].(s,z,material.(a,b))
    req       = @request gp[].(s)
    req       = @request gp[].(material.(a))

The first expression can be read as follows: "In the function, there is a `for` loop over variable `igp`,
and the results are wanted as a vector (one element for each cycle of the loop).  Each element of the vector
shall be a `type` (a structure) with a field `material`, because a function of that name is called in the
for loop.  Within that function, a variable `a` is to be retrieved.

Note the need to use parentheses also for single-element lists, as in `(s)`.

See also: [`@espy`](@ref), [`@espydbg`](@ref), [`makekey`](@ref)
"""
macro request(ex)
    return QuoteNode(ex)
end

######################## Generate new function code
## Clean code
function code_clean_function(ex::Expr)
    if isexpr(:(.),ex)
        exo = makesub(code_clean_function(getleft(ex)),code_clean_function(getright(ex)))
    else
        exo = Expr(ex.head,[code_clean_function(a) for a ∈ ex.args]...)
    end
    return exo
end
code_clean_function(ex::QuoteNode) = ex.value
code_clean_function(ex)            = ex

#@espy_loop key gp igp → key_gp=key.gp[igp]

#@espy_record out key var tmp → out[key.var] = var
function code_write_to_out(out,key,var)
    return quote
        if haskey($key,$(QuoteNode(var)))                                       # if haskey(key,:x)
            if typeof($key.$var) == Int64
                $out[$key.$var] = $var                                          # out[key.x] = x
            else
                $out[$key.$var] .= $var                                         # out[key.x] = x
            end
        end
    end
end
function code_call_to_function(left,right,out,key,trace=false)
    if @capture(right,  :foo_(args__)   )         # if rhs is call with :foo
        quote
            if haskey($key,$(QuoteNode(foo)))
                $left = $foo($out,$key.$foo,$(args...))  
            else
                $left = $foo($(args...))
            end    
        end 
    else
        quote
            $left = $right
        end 
    end
end

## @espy
function code_extractor_function(ex::Expr,out,key,trace=false)
    return if @capture(ex,    function foo_(args__) body_ end    )                                       # foo(a,b,c)...end
        trace && println("function")
        quote 
            function $foo($out,$key,$(args...)) 
                $(code_extractor_function(body,out,key,trace))
            end
        end                              # foo(out,key,a,b,c)...end
    elseif @capture(ex,    for var_=lo_ : hi_ body_ end   )
        trace && println("for")
        loopname = Symbol(string(var)[2:end])                                   # gp
        subkey   = Symbol(key,"_",loopname)   
        quote
            for $var=$lo:$hi
                if haskey($key,$(QuoteNode(loopname)))                                  # if haskey(key,:gp) QuoteNode to produce a Symbol in the code
                    $subkey   = $key.$loopname[$var]                                     # key_gp=key.gp[igp]
                else
                    $subkey =  Nothing
                end
                $(code_extractor_function(body,out,subkey,trace))
            end
        end
    elseif @capture(ex,  left_ = right_   )
        trace && println("assign")
        if @capture(left,   :name_  )                                                    #:a = ...
            quote
                $(code_call_to_function(name,right,out,key,trace))
                $(code_write_to_out(out,key,name))
            end
        elseif @capture(left,   (args__,)    )                                               # (a,:b) = ...
            postfix  = Vector{Expr}(undef,0)                                          # will contain the macros to insert
            left = ()                                                           # will contain (a,b)
            for arg ∈ args
                if @capture(arg,  :name_  )                                     # ...,:b =
                    left = (left...,name)                                    # (a,b)
                    push!(postfix    ,code_write_to_out(out,key,name))                 # @espy_record out key b
                else                                                            #  a,... =
                    left = (left...,arg)                                          # (a,...)
                end
            end
            quote
             #   $(maketuple(left...)) = $right
                $(code_call_to_function(maketuple(left...),right,out,key,trace))
                $(postfix...)
            end
        else
            code_call_to_function(left,right,out,key,trace)
        end
    else
        trace && println("recursion")
        Expr(ex.head,[code_extractor_function(a,out,key,trace) for a ∈ ex.args]...)
    end
end
function code_extractor_function(ex,out,key,trace=false) # TODO what is this for?
    trace && println("default")
    return ex
end
"""

    @espy function residual(x,y)
        ngp=2
        r = 0
        for igp=1:ngp
            :z = x[igp]+y[igp]
            :s,dum  = :material(z)
            r += s
        end
        return r
    end
    @espy function material(z)
        :a = z+1
        :b = a*z
        return b,3.
    end

Transform the code of a function in which variables and function calls have been annotated with `:`
in order to allow the extraction of intermediate results.

The above annotated code will result in the generation of "clean" code in which the `:` annotations
have been taken out

    function residual(x,y)
        ngp=2
        r = 0
        for igp=1:ngp
            z = x[igp]+y[igp]
            s,dum  = material(z)
            r += s
        end
        return r
    end
    function material(z)
        a = z+1
        b = a*z
        return b,3.
    end

The macro will also generate code with additional `out` and `key` arguments:

    function residual(out,key,x,y)
        ngp = 2
        r   = 0
        for igp = 1:ngp
            @espy_loop key gp                     # key_gp = key.gp[igp]
            z = x[igp]+y[igp]
            @espy_record out key_gp z             # out[key_gp.z] = z
            s = @espy_call out key_gp material(z) # s = material(out,key_gp.material,z)
            @espy_record out key_gp s             # out[key_gp.s] = s
            r += s
        end
        return r
    end
    function material(out,key,z)
        a = z+1
        @espy_record out key a                    # out[key.a] = a
        b = a*z
        @espy_record out key b                    # out[key.b] = b
        return b
    end

The above code contains more macros, which in turn evaluate as
shown in the comments.  More precisely,

    @espy_record out key a

evaluates to

    if haskey(key,a)
        out[key.a] = a
    end

`key` is a data structure generated by [`makekey`](@ref) based on a [`@request`](@ref).

When the version of `residual` with additional parameter `out` has been called, the content
of this output is accessed using `key`:

    requestable  = (gp=forloop(2, (z=scalar,s=scalar, material=(a=scalar,b=scalar))),)
    requested    = @request gp[].(s,z,material.(a,b))
    key,nkey     = makekey(requested,requestable)
    residual(out,key,x,y)
    igp          = 2
    z            = out[key.gp[igp].z]

See also: [`@espydbg`](@ref), [`@request`](@ref), [`makekey`](@ref)
"""
macro espy(ex)
    cleancode = code_clean_function(ex)
    extractorcode  = code_extractor_function(ex,newsym(:out),newsym(:key),false)
    return makeblock(esc(extractorcode),esc(cleancode))
end
"""
    @espydbg function ...
    end

Run [`@espy`](@ref) and to generate code and print the output code (for debug purposes).

See also: [`@espy`](@ref), [`@request`](@ref), [`makekey`](@ref), [`forloop`](@ref), [`scalar`](@ref)"""
macro espydbg(ex)
    println("###### espy clean code")
    cleancode = code_clean_function(ex)
    println(prettify(cleancode))
    println("###### espy extractor code")
    extractorcode  = code_extractor_function(ex,newsym(:out),newsym(:key),false)
    println(prettify(extractorcode))
    println("######")
    return makeblock(esc(extractorcode),esc(cleancode))
end

