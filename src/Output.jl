## Nodal results
"""
    dofres = getdof(state;[class=:X],field=:somefield,nodID=[nodids...],[order=0])

Obtain the value of dofs of the same class and field, at various nodes and for various states.

If `state` is a vector, the output `dofres` has size `(ndof,nstate)`.
If `state` is a scalar, the output `dofres` has size `(ndof,)`.

See also: [`getresult`](@ref), [`addnode!`](@ref), [`solve`](@ref)
"""
function getdof(state::State;kwargs...)  
    dofres = getdof([state];kwargs...)
    return reshape(dofres,size(dofres)[1]) 
end
function getdof(state::Vector{S};class::Symbol=:X,field::Symbol,nodID::Vector{NodID}=NodID[],order::ℤ=0)where {S<:State}
        class ∈ [:Λ,:X,:U,:A] || muscadeerror(sprintf("Unknown dof class %s",class))
    c       = class==:Λ      ? :X                                   : class
    dofID   = nodID==NodID[] ? getdofID(state[begin].model,c,field) : getdofID(state[begin].model,c,field,nodID)
    dofres  = 𝕣2(undef,length(dofID),length(state)) 
    for istate ∈ eachindex(state)
        sc = if class==:Λ state[istate].Λ
        elseif  class==:X state[istate].X
        elseif  class==:U state[istate].U
        elseif  class==:A state[istate].A
        end
        if class == :A
            for (idof,d) ∈ enumerate(dofID)
                dofres[idof,istate] = sc[d.idof] 
            end
        else
            if order+1 ≤ length(sc)
                s = sc[order+1]
                for (idof,d) ∈ enumerate(dofID)
                    dofres[idof,istate] = s[d.idof] 
                end
            else
                for (idof,d) ∈ enumerate(dofID)
                    dofres[idof,istate] = 0. 
                end
            end
        end
    end
    return dofres
end
"""
    state = setdof!(state,value        ;[class=:X],field=:somefield,                  [order=0])
    state = setdof!(state,value::Vector;[class=:X],field=:somefield,nodID=[nodids...],[order=0])

Set the value of dofs of the same class and field, at various nodes and for various states.
There are two methods:
1. A single `value` is applied to all relevant nodes in the model
2. `value` and `nodID` are vectors of the same lengths, and each element in `value` is applied to the corresponding node.

`setdof!` is peculiar in that it modifies its input `state` variable, but must be used as a function. A call like
`stateout = setdof!(statein,value;class=:X,field=:somefield,order=1)`
can turn out in two ways:
If the `state` already stores derivatives in `X` to order 1, then `statein` is mutated and `statein===stateout`.
Otherwise, `statein` is unchanged, `stateout` is a new object, sharing as much memory as possible with `statein`.
To avoid confusion, always use the syntax shown above.

See also: [`getresult`](@ref), [`addnode!`](@ref), [`solve`](@ref)
"""
function setdof!(state::State,dofval::𝕣1;class::Symbol=:X,field::Symbol,nodID::Vector{NodID},order::ℤ=0)
    class ∈ [:Λ,:X,:U,:A] || muscadeerror(sprintf("Unknown dof class %s",class))
    c     = class==:Λ ? :X : class
    dofID = getdofID(state.model,c,field,nodID)
    if class == :A 
        for (idof,d) ∈ enumerate(dofID)
            state.A[d.idof] = dofval[idof]  
        end
    else
        sc = if class==:Λ state.Λ
        elseif  class==:X state.X
        elseif  class==:U state.U
        end
        nder = length(sc)
        if order+1>nder
            sc    = (sc...,ntuple(i->zeros(getndof(state.model,class)),order+1-nder)...)
            state = if class==:Λ State(state.time,sc     ,state.X,state.U,state.A,state.SP,state.model,state.dis)
            elseif     class==:X State(state.time,state.Λ,sc     ,state.U,state.A,state.SP,state.model,state.dis)
            elseif     class==:U State(state.time,state.Λ,state.X,sc     ,state.A,state.SP,state.model,state.dis)
            end
        end
        s = sc[order+1] 
        for (idof,d) ∈ enumerate(dofID)
            s[d.idof] = dofval[idof]  
        end
    end
    return state
end
function setdof!(state::State,dofval::𝕣;class::Symbol=:X,field::Symbol,order::ℤ=0)
    class ∈ [:Λ,:X,:U,:A] || muscadeerror(sprintf("Unknown dof class %s",class))
    c     = class==:Λ ? :X : class
    dofID = getdofID(state.model,c,field)
    if class == :A 
        for (idof,d) ∈ enumerate(dofID)
            state.A[d.idof] = dofval[idof]  
        end
    else
        sc = if class==:Λ state.Λ
        elseif  class==:X state.X
        elseif  class==:U state.U
        end
        nder = length(sc)
        if order+1>nder
            sc    = (sc...,ntuple(i->zeros(getndof(state.model,class)),order+1-nder)...)
            state = if class==:Λ State(state.time,sc     ,state.X,state.U,state.A,state.SP,state.model,state.dis)
            elseif     class==:X State(state.time,state.Λ,sc     ,state.U,state.A,state.SP,state.model,state.dis)
            elseif     class==:U State(state.time,state.Λ,state.X,sc     ,state.A,state.SP,state.model,state.dis)
            end
        end
        s = sc[order+1] 
        for (idof,d) ∈ enumerate(dofID)
            s[d.idof] = dofval  
        end
    end
    return state
end
# Elemental results

function extractkernel!(iele::AbstractVector{𝕫},eleobj::Vector{E},dis::EletypDisassembler,state::Vector{S},dbg,req) where{E,S<:State}# typestable kernel
    return [begin
        index = dis.index[i]
        Λ     = s.Λ[1][index.X]                 
        X     = Tuple(Xᵢ[index.X] for Xᵢ∈s.X)
        U     = Tuple(Uᵢ[index.U] for Uᵢ∈s.U)
        A     = s.A[index.A]
        L,FB,e = getlagrangian(eleobj[i],Λ,X,U,A,s.time,s.SP,(dbg...,istep=istep,iele=i),req)
        e
    end for i∈iele, (istep,s)∈enumerate(state)]
end
"""
    eleres = getresult(state,req,els)

Obtain an array of nested `NamedTuples` and `NTuples` of element results.
`req` is a request defined using `@request`.
`state` a vector of `State`s or a single `State`.
`els` can be either
- a vector of `EleID`s (obtained from `addelement!`) all corresponding
  to the same concrete element type
- a concrete element type (see [`describe`](@ref)).

If `state` is a vector, the output `dofres` has size `(nele,nstate)`.
If `state` is a scalar, the output `dofres` has size `(nele)`.

See also: [`getdof`](@ref), [`@request`](@ref), [`@espy`](@ref), [`addelement!`](@ref), [`solve`](@ref), [`describe`](@ref)
"""
function getresult(state::Vector{S},req,eleID::Vector{EleID})where {S<:State}
    # Some elements all of same type, multisteps
    # eleres[iele,istep].gp[3].σ
    ieletyp             = eleID[begin].ieletyp
    all(e.ieletyp== ieletyp for e∈eleID) || muscadeerror("All elements must be of the same element type")
    eleobj              = state[begin].model.eleobj[ieletyp]
    dis                 = state[begin].dis.dis[ieletyp]
    iele                = [e.iele for e∈eleID]
    return extractkernel!(iele,eleobj,dis,state,(func=:getresult,ieletyp=ieletyp),req)
end

function getresult(state::Vector{S},req,::Type{E}) where{S<:State,E<:AbstractElement}
    # All elements within the type, multisteps
    # eleres[iele,istep].gp[3].σ
    ieletyp = findfirst(E.==eletyp(state[begin].model))
    isnothing(ieletyp) && muscadeerror("This type of element is not in the model.")
    eleobj              = state[begin].model.eleobj[ieletyp]
    dis                 = state[begin].dis.dis[ieletyp]
    iele                = eachindex(eleobj)
    return extractkernel!(iele,eleobj,dis,state,(func=:getresult,eletyp=E),req)
end    
# single step
# eleres[iele].gp[3].σ
getresult(state::State,req,args...) = flat(getresult([state],req,args...)) 

"""
    ilast = findlastassigned(state)

Find the index `ilast` of the element before the first non assigment element in a vector `state`.

In multistep analyses, `solve` returns a vector `state` of length equal to the number of steps
requested by the user.  If the analysis is aborted, `solve` still returns any available results
at the begining of `state`, and the vector `state[1:ilast]` is fully assigned.

See also: [`solve`](@ref)
"""     
findlastassigned(v::Vector) = findlast([isassigned(v,i) for i=1:length(v)])


