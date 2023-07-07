## Nodal results
"""
    dofres,dofID = getdof(state;[class=:X],field=:somefield,nodID=[nodids...],[iders=0|1|2])

Obtain the value of dofs of the same class and field, at various nodes and for various states.

If `state` is a vector, the output `dofres` has size `(ndof,nder+1,nstate)`.
If `state` is a scalar, the output `dofres` has size `(ndof,nder+1)`.

See also: [`getresult`](@ref), [`addnode!`](@ref), [`solve`](@ref)
"""
function getdof(state::State;kwargs...)  
    dofres,dofID = getdof([state];kwargs...)
    return reshape(dofres,size(dofres)[1:2]),dofID 
end
function getdof(state::Vector{S};class::Symbol=:X,field::Symbol,nodID::Vector{NodID}=NodID[],iders::ℤ1=[0])where {S<:State}
    class ∈ [:Λ,:X,:U,:A] || muscadeerror(sprintf("Unknown dof class %s",class))
    c     = class==:Λ      ? :X                                   : class
    dofID = nodID==NodID[] ? getdofID(state[begin].model,c,field) : getdofID(state[begin].model,c,field,nodID)
    iders = class∈[:Λ,:A]  ? [0]                                  : iders
    dofres   = Array{𝕣,3}(undef,length(dofID),length(iders),length(state)) # dofres[inod,ider+1]
    for istate ∈ eachindex(state)
        for ider∈iders
            s = if class==:Λ; state[istate].Λ 
            elseif class==:X; state[istate].X[ider+1]    
            elseif class==:U; state[istate].U[ider+1]    
            elseif class==:A; state[istate].A    
            end
            for (idof,d) ∈ enumerate(dofID)
                dofres[idof,ider+1,istate] = s[d.idof] 
            end
        end 
    end
    return dofres,dofID
end

# Elemental results

function extractkernel!(iele::AbstractVector{𝕫},eleobj::Vector{E},dis::EletypDisassembler,state::Vector{S},dbg,req) where{E,S<:State}# typestable kernel
    return [begin
        index = dis.index[i]
        Λ     = s.Λ[index.X]                 
        X     = Tuple(x[index.X] for x∈s.X)
        U     = Tuple(u[index.U] for u∈s.U)
        A     = s.A[index.A]
        L,χn,FB,e = getlagrangian(implemented(eleobj[i])...,eleobj[i],Λ,X,U,A,s.time,nothing,nothing,s.SP,(dbg...,istep=istep,iele=i),req)
        e
    end for i∈iele, (istep,s)∈enumerate(state)]
end
"""
    eleres = getresult(state,req,eleids)

Obtain an array of nested NamedTuples and NTuples of element results.
`req` is a request defined using `@request`.
`state` a vector of `State`s or a single `State`.
`eleids` can be either
- a vector of `EleID`s (obtained from `addelement!`) all corresponding
  to the same concrete element type
- a concrete element type.

If `state` is a vector, the output `dofres` has size `(nele,nstate)`.
If `state` is a scalar, the output `dofres` has size `(nele)`.

See also: [`getdof`](@ref), [`@request`](@ref), [`@espy`](@ref), [`addelement!`](@ref), [`solve`](@ref)
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
    isnothing(ieletyp) && muscadeerror("This type of element is not in the model. See 'eletyp(model)'")
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

############## describe state to the user
function describeX(state::State)
    model = state.model
    nX    = getndof(model,:X)
    nder  = length(state.X)
    for iX = 1:nX
        dofID   = DofID(:X,iX)
        dof     = model.dof[dofID] 
        doftyp  = model.doftyp[dof.idoftyp]
        @printf "NodID(%i), class=:%s, field=:%-15s   " dof.nodID.inod dofID.class doftyp.field
        for ider = 1:nder
            @printf "%15g " state.X[ider][iX]
        end
        @printf "\n" 
    end
end
function describeΛX(state::State)
    model = state.model
    nX    = getndof(model,:X)
    nder  = length(state.X)
    for iX = 1:nX
        dofID   = DofID(:X,iX)
        dof     = model.dof[dofID] 
        doftyp  = model.doftyp[dof.idoftyp]
        @printf "NodID(%i), class=:%s, field=:%-15s   %15g " dof.nodID.inod dofID.class doftyp.field state.Λ[iX]
        for ider = 1:nder
            @printf "%15g " state.X[ider][iX]
        end
        @printf "\n" 
    end
end
function describeU(state::State)
    model = state.model
    nU    = getndof(model,:U)
    nder  = length(state.U)
    for iU = 1:nU
        dofID   = DofID(:U,iU)
        dof     = model.dof[dofID] 
        doftyp  = model.doftyp[dof.idoftyp]
        @printf "NodID(%i), class=:%s, field=:%-15s   " dof.nodID.inod dofID.class doftyp.field
        for ider = 1:nder
            @printf "%15g " state.U[ider][iU]
        end
        @printf "\n"
    end
end
function describeA(state::State)
    model = state.model
    nA    = getndof(model,:A)
    for iA = 1:nA
        dofID   = DofID(:A,iA)
        dof     = model.dof[dofID] 
        doftyp  = model.doftyp[dof.idoftyp]
        @printf "NodID(%i), class=:%s, field=:%-15s   %15g\n" dof.nodID.inod dofID.class doftyp.field state.A[iA] 
    end
end
function describeScale(state::State)
    model = state.model
    scale = Dict{Symbol,Dict{Symbol,𝕣}}() # scale[class][field]
    for class ∈ (:Λ,:X,:U,:A)
        scale[class] = Dict{Symbol,𝕣}()
        for dof ∈ model.dof[class==:Λ ? :X : class]
            idof    = dof.ID.idof
            field   = model.doftyp[dof.idoftyp].field
            val     = if class==:Λ; state.Λ[idof] 
                  elseif class==:X; state.X[1][idof]    
                  elseif class==:U; state.U[1][idof]    
                  elseif class==:A; state.A[idof]    
            end
            top     = max(get(scale,field,0.),abs(val))
            scale[class][field] = top
        end
        for field ∈ keys(scale[class])
            @printf "class= :%-5s field= :%-15s  max(|dof|)= %g\n" class field scale[class][field]
        end
    end
end    

"""
    describe(state;class=:all)

Provide a description of the dofs stored in `state`.
`class` can be either `:all`, `:Λ`, `:ΛX`, `:X`, `:U`, `:A` or `:scale`

See also: [`solve`](@ref)
"""
function describe(state::State;class::Symbol=:all)
    if class ==:all
        describeΛX(state)
        describeU(state)
        describeA(state)
    elseif class==:Λ || class==:ΛX
        describeΛX(state)
    elseif class ==:X    
        describeX(state)
    elseif class ==:U    
        describeU(state)
    elseif class ==:A    
        describeA(state)
    elseif class == :scale 
        describeScale(state)    
    else
        printstyled("Not a valid class\n",color=:red,bold=true)
    end
end

