using StaticArrays,Printf

# model datastructure - private, structure may change, use accessor functions
"""
    AbstractElement

An abstract data type.  An element type `MyElement` must be
declared as a subtype of `AbstractElement`.

`MyELement`must provide a constructor with interface

    `eleobj = MyElement(nod::Vector{Node}; kwargs...)`

See also: [`coord`](@ref), [`Node`](@ref), [`lagrangian`](@ref)    
"""    
abstract type AbstractElement  end
abstract type ID end 
struct DofID <: ID
    class       :: Symbol      # either :X,:U or :A
    idof        :: 𝕫           # index into model X,U or A
end
struct EleID <: ID
    ieletyp     :: 𝕫      
    iele        :: 𝕫           
end
struct NodID <: ID
    inod        :: 𝕫
end
struct Dof
    ID          :: DofID       # global dof number, unique ID
    nodID       :: NodID       # global number of associated node within class
    idoftyp     :: 𝕫
    eleID       :: Vector{EleID}
end
Dof1 = Vector{Dof}
"""
    Node

The eltype of vectors handed by Muscade as first argument to element constructors.

Example:
`function SingleDofCost(nod::Vector{Node};class::Symbol, ... )`

See also: [`coord`](@ref) 
"""
struct Node 
    ID          :: NodID             # global node number, unique ID
    coord       :: 𝕣1                # all nodes in a model have coordinates in same space
    dofID       :: Vector{DofID}     # list of dofs on this node
    eleID       :: Vector{EleID}     # list of elements connected to the node
end
struct Element
    ID          :: EleID             # global element number, unique ID
    ieletyp     :: 𝕫
    iele        :: 𝕫                 # number of element within type
    nodID       :: Vector{NodID}                
    dofID       :: Vector{DofID}                  
end
mutable struct DofTyp   
    class       :: Symbol            # either :X,:U or :A
    field       :: Symbol            # user defined. e.g. :rx1 :tx1
    scale       :: 𝕣        
    dofID       :: Vector{DofID}    
end
mutable struct Model  
    ID          :: Symbol                  # analyses could have multiple models   
    nod         :: Vector{Node}            # model.nod[nodID]
    ele         :: Vector{Vector{Element}} # model.ele[eleID] or model.ele[ieletyp][iele]
    dof         :: @NamedTuple begin       # model.dof[dofID] or model.dof.X[idof] 
                        X::Dof1
                        U::Dof1
                        A::Dof1
                    end           
    eleobj      :: Vector{Vector{E} where{E<:AbstractElement}}             # model.ele[eleID]or  model.eleobj[ieletyp][iele]
    doftyp      :: Vector{DofTyp}          # model.doftyp[idoftyp]
    scaleΛ      :: 𝕣
    locked      :: 𝕓
end

# Model construction - private

firstindex(x)                     = any(x) ? findfirst(x) : 0
getidoftyp(model,class,field)     = firstindex(doftyp.class==class && doftyp.field==field for doftyp∈model.doftyp)
getidoftyp(model,dofID::DofID)    = model.dof[dofID].idoftyp
getieletyp(model,E)               = firstindex(eltype(eleobj) == E                        for eleobj∈model.eleobj)
getdoftyp(model,args...)          = model.doftyp[getidoftyp(model,args...)]  
getneletyp(model)                 = length(model.ele)

Base.getindex(nod::AbstractArray,nodID::NodID)   = nod[nodID.inod   ]
Base.getindex(dof::NamedTuple{(:X,:U,:A), Tuple{Dof1, Dof1, Dof1}},dofID::DofID)   = dof[dofID.class  ][dofID.idof]
Base.getindex(ele::AbstractArray,eleID::EleID)   = ele[eleID.ieletyp][eleID.iele]
Base.getindex(A  ::AbstractArray,id::AbstractArray{ID})   = [A[i] for i ∈ id]
getidof(E::DataType,class)        = findall(doflist(E).class.==class)  
"""
    getndof(model|Element)
    getndof(model|Element,class)
    getndof(model|Element,(class1,class2,[,...]))

where `class` can be any of `:X`, `:U`, `:A`: get the number of dofs of the
specified dof-classes (default: all classes) for the variable `model` or the type
`Element`.

See also: [`describe`](@ref)
"""
getndof(E::DataType)              = length(doflist(E).inod)
getndof(E::DataType,class)        = length(getidof(E,class))  
getndof(E::DataType,class::Tuple) = ntuple(i->getndof(E,class[i]),length(class))
getndof(model::Model,class::Symbol) = length(model.dof[class])
getndof(model::Model,class::Tuple) = ntuple(i->getndof(model,class[i]),length(class))
getndof(model::Model)             = sum(length(d) for d∈model.dof)
getnele(model::Model,ieletyp)     = length(model.ele[ieletyp])
getnele(model::Model)             = sum(length(e) for e∈model.ele)
newdofID(model::Model,class)      = DofID(class  ,getndof(model,class)+1)
neweleID(model::Model,ieletyp)    = EleID(ieletyp,getndof(model,class)+1)
getnnod(E::DataType)              = maximum(doflist(E).inod) 
getdoflist(E::DataType)           = doflist(E).inod, doflist(E).class, doflist(E).field
function getdoftyp(model::Model,class::Symbol,field::Symbol)
    idoftyp = findfirst(doftyp.class==class && doftyp.field==field for doftyp∈model.doftyp)
    isnothing(idoftyp) && muscadeerror(@sprintf("The model has no dof of class %s and field %s.",class,field))    
    return model.doftyp[idoftyp]    
end
getdofID(model::Model,class::Symbol,field::Symbol) = getdoftyp(model,class,field).dofID
function getdofID(model::Model,class::Symbol,field::Symbol,nodID::AbstractVector{NodID})
    dofID  = getdofID(model,class,field) 
    i      = [model.dof[d].nodID ∈ nodID for d ∈ dofID]
    return dofID[i]
end


# Model construction - API

"""
    model = Model([ID=:my_model])

Construct a blank `model`, which will be mutated to create a finite element [optimization] problem.

See also: [`addnode!`](@ref), [`addelement!`](@ref), [`describe`](@ref), [`solve`](@ref)  
"""
Model(ID=:muscade_model::Symbol) = Model(ID, Vector{Node}(),Vector{Vector{Element}}(),(X=Dof1(),U=Dof1(),A=Dof1()),Vector{Any}(),Vector{DofTyp}(),1.,false)

"""
    nodid = addnode!(model,coord)

If `coord` is an `AbstractVector` of `Real`: add a single node to the model.  
Muscade does not prescribe what coordinate system to use.  Muscade will handle 
to each element the `coord` of the nodes of the element, and the element 
constructor must be able to make sense of it. `nodid` is a node identifier, that
is used as input to `addelement!`.

If `coord` is an `AbstractMatrix`, its rows are treated as vectors of coordinates.
`nodid` is then a vector of node identifiers.

See also: [`addelement!`](@ref), [`coord`](@ref) , [`describe`](@ref) 
"""
function addnode!(model::Model,coord::ℝ2) 
    assert_unlocked(model)
    Δnnod = size(coord,1)
    nodID = [NodID(length(model.nod)+inod) for inod ∈ 1:Δnnod]
    append!(model.nod, [Node(nodID[inod],coord[inod,:],Vector{DofID}(),Vector{EleID}()) for inod = 1:Δnnod] )
    return nodID 
end
addnode!(model::Model,coord::ℝ1)  = addnode!(model,reshape(coord,(1,length(coord))))[1]

Base.adjoint(nodid::NodID)=nodid

"""
    eleid = addelement!(model,ElType,nodid;kwargs...)

Add one or several elements to `model`, connecting them to the nodes specified 
by `nodid`.

If `nodid` is an `AbstractVector`: add a single element to the model. `eleid`
is then a single element identifier.

If `nodid` is an `AbstractMatrix`: add multiple elements to the model. Each 
row of `nodid` identifies the node of a single element. `eleid` is then 
a vector of element identifiers.

For each element, `addelement!` will call `eleobj = ElType(nodes;kwargs...)` where
`nodes` is a vector of nodes of the element.

See also: [`addnode!`](@ref), [`describe`](@ref), [`coord`](@ref)
"""
function addelement!(model::Model,::Type{T},nodID::AbstractMatrix{NodID};kwargs...) where{T<:AbstractElement}
    assert_unlocked(model)
    # new element type? make space in model.eletyp and model.eleobj for that
    nod      = [model.nod[nodID[1,i]] for i∈eachindex(nodID[1,:])]
    ele1     = T(nod;kwargs...)
    E        = typeof(ele1)

    # add eletyp to model (if new)
    ieletyp = getieletyp(model,E)
    if ieletyp == 0 # new element type!
        ieletyp = length(model.ele)+1
        push!(model.eleobj, Vector{E      }()       )
        push!(model.ele   , Vector{Element}()       )
    end        

    iele_sofar       = length(model.ele[ieletyp])
    nele_new,nnod    = size(nodID)
    inod,class,field = getdoflist(E)
    neledof          = getndof(E)
    nnod == getnnod(E) || muscadeerror(@sprintf "Connecting element of type %s: Second dimension of inod (%i) must be equal to element's nnod (%i)" T nnod getnnod(E) ) 
    dofID  = Vector{DofID  }(undef,neledof )         # work array
    eleID  = Vector{EleID  }(undef,nele_new)         # allocate return variable
    ele    = Vector{Element}(undef,nele_new)         # work array - will be appended to model.ele
    eleobj = Vector{E      }(undef,nele_new)         # work array - will be appended to model.eleobj
    for iele_new = 1:nele_new
        # add eleID to nodes
        eleID[iele_new] = EleID(ieletyp,iele_new+iele_sofar) 
        for nod ∈ [model.nod[i] for i∈unique(nodID[iele_new,:])] # unique: if several nodes of an element are connected to the same model node, mention element only once
            push!(nod.eleID,eleID[iele_new])
        end
        # for all dofs of the element, make sure they are represented in the nod object and get an ID
        for ieledof  = 1:neledof
            nodid    = nodID[iele_new,inod[ieledof]]  # nodID of current eledof
            idoftyp  = getidoftyp(model,class[ieledof],field[ieledof])
            # add doftyp to model (if new)  
            if idoftyp == 0 # new dof type
                push!(model.doftyp, DofTyp(class[ieledof],field[ieledof],1.,DofID[])) # DofID[]: do not add dof to doftyp, though
                idoftyp = length(model.doftyp)
            end
            # add dof to model (if new)
            idofID = firstindex(model.dof[dofID].idoftyp == idoftyp for dofID ∈ model.nod[nodid].dofID) 
            if idofID == 0 # new dof
                dofID[ieledof] = newdofID(model,class[ieledof]) 
                push!(model.dof[class[ieledof]], Dof(dofID[ieledof],nodid,idoftyp,𝕫[]) ) # do not add element to dof, though
                push!(model.doftyp[idoftyp].dofID, dofID[ieledof])
            else
                dofID[ieledof] = model.nod[nodid].dofID[idofID]
            end
            # add element to dof (always)
            push!(model.dof[dofID[ieledof]].eleID, eleID[iele_new])  
            # add dof to node (if new)
            if dofID[ieledof] ∉ model.nod[nodid].dofID
                push!(model.nod[nodid].dofID, dofID[ieledof])
            end
        end
        # add element to model (always)
        ele[   iele_new] = Element(eleID[iele_new], ieletyp, iele_sofar+iele_new, nodID[iele_new,:],deepcopy(dofID))
        if iele_new==1
            eleobj[iele_new] = ele1
        else
            nod = [model.nod[nodID[iele_new,i]] for i∈eachindex(nodID[iele_new,:])]
            eleobj[iele_new] = T(nod;kwargs...)
        end
    end
    append!(model.ele[   ieletyp],ele   )
    append!(model.eleobj[ieletyp],eleobj)
    return eleID 
end
addelement!(model::Model,::Type{E},nodID::AbstractVector{NodID};kwargs...) where{E<:AbstractElement} = addelement!(model,E,reshape(nodID,(1,length(nodID)));kwargs...)[1] 

"""
    setscale!(model;scale=nothing,Λscale=nothing)

Provides an order of magnitude for each type of dof in the model.  
This is usued to improve the conditioning of the incremental problems and
for convergence criteria. `scale` is a `NamedTuple` with fieldnames within 
`X`, `U` and `A`.  Each field is itself a `NamedTuple` with fieldnames being
dof fields, and value being the expected order of magnitude.

For example `scale = (X=(tx=10,rx=1),A=(drag=3.))` should be read as: X-dofs
of field `:tx` are in tens of meters, `:rx` are in radians, and A-dofs of 
field `drag` are of the order of 3. All other degrees of freedom are expected 
to be roughly of the order of 1.

`Λscale` is a scalar. The scale of a Λ-dof will be deemed to be the scale of the 
corresponding X-dof, times `Λscale`.

See also: [`addnode!`](@ref), [`describe`](@ref), [`coord`](@ref)
"""
function setscale!(model::Model;scale=nothing,Λscale=nothing)  # scale = (X=(tx=10,rx=1),A=(drag=3.))
    assert_unlocked(model)
    if ~isnothing(scale)
        for doftyp ∈ model.doftyp
            if doftyp.class ∈ keys(scale) && doftyp.field ∈ keys(scale[doftyp.class])
                doftyp.scale = scale[doftyp.class][doftyp.field] # otherwise leave untouched
            end
        end
    end
    if ~isnothing(Λscale)
        model.scaleΛ = Λscale
    end
end

assert_unlocked(model::Model) = model.locked && muscadeerror(@sprintf("model %s is initialized and can no longer be edited",model.ID))
"""
    initialstate = initialize!(model)

Finalize a model (invoquing `addnode!` and `addelement!` after `initialize!` will result in an error) 
and return an initial `State` (with all dofs set to zero, as starting point for an analysis.)

See also: [`addnode!`](@ref), [`addelement!`](@ref), [`solve`](@ref)
"""
function initialize!(model::Model)
    assert_unlocked(model)
    model.locked = true
    return State(model,Disassembler(model))
end
function unlock(model::Model,ID::Symbol)
    ID == model.ID && muscadeerror("ID must be distinct from model.ID")
    newmodel        = deepcopy(model)
    newmodel.locked = false
    newmodel.ID     = ID
    return newmodel
end    

### extracting data from the model
eletyp(model::Model) = eltype.(model.eleobj) 

### getting printout of the model
using Printf
"""
    describe(model,spec)

Print out information about `model`.
`spec` can be 
- an `EleID` to describe an element,
- a `DofID` to describe a dof.
- a `NodID` to describe a node,
- `:doftyp` to obtain a list of doftypes, 
- `:dof` to obtain a list of dofs or 
- `:eletyp` for a list of element types.

See also: [`addelement!`](@ref), [`addnode!`](@ref)
"""
function describe(model::Model,eleID::EleID)
    try 
        dof = model.ele[eleID] 
    catch
        printstyled("Not a valid EleID\n",color=:red,bold=true)
        return
    end
    e  = model.ele[eleID]
    eo = model.eleobj[eleID]
    @printf "Element with EleID(%i,%i)\n" eleID.ieletyp eleID.iele 
    @printf "   model.eleobj[%i][%i]::" eleID.ieletyp eleID.iele 
    printstyled(@sprintf("%s\n",typeof(eo)),color=:cyan)
    @printf "   model.ele[%i][%i]:\n" eleID.ieletyp eleID.iele
    for dofid ∈ e.dofID
        dof    = model.dof[dofid]
        nod    = model.nod[dof.nodID]
        doftyp = model.doftyp[dof.idoftyp]
        @printf "      NodID(%i), class=:%s, field=:%-12s\n" dof.nodID.inod doftyp.class doftyp.field 
    end
end
function describe(model::Model,dofID::DofID)
    try 
        dof = model.dof[dofID] 
    catch
        printstyled("Not a valid DofID\n",color=:red,bold=true)
        if dofID.class==:Λ
            @printf "Optimisation solvers introduce a one-to-one correspondance between :Λ-dofs and :X-dofs, \nbut :Λ-dofs are not part of the model description: try DofID(:X,...)\n"
        end
        return
    end
    dof     = model.dof[dofID] 
    doftyp  = model.doftyp[dof.idoftyp]
    @printf "Degree of freedom with DofID(:%s,%i)\n" dofID.class dofID.idof
    @printf "   model.dof.%s[%i]:\n" dofID.class dofID.idof
    @printf "   NodID(%i), class=:%s, field=:%-12s\n" dof.nodID.inod dofID.class doftyp.field 
    @printf "   elements:\n"
    for eleid ∈ dof.eleID
        @printf "      EleID(%i,%i), " eleid.ieletyp eleid.iele 
        printstyled(@sprintf("%s\n",eltype(model.eleobj[eleid.ieletyp])),color=:cyan)
    end
    if dofID.class == :X
        @printf "   Output in state[istep].X[ider+1][%i] and state[istep].Λ[%i]\n" dofID.idof dofID.idof    
    elseif dofID.class ==:U
            @printf "   Output in state[istep].U[ider][%i]\n" dofID.idof   
        elseif dofID.class == :A
        @printf "   Output in state[istep].A[%i]\n" dofID.idof   
    end            
end
function describe(model::Model,nodID::NodID)
    try 
        nod = model.nod[nodID] 
    catch
        printstyled("Not a valid NodID\n",color=:red,bold=true)
        return
    end
    nod = model.nod[nodID]
    @printf "Node with NodID(%i)\n" nodID.inod
    @printf "   model.nod[%i]:\n" nodID.inod
    nc = length(nod.coord)
    @printf "   coord=[" 
    for ic=1:nc-1
        @printf "%g," nod.coord[ic] 
    end
    if nc>0
        @printf "%g" nod.coord[nc] 
    end
    @printf "]\n" 
    @printf "   dof (degrees of freedom):\n"
    for dofid ∈ nod.dofID
        dof = model.dof[dofid]
        doftyp = model.doftyp[dof.idoftyp]
        @printf "      DofID(:%s,%i), class=:%s, inod=%i, field=:%-12s\n" dofid.class dofid.idof dofid.class nodID.inod doftyp.field    
    end
    @printf "   elements:\n"
    for eleID ∈ nod.eleID
        @printf "      EleID(%i,%i), " eleID.ieletyp eleID.iele 
        printstyled(@sprintf("%s\n",typeof(model.eleobj[eleID])),color=:cyan)
    end
 end
function describe(model::Model,s::Symbol)
    @printf "\nModel '%s'\n" model.ID
    if s==:dof
        for class ∈ (:X,:U,:A)
            ndof    = getndof(model,class)
            ndof>0 && @printf "\n   Dofs of class :%s\n" class
            for idof = 1:ndof
                dof     = model.dof[class][idof] 
                doftyp  = model.doftyp[dof.idoftyp]
                @printf "      %6d. field= :%-15s NodID(%i)\n" idof doftyp.field dof.nodID.inod 
            end
        end
    elseif s==:doftyp
        for doftyp ∈ model.doftyp
            doftyp.class==:X && @printf "      class= :%-5s field= :%-15s\n" doftyp.class doftyp.field 
        end
        for doftyp ∈ model.doftyp
            doftyp.class==:U && @printf "      class= :%-5s field= :%-15s\n" doftyp.class doftyp.field 
        end
        for doftyp ∈ model.doftyp
            doftyp.class==:A && @printf "      class= :%-5s field= :%-15s\n" doftyp.class doftyp.field 
        end
    elseif s==:eletyp
        et = eletyp(model)
        for i∈eachindex(et)
            @printf "    %6d elements of type %s\n" length(model.eleobj[i]) et[i]
        end
    end
end
 
