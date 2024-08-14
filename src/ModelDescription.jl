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

Provide scaling value for each type (class and field) of dof in the model.  
This is usued to improve the conditioning of the incremental problems and
for convergence criteria. `scale` is a `NamedTuple` with fieldnames within 
`X`, `U` and `A`.  Each field is itself a `NamedTuple` with fieldnames being
dof fields, and value being the expected order of magnitude.

For example `scale = (X=(tx=10,rx=1),A=(drag=3.))` should be read as: X-dofs
of field `:tx` are scaled by 10 moters, `:rx` by radians, and A-dofs of 
field `drag` by 3. All other degrees of freedom are scaled by 1.

Determining scaling coefficients that improve the condition number of
incremental matrices is a hard problem.

`Λscale` is a scalar. The scale of a Λ-dof will be deemed to be the scale of the 
corresponding X-dof, times `Λscale`.

See also: [`addnode!`](@ref), [`describe`](@ref), [`coord`](@ref), [`studyscale`](@ref)
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

Return an initial `State` for the model with all dofs set to zero

Modifying a model (invoquing `addnode!` and `addelement!` after `initialize!` will result in an error) 

Optional keyword arguments:
`nΛder=1`, `nXder=1`, `nUder=1` to specify the number of "derivatives" to store in the `State`.
  note that "`n⋅der==order+1`", that is, for a dynamic problem (with accelerations), `nXder=3` so that
  `order==2`.  Setting  `n⋅der` is only required if [`setdof!`](@ref) will be used.  A dynamic solver
  handles a "static" initial state perfectly well.
`time=-∞` the time associated to the initial state.  Note that for example setting `time=0.`, and then
calling a solver with a first time step also at `0.` causes an error.  

See also: [`setdof!`](@ref), [`Model`](@ref), [`addnode!`](@ref), [`addelement!`](@ref), [`solve`](@ref)
"""
function initialize!(model::Model;nΛder=1,nXder=1,nUder=1,kwargs...)
    assert_unlocked(model)
    model.locked = true
    dis = Disassembler(model)
    return State{nΛder,nXder,nUder}(model,dis;kwargs...)
end
function unlock(model::Model,ID::Symbol)
    ID == model.ID && muscadeerror("ID must be distinct from model.ID")
    newmodel        = deepcopy(model)
    newmodel.locked = false
    newmodel.ID     = ID
    return newmodel
end    
 


 
