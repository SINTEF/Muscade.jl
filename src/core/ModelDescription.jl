using StaticArrays,Printf

# model datastructure - private, structure may change, use accessor functions

# TODO REFACTOR: model.ele and model.eleobj are both Vector{Vector} and can both be accessed directly or using model.ele[eleID]

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
struct Node 
    ID          :: NodID              # global node number, unique ID
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
    eleobj      :: Vector{Any}             # model.ele[eleID]or  model.eleobj[ieletyp][iele]
    doftyp      :: Vector{DofTyp}          # model.doftyp[idoftyp]
    Λscale      :: 𝕣
end

# Model construction - private

firstindex(x)                     = any(x) ? findfirst(x) : 0
getidoftyp(model,class,field)     = firstindex(doftyp.class==class && doftyp.field==field for doftyp∈model.doftyp)
getidoftyp(model,dofID::DofID)    = model.dof[dofID].idoftyp
getieletyp(model,E)               = firstindex(eltype(eleobj) == E                        for eleobj∈model.eleobj)
getdoftyp(model,args...)          = model.doftyp[getidoftyp(model,args...)]  

Base.getindex(nod::AbstractArray,nodID::NodID)   = nod[nodID.inod   ]
Base.getindex(dof::NamedTuple{(:X,:U,:A), Tuple{Dof1, Dof1, Dof1}},dofID::DofID)   = dof[dofID.class  ][dofID.idof]
Base.getindex(ele::AbstractArray,eleID::EleID)   = ele[eleID.ieletyp][dofID.iele]
Base.getindex(A  ::AbstractArray,id::AbstractArray{ID})   = [A[i] for i ∈ id]
getndof(model::Model,class)       = length(model.dof[class])
getndof(model::Model)             = sum(length(d) for d∈model.dof)
getnele(model::Model,ieletyp)     = length(model.ele[ieletyp])
getnele(model::Model)             = sum(length(e) for e∈model.ele)
newdofID(model::Model,class)      = DofID(class  ,getndof(model,class)+1)
neweleID(model::Model,ieletyp)    = EleID(ieletyp,getndof(model,class)+1)

# Model construction - API

# TODO sizehint!(vec,n) to accelerate push! .

Model(ID=:muscade_model::Symbol) = Model(ID, Vector{Node}(),Vector{Vector{Element}}(),(X=Dof1(),U=Dof1(),A=Dof1()),Vector{Any}(),Vector{DofTyp}(),1.)

function addnode!(model::Model,coord::ℝ2) 
    Δnnod = size(coord,1)
    nodID = [NodID(length(model.nod)+inod) for inod ∈ 1:Δnnod]
    append!(model.nod, [Node(nodID[inod],coord[inod,:],Vector{DofID}(),Vector{EleID}()) for inod = 1:Δnnod] )
    return nodID 
end
addnode!(model::Model,coord::ℝ1)  = addnode!(model,reshape(coord,(1,length(coord))))[1]

coord(nod::AbstractVector{Node}) = [n.coord for n∈nod]

function addelement!(model::Model,::Type{T},nodID::Matrix{NodID};kwargs...) where{T<:AbstractElement}
    # new element type? make space in model.eletyp and model.eleobj for that
#    ele1     = T(collect(model.nod[nodID[1,:]]);kwargs...)
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
        nod = [model.nod[i] for i∈unique(nodID[iele_new,:])]
        for nod ∈ [model.nod[i] for i∈unique(nodID[iele_new,:])] # unique: if several nodes of an element are connected to the same model node, mention element only once
            push!(nod.eleID,eleID[iele_new])
        end
        # for nod ∈ model.nod[unique(nodID[iele_new,:])] # unique: if several nodes of an element are connected to the same model node, mention element only once
        #     push!(nod.eleID,eleID[iele_new])
        # end
        # for all dofs of the element, make sure they are represented in the nod object and get an ID
        for ieledof  = 1:neledof
            nodid    = nodID[iele_new,inod[ieledof]]  # nodID of current eledof
            idoftyp  = getidoftyp(model,class[ieledof],field[ieledof])
            # add doftyp to model (if new)  
            if idoftyp == 0 # new!
                idoftyp = length(model.doftyp)+1
                push!(model.doftyp, DofTyp(class[ieledof],field[ieledof],1.,DofID[])) # DofID[]: do not add dof to doftyp, though
            end
            # add dof to model (if new)
            idofID = firstindex(model.dof[dofID].idoftyp == idoftyp for dofID ∈ model.nod[nodid].dofID) 
            if idofID == 0 # new
                dofID[ieledof] = newdofID(model,class[ieledof]) 
                push!(model.dof[class[ieledof]], Dof(dofID[ieledof],nodid,idoftyp,𝕫[]) ) # do not add element to dof, though
            else
                dofID[ieledof] = model.nod[nodid].dofID[idofID]
            end
            # add dof to doftyp (always)
            push!(model.doftyp[idoftyp].dofID, dofID[ieledof])
            # add element to dof (always)
            push!(model.dof[dofID[ieledof]].eleID, eleID[iele_new])  
            # add dof to node (if new)
            if dofID[ieledof] ∉ model.nod[nodid].dofID
                push!(model.nod[nodid].dofID, dofID[ieledof])
            end
        end
        # add element to model (always)
        ele[   iele_new] = Element(eleID[iele_new], ieletyp, iele_sofar+iele_new, nodID[iele_new,:],dofID)
#        eleobj[iele_new] = iele_new==1 ? ele1 : T(collect(model.nod[nodID[iele_new,:]]);kwargs...)   # call element constructor
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
addelement!(model::Model,::Type{E},nodID::Vector{NodID};kwargs...) where{E<:AbstractElement} = addelement!(model,E,reshape(nodID,(1,length(nodID)));kwargs...)[1] 

function setscale!(model;scale=nothing,Λscale=nothing)  # scale = (X=(tx=10,rx=1),A=(drag=3.))
    if ~isnothing(scale)
        for doftyp ∈ model.doftyp
            if doftyp.class ∈ keys(scale) && doftyp.field ∈ keys(scale[doftyp.class])
                doftyp.scale = scale[doftyp.class][doftyp.field] # otherwise leave untouched
            end
        end
    end
    if ~isnothing(Λscale)
        model.Λscale = Λscale
    end
end
