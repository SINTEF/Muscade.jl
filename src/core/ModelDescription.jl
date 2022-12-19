using StaticArrays,Printf

# model datastructure - private, structure may change, use accessor functions

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
Base.getindex(ele::AbstractArray,eleID::EleID)   = ele[eleID.ieletyp][eleID.iele]
Base.getindex(A  ::AbstractArray,id::AbstractArray{ID})   = [A[i] for i ∈ id]
getndof(model::Model,class::Symbol) = length(model.dof[class])
getndof(model::Model,class::Tuple) = (getndof(model,c) for c∈class)
getndof(model::Model)             = sum(length(d) for d∈model.dof)
getnele(model::Model,ieletyp)     = length(model.ele[ieletyp])
getnele(model::Model)             = sum(length(e) for e∈model.ele)
newdofID(model::Model,class)      = DofID(class  ,getndof(model,class)+1)
neweleID(model::Model,ieletyp)    = EleID(ieletyp,getndof(model,class)+1)
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
        ele[   iele_new] = Element(eleID[iele_new], ieletyp, iele_sofar+iele_new, nodID[iele_new,:],dofID)
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

### Obtain printouts describing elements, nodes or dofs of a model
using Printf
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
        @printf "      DofID(:%s,%i), class=:%s, idof=%i, field=:%-12s\n" dofid.class dofid.idof dofid.class dofid.idof doftyp.field    
    end
    @printf "   elements:\n"
    for eleID ∈ nod.eleID
        @printf "      EleID(%i,%i), " eleID.ieletyp eleID.iele 
        printstyled(@sprintf("%s\n",typeof(model.eleobj[eleID])),color=:cyan)
    end
 end
