using StaticArrays,Printf

# model datastructure - private, structure may change, use accessor functions

struct Dof
    ID          :: 𝕫       # global dof number, unique ID
    nodID       :: 𝕫       # global number of associated node
    doftypID    :: 𝕫
    eleID       :: 𝕫1
end
struct Node 
    ID          :: 𝕫                 # global node number, unique ID
    coord       :: 𝕣1                # all nodes in a model have coordinates in same space
    dofID       :: 𝕫1                # list of dofs on this node
    eleID       :: 𝕫1                # list of elements connected to the node
end
struct Element
    ID          :: 𝕫                 # global element number, unique ID
    eletypID    :: 𝕫
    iele        :: 𝕫                 # number of element within type
    nodID       :: 𝕫1                
    dofID       :: 𝕫1                  
end
struct DofTyp   
    ID          :: 𝕫
    class       :: Symbol  # either :X,:U or :A
    field       :: Symbol  # user defined. e.g. :rx1 :tx1
    scale       :: 𝕣        
    dofID       :: 𝕫1    
end
struct EleTyp
    ID          :: 𝕫
    type        :: DataType
end
mutable struct Model  
    ID          :: Symbol                # analyses could have multiple models   
    dof         :: Vector{Dof}           # model.dof[dofID]
    nod         :: Vector{Node}          # model.nod[nodID]
    ele         :: Vector{Element}       # model.ele[eleID]
    doftyp      :: Vector{DofTyp}        # model.doftyp[doftypID]
    eletyp      :: Vector{EleTyp}        # model.eletyp[eletypID]
    eleobj      :: Vector{Any}           # el = model.eleobj[eletypID][iele]
    Λscale      :: 𝕣
end

# Model construction - private

maxdofID(model)                   = length(model.dof)
maxnodID(model)                   = length(model.nod)
maxeleID(model)                   = length(model.ele)
maxdoftypID(model)                = length(model.doftyp)
maxeletypID(model)                = length(model.eletyp)

firstindex(x)                     = any(x) ? findfirst(x) : 0
safeindirection(a,i)              = i==0 ? zero(eltype(a)) : a[i]
getdoftypID(model,class,field)    = firstindex(doftyp.class==class && doftyp.field==field for doftyp∈model.doftyp)
geteletypID(model,E)              = firstindex(eletyp.type ==E                            for eletyp∈model.eletyp)
getdoftyp(model,class,field)      = model.doftyp[getdoftypID(model,class,field)]  
getdofID(model,nodID,class,field) = safeindirection(model.nod[nodID].dofID , 
                                    firstindex(model.dof[dofID].doftypID == getdoftypID(model,class,field) for dofID ∈ model.nod[nodID].dofID)) 
# Model construction - API

Model(ID=:muscade_model::Symbol) = Model(ID, Vector{Dof}(),Vector{Node}(),Vector{Element}(),Vector{DofTyp}(),Vector{EleTyp}(),Vector{Any}(),1.)

function addnode!(model::Model,coord::ℝ2) 
    Δnnod = size(coord,1)
    nodID = maxnodID(model).+(1:Δnnod)
    append!(model.nod, [Node(maxnodID(model)+inod,coord[inod,:],Vector{𝕫}(),Vector{𝕫}()) for inod = 1:Δnnod] )
    return nodID 
end
addnode!(model::Model,coord::ℝ1)  = addnode!(model,reshape(coord,(1,length(coord))))[1]

coord(nod::AbstractVector{Node}) = [n.coord for n∈nod]

function addelement!(model::Model,::Type{T},nodID::ℤ2;kwargs...) where{T<:AbstractElement}
    # new element type? make space in model.eletyp and model.eleobj for that
    ele1     = T(collect(model.nod[nodID[begin,:]]);kwargs...)
    E        = typeof(ele1)

    # add eletyp to model (if new)
    eletypID = geteletypID(model,E)
    if eletypID == 0 # new element type!
        eletypID = maxeletypID(model)+1
        push!(model.eletyp, EleTyp(eletypID,E))
        push!(model.eleobj, Vector{E}()       )
    end        

    iele_sofar       = length(model.eleobj[eletypID])
    nele_new,nnod    = size(nodID)
    inod,class,field = getdoflist(E)
    ndof             = getndof(E)
    nnod == getnnod(E) || muscadeerror(@sprintf "Connecting element of type %s: Second dimension of inod (%i) must be equal to element's nnod (%i)" T nnod getnnod(E) ) 

    eleID  = maxeleID(model) .+ collect(1:nele_new)  # allocate return variable
    dofID  = Vector{𝕫      }(undef,ndof    )         # work array
    ele    = Vector{Element}(undef,nele_new)         # work array - will be appended to model.ele
    eleobj = Vector{E      }(undef,nele_new)         # work array - will be appended to model.eleobj
    for iele_new = 1:nele_new
         # add eleID to nodes
        for nod ∈ model.nod[unique(nodID[iele_new,:])] # unique: if several nodes of an element are connected to the same model node, mention element only once
            push!(nod.eleID,eleID[iele_new])
        end
        # for all dofs of the element, make sure they are represented in the nod object and get an ID
        for idof = 1:ndof
            nodid    = nodID[iele_new,inod[idof]]  
            # add doftyp to model (if new)  
            doftypID = getdoftypID(model,class[idof],field[idof])
            if doftypID == 0 
                doftypID = maxdoftypID(model)+1
                push!(model.doftyp, DofTyp(doftypID,class[idof],field[idof],1.,𝕫[])) # do not add dof to doftyp, though
            end
            # add dof to model (if new)
            dofID[idof]    = getdofID(model,nodid,class[idof],field[idof])
            if dofID[idof] == 0 
                dofID[idof] = maxdofID(model)+1
                push!(model.dof, Dof(dofID[idof],nodid,doftypID,𝕫[]) ) # do not add element to dof, though
            end
            # add dof to doftyp (always)
            push!(model.doftyp[doftypID].dofID, dofID[idof])
            # add element to dof (always)
            push!(model.dof[dofID[idof]].eleID, eleID[iele_new])
            # add dof to node (if new)
            if dofID[idof] ∉ model.nod[nodid].dofID
                push!(model.nod[nodid].dofID, dofID[idof])
            end
        end
        # add element to model (always)
        ele[   iele_new] = Element(eleID[iele_new], eletypID, iele_sofar+iele_new, nodID[iele_new,:],dofID)
        eleobj[iele_new] = iele_new==1 ? ele1 : T(collect(model.nod[nodID[iele_new,:]]);kwargs...)   # call element constructor
    end
    append!(model.ele,ele   )
    append!(model.eleobj[eletypID],eleobj)
    return eleID 
end
addelement!(model::Model,::Type{E},nodID::ℤ1;kwargs...) where{E<:AbstractElement} = addelement!(model,E,reshape(nodID,(1,length(nodID)));kwargs...)[1] 

# Model inspection

function show_dof(model::Model,dofID::𝕫)
    d = model.dof[dofID]
    t = model.doftyp[d.doftypID]
    @printf("Degreee of freedom %i\n   class = :%s, field = :%s, scale = %g (doftypID = %i)\n   node = %i, eleID = ",d.ID,t.class,t.field,t.scale,d.doftypID,d.nodID)
    print(d.eleID)
    @printf("\n")
end
function show_nod(model::Model,nodID::𝕫)
    n = model.nod[nodID]
    @printf("Node %i\n   coord = ",n.ID)
    print(n.coord)
    @printf("\n   dof   = ")
    print(n.dofID)
    @printf("\n   ele   = ")
    print(n.eleID)
    @printf("\n")
end
function show_ele(model::Model,eleID::𝕫)
    e = model.ele[eleID]
    @printf("Element %i   (eletypID,iel)=(%i,%i)\n   nodes = ",e.ID,e.eletypID,e.iele)
    print(e.nodID)
    @printf("\n   dof   = ")
    print(e.dofID)
    @printf("\n")
end

