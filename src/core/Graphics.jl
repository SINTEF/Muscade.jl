## NOT OPERATIONAL




draw(axe,key,out,::AbstractElement,args...) = nothing # by default, an element draws nothing

function draw(axe,key,nkey,elca,iels,Y,χo,t) # typestable kernel
    el              = elca[begin].el
    rΛ              = 𝕣1(undef,length(el.iΛdof))
    rX              = 𝕣1(undef,length(el.iXdof))
    rU              = 𝕣1(undef,length(el.iUdof))
    rA              = 𝕣1(undef,length(el.iAdof))
    χn              = rearget(χo,1)  # just allocating
    out             = 𝕣1(undef,nkey) # just allocating
    χcv             = identity
    dbg             = NamedTuple()
    for iel ∈ iels
        ca,χ = elca[iel],rearview(χo,iel)
        draw(axe,key,out, ca.el, [y[ca.iΛdof] for y∈Y],[y[ca.iXdof] for y∈Y],[y[ca.iUdof] for y∈Y],[y[ca.iAdof] for y∈Y], rΛ,rX,rU,rA ,χ,χn,χcv, t,dbg)
    end
end
function draw(axe,state::State,eltyp::DataType; iels  ::ℤ1=1:Muscade.nel(state.model,eltyp))
    # One element type, some or all elements within the types
    elca                = state.model.elca[eltyp]
    el                  = elca[begin].el
    key,nkey            = makekey(request2draw(el),requestable(el))
    y,χ,t               = state.y, state.χ, state.time
    draw(axe,key,nkey,elca,iels,y,χ[eltyp],t)
end
function draw(axe,state,eltyps::AbstractVector;args...)
    # Several element types and/or abstract element types.  cannot specify iel
    eltypmods             = collect(keys(state.model.elca)) # all the concrete element types in model
    eltyp                 = subtypeof(eltypmods,eltyps)     # all the concrete element types wanted by user (input eltyps may contain abstract types)
    for elt ∈ eltyp
        draw(axe,state,elt;args...)
    end
end
