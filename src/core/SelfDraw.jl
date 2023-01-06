export draw,request2draw

request2draw(::AbstractElement) = () # by default, an element does not need access to element-results to plot
draw(axe,key,out,::AbstractElement;kwargs...) = nothing # by default, an element draws nothing
function draw(axe,key,out,dis::EletypDisassembler,eleobj,iele,state,dbg;kwargs...) 
    # typestable kernel
    for ie ∈ iele
        index = dis.index[ie]
        Λe    = state.Λ[index.X]                 
        Xe    = Tuple(x[index.X] for x∈state.X)
        Ue    = Tuple(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        eo    = eleobj[ie]
        draw(axe,key,out, eo, Λe,Xe,Ue,Ae, state.time,state.ε,(dbg...,iele=ie);kwargs...)
    end
    return
end
function draw(axe,state::State,ieletyp::𝕫; iele::ℤ1=1:length(state.model.ele[ieletyp]),kwargs...)
    # User syntax 2: One element type, some or all elements within the types
    eleobj              = state.model.eleobj[ieletyp]
    dis                 = state.dis.dis[ieletyp]  
    key,nkey            = makekey(request2draw(eltype(eleobj)),espyable(eltype(eleobj)))
    out                 = 𝕣1(undef,nkey) # allocating
    draw(axe,key,out,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    return
end
function draw(axe,state::State;ieletyp::ℤ1=1:length(state.model.ele),kwargs...)
    # User syntax 1: Draw several element types -  cannot specify iele
    for et ∈ ieletyp
        draw(axe,state,et;kwargs...) # call user syntax 2
    end
    return
end
