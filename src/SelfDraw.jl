draw(axe,::AbstractElement,args...;kwargs...) = nothing # by default, an element draws nothing
function draw_(axe,dis::EletypDisassembler,eleobj,iele,state,dbg;kwargs...) 
    # typestable kernel
    for ie ∈ iele
        index = dis.index[ie]
        Λe    = state.Λ[1][index.X]                 
        Xe    = Tuple(x[index.X] for x∈state.X)
        Ue    = Tuple(u[index.U] for u∈state.U)
        Ae    = state.A[index.A]
        eo    = eleobj[ie]
        draw(axe, eo, Λe,Xe,Ue,Ae, state.time,state.SP,(dbg...,iele=ie);kwargs...)
    end
    return
end
function draw(axe,state::State,ieletyp::𝕫; iele::ℤ1=1:length(state.model.ele[ieletyp]),kwargs...)
    # User syntax 2: One element type, some or all elements within the types
    eleobj              = state.model.eleobj[ieletyp]
    dis                 = state.dis.dis[ieletyp]  
    draw_(axe,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    return
end
function draw(axe,state::State;ieletyp::ℤ1=1:length(state.model.ele),kwargs...)
    # User syntax 1: Draw several element types -  cannot specify iele
    for et ∈ ieletyp
        draw(axe,state,et;kwargs...) # call for one element type
    end
    return
end
