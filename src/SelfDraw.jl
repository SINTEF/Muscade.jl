# Default: an element asked to draw itself does nothing
draw(::Type{<:AbstractElement}, axe,eleobj, Λ,X,U,A, t,SP,dbg;kwargs...) = nothing

# Model drawing
# typestable kernel
function draw_(axe,dis::EletypDisassembler,eleobj,iele,state,dbg;kwargs...)  
    nel      = length(iele)
    nXder    = length(state.X)
    nUder    = length(state.U) 
    nX,nU,nA = getndof(eltype(eleobj),(:X,:U,:A))
    Λ        = 𝕣2(undef,nX,nel)
    X        = ntuple(i->𝕣2(undef,nX,nel),nXder)
    U        = ntuple(i->𝕣2(undef,nU,nel),nUder)
    A        = 𝕣2(undef,nA,nel)
    for (i,ieleᵢ) ∈ enumerate(iele)
        index  = dis.index[ieleᵢ]
        Λ[:,i] = state.Λ[1][index.X]
        for jder ∈ eachindex(state.X)               
            X[jder][:,i] = state.X[jder][index.X]
        end
        for jder ∈ eachindex(state.U)             
            U[jder][:,i] = state.X[jder][index.U]
        end
        A[:,i]  = state.A[index.A]
    end
    draw(eltype(eleobj),axe, eleobj, Λ,X,U,A, state.time,state.SP,(dbg...,iele=iele);kwargs...)
    return
end

"""
    draw(state[,els];kwargs...)

`state` a single `State`.
`els` specifies which elements to draw and can be either
- a vector of `EleID`s (obtained from [`addelement!`](@ref)`), all corresponding
  to the same concrete element type
- a concrete element type (see [`eletyp`](@ref)).
- omitted: all the element of the model are drawn.
`kwargs...` is any additional key words arguments that will be passed to the `draw` method of each element, 
for example to specify colors, etc.  See the elements' documentation.

See also: [`getdof`](@ref), [`@request`](@ref), [`@espy`](@ref), [`addelement!`](@ref), [`solve`](@ref)
"""
function draw(axe,state::State,eleID::Vector{EleID};kwargs...)     # Some elements, all of same concrete type
    ieletyp             = eleID[begin].ieletyp
    all(e.ieletyp== ieletyp for e∈eleID) || muscadeerror("All elements must be of the same element type.")
    dis                 = state.dis.dis[ieletyp]
    iele                = [e.iele for e∈eleID]
    eleobj              = view(state.model.eleobj[ieletyp],iele)
    draw_(axe,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
end
function draw(axe,state::State,::Type{E};kwargs...) where{E<:AbstractElement}  # All elements of given concrete type
    ieletyp = findfirst(E.==eletyp(state.model))
    isnothing(ieletyp) && muscadeerror("This type of element is not in the model.'")
    eleobj              = state.model.eleobj[ieletyp]
    dis                 = state.dis.dis[ieletyp]
    iele                = eachindex(eleobj)
    draw_(axe,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
end    
function draw(axe,state::State;kwargs...)   # whole model
    for ieletyp ∈ eachindex(state.model.eleobj)
        eleobj              = state.model.eleobj[ieletyp]
        dis                 = state.dis.dis[ieletyp]
        iele                = eachindex(eleobj)
        draw_(axe,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    end
end   

