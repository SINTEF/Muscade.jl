# mut,opt = allocategraphicdata(axis::Axis,o::Vector{EulerBeam3D{Tmat,Udof}};kwargs...)
# mut = updategraphicdata(axis::Axis,o::Vector{EulerBeam3D{Tmat,Udof}}, Λ,X,U,A,t,SP,dbg,oldmut,opt)
# draw!(axis::Axis,::Type{EulerBeam3D{Tmat,Udof}},mut,opt)

# Model drawing
# typestable kernel

using Observables: Observable

Graphic{Tobs,Topt,Taxis} = @NamedTuple{obs::Tobs, opt::Topt, axis::Taxis}
# Single call draw! for one element type
function draw!(axis,eleobj::AbstractVector{Eletyp}, Λ,X,U,A,t,SP,dbg;kwargs...) where{Eletyp<:AbstractElement}
    mut,opt = allocate_drawdata(axis,eleobj;kwargs...)
    mut     = update_drawdata(  axis,eleobj,mut,opt, Λ,X,U,A,t,SP,dbg)
    obs     = isnothing(mut) ? nothing : map(Observable,mut)
    draw!(                      axis,Eletyp,obs,opt)
    return obs,opt
end
function draw_!(axis,dis::EletypDisassembler,eleobj::AbstractVector{Eletyp},iele,state,dbg;kwargs...) where{Eletyp} 
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
            U[jder][:,i] = state.U[jder][index.U]
        end
        A[:,i]  = state.A[index.A]
    end
    obs,opt = draw!(axis,eleobj, Λ,X,U,A,state.time,state.SP,dbg;kwargs...)
    return (obs=obs,opt=opt,axis=axis)
end

"""
    graphic = draw!(axis   ,state[,els];kwargs...)
              draw!(graphic,state[,els];kwargs...)

Plot all or part of a `Model`.

Currently, only `GLMakie.jl` is supported and tested, but `Muscade` is designed to allow application developers to 
chose other graphic system, including exporting data to Paraview. `GLMakie.jl` is thus not a dependency of `Muscade`,
and must be installed and invoked (`using`) separately to run demos provided with `Muscade`.

Application developers can implement methods [`Muscade.allocate_drawdata`](@ref), [`Muscade.update_drawdata`](@ref) and
[`Muscade.draw!`](@ref) to make their element "drawable".

`axis` a `GLMakie.jl` `Axis`, a `Muscade.SpyAxis` (for automated testing of graphic generation), and in the future 
         a HDF5/VTK file handle for export of data to Paraview.   
`state` a single `State`.
`els` specifies which elements to draw and can be either
- a vector of `EleID`s (obtained from [`addelement!`](@ref)`), all corresponding
  to the same concrete element type
- a concrete element type (see [`eletyp`](@ref)).
- omitted: all the element of the model are drawn.
`kwargs...` is any additional key words arguments that will be passed to the `draw` method of each element, 
for example to specify colors, etc.  See the elements' documentation.

When a plot of the `Model` is first generated, `axis` must be provided, and `draw!` returns `graphic`. `graphic` can
then be provided for further calls to `draw!` to update the graphic.

See also: [`getdof`](@ref), [`@request`](@ref), [`@espy`](@ref), [`addelement!`](@ref), [`solve`](@ref)
"""
function draw!(axis,state::State,eleID::Vector{EleID};kwargs...)     # Some elements, all of same concrete type
    ieletyp              = eleID[begin].ieletyp
    all(e.ieletyp== ieletyp for e∈eleID) || muscadeerror("All elements must be of the same element type.")
    dis                  = state.dis.dis[ieletyp]
    iele                 = [e.iele for e∈eleID]
    eleobj               = view(state.model.eleobj[ieletyp],iele)
    graphic              = draw_(axis,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    return graphic
end
function draw!(axis,state::State,::Type{E};kwargs...) where{E<:AbstractElement}  # All elements of given concrete type
    ieletyp              = findfirst(E.==eletyp(state.model))
    isnothing(ieletyp)  && muscadeerror("This type of element is not in the model.'")
    eleobj               = state.model.eleobj[ieletyp]
    dis                  = state.dis.dis[ieletyp]
    iele                 = eachindex(eleobj)
    graphic              = draw_!(axis,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    return graphic
end    
function draw!(axis,state::State;kwargs...)   # whole model
    neletyp              = getneletyp(state.model)
    graphic              = Vector{Graphic}(undef,neletyp)
    for ieletyp          = 1:neletyp
        eleobj           = state.model.eleobj[ieletyp]
        dis              = state.dis.dis[ieletyp]
        iele             = eachindex(eleobj)
        graphic[ieletyp] = draw_!(axis,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    end
    return graphic
end   

# to update an existing graphic
function to_observable!(obs::Observable,v)
    obs[] = v
end
function draw_!(graphic::Graphic,dis::EletypDisassembler,eleobj::AbstractVector{Eletyp},iele,state,dbg;kwargs...) where{Eletyp} 
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
            U[jder][:,i] = state.U[jder][index.U]
        end
        A[:,i]  = state.A[index.A]
    end
    mut     = map(Observables.to_value,graphic.obs)
    mut     = update_drawdata(graphic.axis,eleobj,mut,graphic.opt, Λ,X,U,A,state.time,state.SP,(dbg...,iele=iele))
    foreach(to_observable!,graphic.obs,mut) # TODO do syntax
    draw!(                    graphic.axis,Eletyp,graphic.obs,graphic.opt)
end
function draw!(graphic::Vector{Graphic},state::State;kwargs...)   # whole model
    for ieletyp ∈ eachindex(state.model.eleobj)
        eleobj          = state.model.eleobj[ieletyp]
        dis             = state.dis.dis[ieletyp]
        iele            = eachindex(eleobj)
        draw_!(graphic[ieletyp],dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    end
end   
function draw!(graphic::Graphic,state::State,eleID::Vector{EleID};kwargs...)     # Some elements, all of same concrete type
    ieletyp             = eleID[begin].ieletyp
    all(e.ieletyp== ieletyp for e∈eleID) || muscadeerror("All elements must be of the same element type.")
    dis                 = state.dis.dis[ieletyp]
    iele                = [e.iele for e∈eleID]
    eleobj              = view(state.model.eleobj[ieletyp],iele)
    draw_!(graphic,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
end
function draw!(graphic::Graphic,state::State,::Type{E};kwargs...) where{E<:AbstractElement}  # All elements of given concrete type
    ieletyp             = findfirst(E.==eletyp(state.model))
    isnothing(ieletyp) && muscadeerror("This type of element is not in the model.'")
    eleobj              = state.model.eleobj[ieletyp]
    dis                 = state.dis.dis[ieletyp]
    iele                = eachindex(eleobj)
    draw_!(graphic,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
end    

"""
    axis = Muscade.SpyAxis

Spoof a GLMakie `axis` object so that calls like

    using Muscade: lines!
    lines!(  axis,args...;kwargs...) 
    
result in `args` and `kwargs` being stored in `axis`, allowing to test functions that generate plots.
Results are accessed by for example

    axis.call[3].fun        
    axis.call[3].args[2]

To get the name of the 3rd GLMakie function that was called, and the
2nd input argument of the call.

Only `lines!`, `scatter!` and `mesh!`  are implemented for now, but more functions can
easily be added.

In given Julia session, if GLMakie is used, then Muscade: lines! (etc.) cannot be used,
and conversedly: restart Julia when switching.  

The reason is that Muscade defines `lines!` (etc.) instead of overloading it.
This is deliberate, to allow Muscade to run unit tests of graphical generation, without making
GLMakie a dependency of Muscade.  This again is because element developers remain free
to base `draw` for their suite of elements on other graphic packages.
"""
struct SpyAxis
    call::Vector{Any}
end
SpyAxis() = SpyAxis(Any[])
lines!(  axis::SpyAxis,args...;kwargs...) = push!(axis.call,(fun=:lines!  ,args=args,kwargs=kwargs))
scatter!(axis::SpyAxis,args...;kwargs...) = push!(axis.call,(fun=:scatter!,args=args,kwargs=kwargs))
mesh!(   axis::SpyAxis,args...;kwargs...) = push!(axis.call,(fun=:mesh!   ,args=args,kwargs=kwargs))