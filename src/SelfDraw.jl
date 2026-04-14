# mut,opt = allocategraphicdata(axis::Axis,o::Vector{EulerBeam3D{Tmat,Udof}};kwargs...)
# mut = updategraphicdata(axis::Axis,o::Vector{EulerBeam3D{Tmat,Udof}}, Λ,X,U,A,t,SP,dbg,oldmut,opt)
# draw!(axis::Axis,::Type{EulerBeam3D{Tmat,Udof}},mut,opt)

# Model drawing
# typestable kernel

# create, update and read observables that are the leafs of a NamedTuple-of-NamedTuples tree
recursiveObservable(mut            ) = Observable(mut)
recursiveObservable(mut::NamedTuple) = map(recursiveObservable,mut)
readObservable(obs::Observable) = obs[]
readObservable(obs::NamedTuple) = map(readObservable,obs)
function writeObservable!(obs::Observable,mut) 
    obs[] = mut
end
function writeObservable!(obs::NamedTuple,mut::NamedTuple) 
    foreach(obs,mut) do obsᵢ,mutᵢ
        writeObservable!(obsᵢ,mutᵢ)
    end
end


Graphic{Tobs,Topt,Taxis} = @NamedTuple{obs::Tobs, opt::Topt, axis::Taxis}
# Single call draw! for one element type
Base.keys(::Nothing)= :nothing
function draw_element!(axis,eleobj::AbstractVector{Eletyp}, Λ,X,U,A,t,SP,dbg;kwargs...) where{Eletyp<:AbstractElement}
    mut,opt = allocate_drawing(axis,eleobj;kwargs...)
    mut     = update_drawing(  axis,eleobj,mut,opt, Λ,X,U,A,t,SP,dbg)
    obs     = isnothing(mut) ? nothing : recursiveObservable(mut)
    display_drawing!(                      axis,Eletyp,obs,opt)
    return obs,opt
end
# Typestable draw! for all elements of given type
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
    obs,opt = draw_element!(axis,eleobj, Λ,X,U,A,state.time,state.SP,dbg;kwargs...)
    return (obs=obs,opt=opt,axis=axis)
end

"""
    graphic = draw!(axis   ,state[,els];kwargs...)
              draw!(graphic,state[,els];kwargs...)

Plot all or part of a `Model`.

Currently, only `GLMakie.jl` is supported and tested, but `Muscade` is designed to allow application developers to 
chose other graphic system, including exporting data to Paraview. `GLMakie.jl` is thus not a dependency of `Muscade`,
and must be installed and invoked (`using`) separately to run demos provided with `Muscade`.

Application developers can implement methods [`Muscade.allocate_drawing`](@ref), [`Muscade.update_drawing`](@ref) and
[`Muscade.display_drawing!`](@ref) to make their element "drawable".

`axis` a `GLMakie.jl` `Axis`, a `Muscade.SpyAxis` (for automated testing of graphic generation), and in the future 
         a HDF5/VTK file handle for export of data to Paraview.   
`state` a single `State`.
`els` specifies which elements to draw and can be either
- a vector of `EleID`s (obtained from [`addelement!`](@ref)`), all corresponding
  to the same concrete element type
- a concrete element type (see [`describe`](@ref)).
- omitted: all the element of the model are drawn.
`kwargs...` is any additional key words arguments that will be passed to the `draw` method of each element, 
for example to specify colors, etc.  See the elements' documentation.

When a plot of the [`Model`](@ref) is first generated, `axis` must be provided, and `draw!` returns `graphic`. `graphic` can
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
draw_!(graphic::Graphic{Nothing},dis::EletypDisassembler,eleobj::AbstractVector{Eletyp},iele,state,dbg;kwargs...) where{Eletyp} = nothing
function draw_!(graphic::Graphic,dis::EletypDisassembler,eleobj::AbstractVector{Eletyp},iele,state,dbg;kwargs...) where{Eletyp} 
    graphic.obs
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
    mut     = readObservable(graphic.obs)
    mut     = update_drawing(graphic.axis,eleobj,mut,graphic.opt, Λ,X,U,A,state.time,state.SP,(dbg...,iele=iele))
    writeObservable!(graphic.obs,mut)
    # foreach(graphic.obs,mut) do obsᵢ,mutᵢ
    #     obsᵢ[] = mutᵢ
    # end
    display_drawing!(graphic.axis,Eletyp,graphic.obs,graphic.opt)
    return graphic
end
function draw!(graphic::Vector{Graphic},state::State;kwargs...)   # whole model
    for ieletyp ∈ eachindex(state.model.eleobj)
        eleobj          = state.model.eleobj[ieletyp]
        dis             = state.dis.dis[ieletyp]
        iele            = eachindex(eleobj)
        draw_!(graphic[ieletyp],dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    end
    return graphic
end   
function draw!(graphic::Graphic,state::State,eleID::Vector{EleID};kwargs...)     # Some elements, all of same concrete type
    ieletyp             = eleID[begin].ieletyp
    all(e.ieletyp== ieletyp for e∈eleID) || muscadeerror("All elements must be of the same element type.")
    dis                 = state.dis.dis[ieletyp]
    iele                = [e.iele for e∈eleID]
    eleobj              = view(state.model.eleobj[ieletyp],iele)
    draw_!(graphic,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    return graphic
end
function draw!(graphic::Graphic,state::State,::Type{E};kwargs...) where{E<:AbstractElement}  # All elements of given concrete type
    ieletyp             = findfirst(E.==eletyp(state.model))
    isnothing(ieletyp) && muscadeerror("This type of element is not in the model.'")
    eleobj              = state.model.eleobj[ieletyp]
    dis                 = state.dis.dis[ieletyp]
    iele                = eachindex(eleobj)
    draw_!(graphic,dis,eleobj,iele,state,(ieletyp=ieletyp,);kwargs...) # call kernel
    return graphic
end    


"""
    Muscade.GUI(state,refstate=state[1];dim=3,kwargs...)

Taking `state`, a `Vector` of `State`s output by various solvers, provide
a GUI to explore the results.

This assumes that elements' drawing methods are writen for [`GLMakie`](https://docs.makie.org/).

The GUI allows to intereactively amplify responses (`Λ`,`X`,`U` and `A`-dofs) to 
make then easier to visualise. For `X`-dofs, it is the difference from the `refstate` 
(by default: `state[1]`) that is amplified.

Optional keyword arguements are
- `dim`, 2 or 3 depending on whether elements assume `Axis` or `Axis3` 
- `kwargs` keywords argument, that will be passed to [`draw!`](@ref)

See also [`EigXU`](@ref)
"""
function GUI(state::AbstractVector{S},initialstate=state[1];dim=3,kwargs...) where{S<:State{nΛder,nXder,nUder}} where{nΛder,nXder,nUder}
    ## Organize the window
    fig             = Figure(size = (1500,900))
    GLMakie.activate!(title = "Muscade.jl")
    display(fig) 
    panelSteps      = fig[1,1]        
    panelEmpty      = panelSteps[1,1] 
    panelSlide      = panelSteps[2,1] 
    panelModel      = fig[1,2:3]        
    Box(panelModel, cornerradius = 20,z=1., color = :transparent)
    axisModel = if dim==3
         Axis3(panelModel,aspect=:data,viewmode=:free,perspectiveness=.5,clip=false)
    else
         Axis2(panelModel)
    end
    ## sliders
    time            = [stateᵢ.time for stateᵢ∈state]
    sg              = SliderGrid(panelSlide,
                    (label="t"      , range=time     , startvalue=time[1], snap=true, update_while_dragging=true, format = "{:.1f} s" ),
                    (label="Λ scale", range=-5:0.01:5, startvalue=0      , snap=true, update_while_dragging=true, format = "10^{:.1f}"),
                    (label="X scale", range=-5:0.01:5, startvalue=0      , snap=true, update_while_dragging=true, format = "10^{:.1f}"),
                    (label="U scale", range=-5:0.01:5, startvalue=0      , snap=true, update_while_dragging=true, format = "10^{:.1f}"),
                    (label="A scale", range=-5:0.01:5, startvalue=0      , snap=true, update_while_dragging=true, format = "10^{:.1f}"))
    obs  = (t      = sg.sliders[1].value,
            Λscale = sg.sliders[2].value,
            Xscale = sg.sliders[3].value,
            Uscale = sg.sliders[4].value,
            Ascale = sg.sliders[5].value)
    ## Model
    graphic        = draw!(axisModel,initialstate;kwargs...)  # Create graphic objects
    _ = map(obs.t,obs.Λscale,obs.Xscale,obs.Uscale,obs.Ascale) do t,Λscale,Xscale,Uscale,Ascale                                    # Then observe the sliders to update the graphic objects
        istep      = argmin(abs.(t.-time))
        ampedstate = State{nΛder,nXder,nUder}(copy(initialstate)) 
        for ider∈nΛder; ampedstate.Λ[ider] .+=  state[istep].Λ[ider]                    *exp10(Λscale) end
        for ider∈nXder; ampedstate.X[ider] .+= (state[istep].X[ider]-initialstate.X[1]) *exp10(Xscale) end
        for ider∈nUder; ampedstate.U[ider] .+=  state[istep].U[ider]                    *exp10(Uscale) end
                        ampedstate.A       .+=  state[istep].A                          *exp10(Ascale) 
        draw!(graphic,ampedstate;kwargs...)
    end
end

