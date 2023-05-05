# Creating a model

## Script as input

`Muscade` being a framework for the development of optimization-FEM applications, it only provides a limited number of generic modeling capabilities, like fixing degrees of freedom (dofs) to describe boundary conditions, introducing holonomic constraints or costs on either dofs, or element-results (see [Built-in elements](@ref)). `Muscade` does not provide the elements needed to treat any specific application.  Hence to create a model, one will typicaly be `using` both `Muscade` and another package that provide a `Muscade`-based application (app).  The app provides specific elements for domains like continuum mechanics, marine structures, hydrogen diffusion etc.

Input to such an app is provided in the form of a Julia script containing instructions (calls to `Muscade`, using elements and possibly solvers provided by the app) to define the model, execute analyses, and extract and process analysis results.  This has two advantages: 

1. Scripting a series of analyses, or some specific pre or postprocessing is simply done in the same script.  
2. App developpers do not need to write code pertaining to a user interface.

That said, an app could introduce a GUI that would itself do the calls to `Muscade`.

Here is a simple example of analysis:

```jldoctest; output = false
using Muscade
include("../Test/SomeElements.jl")

@once sea(t,x)  = SVector(1.,0.)*t
@once sky(t,x)  = SVector(0.,10.)
α(i)            = SVector(cos(i*2π/3),sin(i*2π/3))

model           = Model(:TestModel)
n1              = addnode!(model,𝕣[0,0,+100]) # turbine
n2              = addnode!(model,𝕣[])         # Anod for turbine 
n3              = addnode!(model,𝕣[])         # Anod for anchor
e1              =  addelement!(model,Turbine   ,[n1,n2], 
                  seadrag=1e6, sea=sea, skydrag=1e5, sky=sky)
e2              = [addelement!(model,AnchorLine,[n1,n3], 
                  Δxₘtop=vcat(5*α(i),[0.]), xₘbot=250*α(i), 
                  L=290., buoyancy=-5e3) for i∈0:2]

initialstate    = initialize!(model)
state           = solve(StaticX;initialstate,time=[0.,1.],verbose=false)

tx1  ,dofid_tx1 = getdof(state[1],         field=:tx1,nodID=[n1])
ΔL   ,dofid_ΔL  = getdof(state[1],class=:A,field=:ΔL ,nodID=[n2])
req             = @request cr,ltf
eleres          = getresult(state,req,e2) 
iele,istep      = 2,1
cr              = eleres[iele,istep].cr
ltf             = eleres[iele,istep].ltf

# output

121.62396272109176
```

## Model definition

The definition of a model is done in three phases:

1. Creating a blank model.
2. Adding nodes.
3. Adding elements.

One can actually add nodes to the model after elements have been added.  One can however only add an element to the model *after* all the nodes of the element have been added to the model.

`Muscade` does not provide a mesher. There are some general purposes meshers with Julia API that could be used.

## Running the analysis

`initialize!` is used to create an as-meshed state of the system. Typicaly (but this depends on the app), all dofs are set to zero. The resulting variable, here called `initialstate` contains a pointer to the model.  Once a model is thus initialized, one can no longer add nodes or elements to it.  This is to ensure that a model can not be modified during a sequence of analyses.

`solve` is then called with the name of the solver to be used (here `StaticX`), and any named parameters required by the solver. The return value `state` is *for this solver* a vector (over the time steps) of `State`s.

## Extracting results

`State`s (returned by `initialize!` and `solve`). are variables which contents are private (not part of the API, and subject to change), but can be accessed using functions like `getdof` and `getres`.

`getdof` allows to obtain dofs which are directly stored in `state`, by specifying class, field and node.

`getresult` (used in combination with `@request`) allows to obtain element-results which have been marked as requestable inside the function `lagrangian` or `residual` of an element. These element-results are not stored in the `State`: `getresult` will call a modified version of `lagrangian` or `residual` to obtain the `@request`ed results.

TODO Plotting results
