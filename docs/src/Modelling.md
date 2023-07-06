# Creating a model

## Script as input

`Muscade` being a framework for the development of optimization-FEM applications, it only provides a limited number of generic modeling capabilities, like fixing degrees of freedom (dofs) to describe boundary conditions, introducing holonomic constraints or costs on either dofs, or element-results (see [Built-in elements](@ref)). `Muscade` does not provide the elements needed to treat any specific application.  Hence to create a model, one will typicaly be `using` both `Muscade` and another package that provide a `Muscade`-based application (app).  The app provides specific elements for domains like continuum mechanics, marine structures, hydrogen diffusion etc.

Input to such an app is provided in the form of a Julia script containing instructions (calls to `Muscade`, using elements and possibly solvers provided by the app) to define the model, execute analyses, and extract and process analysis results.  This has two advantages: 

1. Scripting a series of analyses, or some specific pre or postprocessing is simply done in the same script.  
2. App developpers do not need to write code pertaining to a user interface.

That said, an app could introduce a GUI that would itself do the calls to `Muscade`.

Here is a simple example of analysis:

```julia
using Muscade
using StaticArrays

model           = Model(:TestModel)
n1              = addnode!(model,[0.]) 
n2              = addnode!(model,[1.])
e1              = addelement!(model,Hold,[n1];field=:tx1)
e2              = addelement!(model,DofLoad,[n2];field=:tx1,value=t->3t)
e3              = addelement!(model,QuickFix,[n1,n2];inod=(1,2),field=(:tx1,:tx1),
                              res=(X,X′,X″,t)->12SVector(X[1]-X[2],X[2]-X[1]))

initialstate    = initialize!(model)
state           = solve(StaticX;initialstate,time=[0.,1.],verbose=false)

tx1,_           = getdof(state[2],field=:tx1,nodID=[n2])
req             = @request F
eleres          = getresult(state,req,[e2]) 
iele,istep      = 1,2
force           = eleres[iele,istep].F
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

`State`s (returned by `initialize!` and `solve`). are variables which contents are private (not part of the API, and subject to change), but can be accessed using functions like `getdof` and `getresult`.

`getdof` allows to obtain dofs which are directly stored in `state`, by specifying class, field and node.

`getresult` (used in combination with `@request`) allows to obtain element-results which have been marked as requestable inside the function `lagrangian` or `residual` of an element. These element-results are not stored in the `State`: `getresult` will call a modified version of `lagrangian` or `residual` to obtain the `@request`ed results.

TODO Plotting results
