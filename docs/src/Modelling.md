# Creating a model

## Script as input

`Muscade` being a framework for the development of optimization-FEM applications, it does not provide the elements needed to treat any specific application. It only provides a limited number of generic modeling capabilities, like fixing degrees of freedom (dofs) to describe boundary conditions, introducing holonomic constraints or costs on either dofs, or element-results (see [Built-in elements](@ref)).

Hence to create a model, one will typicaly be `using` both `Muscade` and another package that provides a `Muscade`-based application (app).  The app provides specific elements for domains like continuum mechanics, marine structures, hydrogen diffusion etc.

Input to such an app is provided in the form of a Julia script containing instructions (calls to `Muscade`, using elements provided by the app) to define the model, execute analyses, and extract and process analysis results.  In other words, scripting a series of analyses, or some specific pre or postprocessing is simply done in the same script, and app developpers do not have to write code pertaining to a user interface. That being said, an app could introduce a GUI that would itself do the calls to `Muscade`.

Here is a simple example of analysis:

```julia
using Muscade
using StaticArrays

model           = Model(:TestModel)
n1              = addnode!(model,[0.]) 
n2              = addnode!(model,[1.])
e1              = addelement!(model,Hold,[n1];field=:tx1)                       # Hold first node
@once id1 load(t) = 3t
e2              = addelement!(model,DofLoad,[n2];field=:tx1,value=load)        # Increase load on second node
@once id2 res(X,X′,X″,t)  = 12SVector(X[1]-X[2],X[2]-X[1])
e3              = addelement!(model,QuickFix,[n1,n2];inod=(1,2),field=(:tx1,:tx1),
                              res=res)  # Linear elastic spring with stiffness 12
initialstate    = initialize!(model)
state           = solve(SweepX{0};initialstate,time=[0.,1.],verbose=false)      # Solve the problem
tx1,_           = getdof(state[2],field=:tx1,nodID=[n2])                        # Extract the displacement of the free node
req             = @request F                                                    # Extract internal results from the spring element
eleres          = getresult(state,req,[e2]) 
iele,istep      = 1,2
force           = eleres[iele,istep].F
@show tx1
@show force
```

## Model definition

The definition of a model is done in three phases:

1. Creating a blank model, with [`Model`](@ref).
2. Adding nodes and elements, with [`addnode!`](@ref) and [`addelement!`](@ref). One can however add an element to the model *after* all the nodes of the element have been added to the model.
3. Initialising the model with [`initialize!`](@ref).  Once this is done, one can no longer add nodes or elements to the model. [`initialize!`](@ref) hashes some tables and generates an initial "as meshed" state of the system. Typicaly (but this depends on the app), all dofs are set to zero. The resulting variable, here called `initialstate` contains a pointer to the model: passing a state to a solver makes the model available to the solver. 
[`setdof!`](@ref) can be used to set the value of specific dofs for more specific initial conditions.

`Muscade` does not provide a mesher. There are some general purposes meshers with Julia API, which outputs could be used to generate calls to [`addnode!`](@ref) and [`addelement!`](@ref).

Note that two `function`s, `load` and `res` are defined in the script, and then passed as argument to element constructors. In the script it is *recommended* (but not compulsory) to annotate the function definition with the macro[`@once`](@ref).  The first argument must be a unique variable name. The second argument is the function definition.  The macro prevents the function to be re-parsed if unchanged, which in turn prevents unnecessary recompilations of Muscade when the script is runned multiple times in a session. 

The model - either finitialized or under construction, can be examined using [`describe`](@ref) and [`getndof`](@ref).  

Optionaly, one can also use [`setscale!`](@ref) (with the help of [`studyscale`](@ref)) to scale the variables and thus improve the conditioning of the problem. 

Information on commands provided by Julia and packages (including ``Muscade``) can be obtained from the help mode in the REPL.  Make sure the command is available by `using Muscade`, then activate the help mode by pressing `?`. 

```julia
help?> Hold

  Hold <: AbstractElement


  An element to set a single X-dof to zero.

  Named arguments to the constructor
  ≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡≡

    •  field::Symbol. The field of the X-dof to constraint.

    •  λfield::Symbol=Symbol(:λ,field). The field of the Lagrange multiplier.  
(...)
```

## Built-in elements

With a few exceptions for testing and demonstration, `Muscade` does not provide physical elements.  However, it provides several general purpose elements  to introduce loads, costs or  constraints.

[`DofLoad`](@ref) adds a time-varying load on a single ``X``-dof.  Elements for more general loads, in particular, consistent loads on element boundaries or domain, or follower loads, need to be implemented if required.

[`DofCost`](@ref) adds a cost as a function of either ``X``-dofs ,``U``-dofs (and/or their derivatives), ``A``-dofs and time, or as a function of ``A``-dofs alone. Elements for costs on unknwn distributed load *fields* (over boundary or domain) must be provided by apps if required.

[`SingleDofCost`](@ref) provides a simplified syntax for costs on a single dof.

[`SingleUdof`](@ref) allows to define an unknown external nodal load and apply a cost to it.

[`ElementCost`](@ref) adds a cost on a combination of one element's dofs and element-results.

[`DofConstraint`](@ref) adds a constraint to a combination of *values* (no time derivatives) of dofs. The constraints can switch over time between equality, inequality and "off". Inequality constraints are handled using a modified interior point method.

[`ElementConstraint`](@ref) adds a constraint to a function of internal results from one element. The constraints can switch over time between equality, inequality and "off". Inequality constraints are handled using a modified interior point method.

[`Hold`](@ref) provides a simplified syntax to set a single ``X``-dof to zero.

[`QuickFix`](@ref) allows to rapidly create a simple element, with limitations in functionality. 

## Running the analysis

[`solve`](@ref) is then called with the name of the solver to be used (here [`SweepX{0}`](@ref)), and any named parameters required by the solver. The return value `state` can have different structures, depending on the solver.  For [`SweepX{0}`](@ref), `state` is a vector (over the time steps) of `State`s.

[`describe`](@ref) can also be used to inspect `State`s.

Analyses may fail due to singular matrix.  The source of the singularity can be challenging to diagnose. [`studysingular`](@ref) can help determine the null-space of an incremental matrix, for small problems.

## Extracting results

`State`s (returned by [`initialize!`](@ref) and [`solve`](@ref)). are variables which contents are private (not part of the API, and subject to change), but can be accessed using [`getdof`](@ref) and [`getresult`](@ref).

 [`getdof`](@ref) allows to obtain dofs which are directly stored in `state`, by specifying class, field and node.

[`getresult`](@ref) (used in combination with [`Muscade.@request`](@ref)) allows to obtain "element-results".  Element-results are intermediate values that are computed within [`Muscade.lagrangian`](@ref) or [`Muscade.residual`](@ref), but are (generaly) not returned, because the API for these functions does not open for this.  In mechanics,  [`Muscade.residual`](@ref) would take displacements as inputs (``X``-dofs) and from them compute the forces that must act on the element to cause these displacements. Element-results woudl then include quantities such as stresses and strains.  To be requestable, a variable must be tagged in [`Muscade.lagrangian`](@ref) or [`Muscade.residual`](@ref), prefixing its name with `☼` (`\sun`) at the right hand of an assigment.

These element-results are not stored in the `State`, and tagging variables does not result in either increase storage or computing time: [`getresult`](@ref) will compute requested values on the fly by calling a modified version of [`Muscade.lagrangian`](@ref) or [`Muscade.residual`](@ref) generated by [`@espy`](@ref).  This also implies that one does not need to decide on what variables to store before an anlysis, a great advantage for heavy analyses.

## Units

`Muscade` provides functionality to transform quantities to and from basic SI units.

```julia
using Muscade, Printf
using Muscade: m, kg, pound, foot
rho          = 3←(pound/foot^3)                      # convert to SI
vieuxquintal = 1000*pound                            # define new unit
@printf("Density [pound/foot^3] %f",rho→pound/foot^3) # convert from SI
```

Arrays can be converted in the same way: `[200,300,24]←mm`.

A guideline for handling units without [problems](https://en.wikipedia.org/wiki/Mars_Climate_Orbiter) is:

- **Element developers** assume inputs with consistent units, and thus never make unit conversions.
- **Element developers** do not assume that the input are expressed in base SI units, and thus require all necessary dimensional constants (acceleration of gravity, gas constant...) as user input.
- **Users** convert all their input values as they define them in the input `rho = 3 ← pound/foot^3`.
- **Users** convert Muscade outputs just before printing them out `printf("stress [MPa] %f",stress → MPa)`.

Excellent packages exist for the handling of units ([`Unitful.jl`](https://painterqubits.github.io/Unitful.jl/stable/) ).  These packages have zero
runtime overhead, and allow to verify code for unit consistency (`Muscade` does not provide this). However, it is arguably not possible to make these packages work with `Muscade`: In `Muscade`, `3←(pound/foot^3)` is of type `Float64`.  A comparable operation in [`Unitful.jl`](https://painterqubits.github.io/Unitful.jl/stable/) would output a variable with a *type* containing data about dimensionality. `Muscade` handles various arrays of quantities with different dimensionality: such a solution would result in arrays of heterogeneous types. `Muscade` does not allow this, as this would result in catastrophic loss of performance due to [type instability](@ref typestab).

## Drawing

TODO 
