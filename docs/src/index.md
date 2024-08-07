```@meta
CurrentModule = Muscade
```

!!! info
    `Muscade.jl` is under development. New solvers will be added.

# [Purpose](@id purpose)

**[Muscade.jl](https://github.com/SINTEF/Muscade.jl): Framework for the description and solution of optimization-FEM models.**

`Muscade` allows to describe and solve optimization problems constrained by the equilibrium of a FEM model.

## Examples of FEM-optimization problems

Starting with some terminology, FEM-optimization problems have more degrees of freedom (refered to as "dofs" in the following) than classical FEM problems. In `Muscade`, dofs are separated into three `classes`:

- `X`-dofs are the classical dofs of a finite element model. `X`-dofs are duals of the FEM equations. 
- `U`-dofs typicaly represent unknown external loads acting on a structure, that vary with time.
- `A`-dofs typicaly represent unknown model parameters or design parameters that remain constant in time.

Several types of problems fall under this description of "FEM-optimisation problems":

**Design optimization**: What is the cheapest way to engineer a system whose strength is described by `A`, that will survive a set of loading conditions?  The target function describes a monetary cost.

**Reliability analysis**: Finding a design point. What is the most probable combination of external loads `U` and strength of the structure `A` that may cause the response `X` to exceed, in one of many ways, an acceptable limit? Here the target function in the optimization problem describes the probability density (or rather, its logarithm) of unknown loads and strength parameters.

**Load identification and monitoring**: Given incomplete and noisy measurements of the response of a system, what are the loads that are most likely to have caused a response close to what has been measured?  The target function describes prior knowledge of the load processes, and the type, value and precision of the measurements. 

**Model identification**: Adjust the model of a system (model calibration, damage detection) given measurements on its response when responding to at least partly unknown load. As in load identification problems, the target function describes prior knowledge of the load processes, and the measurements.  In addition, it expresses prior knowledge of the state of the structure. 

**Optimal control**: how to steer a system with many dofs into a wanted behaviour? The target function describes the cost of actuation and the cost of deviation from a target behaviour.

## Elements in Muscade

`Muscade` provides an API to create new elements. While `Muscade` comes with a few elements for testing and demonstration purpose, the idea is to support users in developing their own element formulations and implementations.  

`Muscade` inherently supports the developement of multiphysics model:  Elements can introduce new *fields* (new types of degrees of freedom: temperature, chenical species concentration, electric potential...). Elements that model the physics of a problem require the implementation a function ([`Muscade.residual`](@ref)) that, given degrees of freedom, returns the element's contribution to the residual of the discretized equations to be solved, representing for example dynamic equilibrium. 

Elements are also used to implement constraints, or, for example, measurements: given degrees of freedom, such an element returns the logarithm of the probability of this combination of degrees of freedom, given the measured data ([`Muscade.lagrangian`](@ref)).  

`Muscade` makes it easier to develop new element types by applying several techniques:

**Automatic differentiation**: elements describing a physical phenomena need only return the contribution of the element to the residual of the equations.  Elements describing a contribution to the target function do not need to provide gradient or Hessian.

**Automatic extraction of element-results**: "Element-results" are intermediate results in the process of computing contributions to the residual from degrees of freedom. An example is stresses and strains, in an element that given nodal displacements computes nodal forces. To make an element-result available, the only requirement is to prefix it with a `☼` (special character `\sun`).  For example: `☼σ = E*ε`.  This automatic result extraction also makes it possible to describe measurements taken on element-results (for example a strain).

**Element constructor for modeling**: Models in `Muscade` are created by writing Julia code that among other things, constructs new element instances.  Thanks to this, an element constructor is all that is required of the element programmer to facilitate model creation by the model user.




