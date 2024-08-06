```@meta
CurrentModule = Muscade
```

!!! info
    `Muscade.jl` is under development. New solvers will be added.

# [Purpose](@id purpose)

**[Muscade.jl](https://github.com/SINTEF/Muscade.jl): Framework for the description and solution of optimization-FEM models.**

`Muscade` allows to describe and solve optimization problems constrained by equilibrium of a FEM model.  For example, we can have a FEM model of a structure subjected to forces that are to be estimated, based on partial and imperfect measurements of its response.  In addition, the stiffness of the structure can also be adjusted.

FEM-optimization problems have more degrees of freedom (in the following: "dofs") than classical FEM problems. They are separated into three `classes`:

- `X`-dofs: duals of the FEM equations. `X`-dofs are the classical dofs of a finite element model.
- `U`-dofs: not dual of the FEM equations and varying with time. Typicaly, `U`-dofs represent unknown external loads activing on the system.
- `A`-dofs: not dual of the FEM equations and constant over time. Typicaly, `A`-dofs represent unknown model parameters.

Several types of problems fall under this description of "FEM-optimisation problems":

**Reliability analysis**: Finding a design point. What is the most probable combination of external loads `U` and strength of the structure `A` that may cause the response `X` to exceed, in one of many ways, the acceptable? Here the target function describes the probability density (or rather, it's logarithm) describing unknown load and strength parameters.

**Design optimization**: What is the cheapest way to engineer a system (for example a structure) that will survive a set of loading conditions?  The target function describes a monetary cost.

**Load identification and monitoring**: Given incomplete and noisy measurements of the response of a system, what are the loads that are most likely to cause a response close to what has been measured?  The target function describes prior knowledge of the load processes, and the type, value and precision of the measurements. 

**Model identification**: given enough measurements on the response of a system responding to at least partly unknown load, is it possible to adjust the model of the system (model calibration, damage detection)? As in load identification problems, the target function describes prior knowledge of the load processes, and the measurements.  In addition, it expresses prior knowledge of the state of the structure. 

**Optimal control**: how to steer a system with many dofs into a wanted behaviour? The target function describes the cost of actuation and the cost of deviation from a target behaviour.

`Muscade` provides an API to create new elements. While `Muscade` comes with a few elements for testing and demonstration purpose, the idea is to support users in developing their own element formulation and implementation.  `Muscade` supports the developement of multiphysics model:  ELements can introduce new types of degrees of freedim

Elements are used to implement, for example measurements: given degrees of freedom, such an element return the logarithm of the probability of this combination of degrees of freedom, given the measured data.  Elements for this purpose are provided by `Muscade`, but new elements of this type can be added.

Elements that model the physics of a problem require a function that, given degrees of freedom, return the element's contribution to the residual of the discretized equations to be solved.  The formulation of such elements will typicaly involve the discretization of differential equation by introducing shape functions. `Muscade` is thus not a *modelling language* like for example UFL which provides a formalism in which a "well posed problem" can be defined with domain, differential equations and boundary conditions.

``Muscade` makes it easier to develop new element types by applying several techniques:

**Automatic differentiation**: elements describing a physical phenomena need only return the contribution of the element to the residual of the equations.  Elements describing a contribution to the target function do not need to provide gradient or Hessian.

**Automatic extraction of element-results**: "Element-results" are intermediate results in the process of computing contributions to the residual from degrees of freedom. An example is stresses and strains, in an element that given nodal displacements computes nodal forces. To make an element-result available, the only requirement is to prefix it with a special character `☼`.  For example: `☼σ = E*ε`.  This automatic result extraction also makes it possible to describe measurements taken on element-results (for example a strain).

**Element constructor for modeling**: Models in `Muscade` are created by writing Julia code that among other things, constructs new element instances.  Thanks to this, an element constructor is all that is required of the element programmer to facilitate model creation by the model user.




