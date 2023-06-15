```@meta
CurrentModule = Muscade
```

!!! info
    `Muscade.jl` is still under development, and not ready for general usage:
    - New versions may introduce breaking changes in user and element APIs.
    - Graphics generation will be revised.
    - New solvers will be added.

# [Purpose](@id purpose)

**[Muscade.jl](https://github.com/SINTEF/Muscade.jl): Create and solve optimization-FEM models.**

**Framework**: `Muscade` provides an API to create new elements and new solvers.  It also provides the "piping" - facilities to create models, assemble matrices and vectors and extract and export results. As a framework, `Muscade` does not provide finite element types modeling any specific domain of physics.  It does however provide some domain-agnostic element types for doing things like introducing boundary conditions and introducing costs. 

The API to create new elements is quite classic (given degrees of freedom, compute their duals, aka. residuals): as such Muscade is not a *modelling language* like for example UFL which provides a formalism in which a "well posed problem" can be defined with domain, differential equations and boundary conditions.

In order to obtain high performance, elements have to be implemented using a specific programming style: functional programming, using immutables.

**Rapid development**: `Muscade` lets application developer focus on the interesting parts: formulating new element types and creating new solvers.  `Muscade` takes care of the tedious parts.

It further accelerates the rapid development of new element types by using automatic differentiation and automatic extraction of element-results (think of stresses and strains, in an element that given nodal displacements computes nodal forces).  This results in shorter and more readable element code.

The idea is thus to be able to produce specialised FEM applications, at low cost, to create highly efficient solution to niche problems.  It is also hope that `Muscade` will help students of the finite element method to put theory into practice.

**Multiphysics**: When creating a new element, one must define a function `doflist` to describe the `class` (see below) `field` and `node` of each degree of freedom (dof). Element developers can introduce new `field`s (node rotation, hydrogen concentration, temperature, components of the magnetic field...).  The user can then associate a `scale` to each field stating that "displacements are the order of `10^-3 m`" or "potentials are of the order of `kV`.  This is used to improve the conditioning of the problem.

**Optimization**: Muscade handles optimization-FEM problems, that is, optimization problems constrained by equilibrium of the FEM model.

Optimization-FEM consists in optimising a target function (here called the cost function) under the constraint that a FEM model be exactly verified. This in turn implies that the problem has more unknowns than the model has equation (or at least: more unknowns than linearly independent equations, as would be the case with “insufficient” boundary conditions).

The unknowns that are dual to the FEM equations are noted X-dofs, and sometimes refered to as “response dofs”. The rest of the unknowns can be separated into U-dofs (varying with time, generaly unknown loads) and A-dofs (comstant over time, generaly unknown model parameters). The conditions for such a constrained optimization problem to be well-posed are the Babushka-Brezzi conditions, which say, in essence “if you do not restrain, then at least measure”.


## Applications of optimization-FEM

Besides solving well-posed FEM problems (obtaining the repsonse of a system, given adequate boundary conditions and known loading terms), many applications can be, or should be possible within `Muscade`.

**Reliability analysis**: Finding a design point. What is the most probable combination of external loads `U` and strength of the structure `A` that may cause the response `X` to exceed, in one of many ways, the acceptable?

**Design optimization**: What is the cheapest way to engineer a system (for example a structure) that will survive a set of loading conditions?

**Load identification and monitoring**: Given incomplete and noisy measurements of the response of a system, what are the loads that are most likely to cause a response close to what has been measured?

**Optimal control**: how to steer a system with many dofs into a wanted behaviour?

**Model identification**: given enough measurements on the response of a system responding to at least partly unknown load, is it possible to adjust the model of the system (model calibration, damage detection)?

**Sensor array optimization**: how best to place sensors in a system in order to support the above applications?

