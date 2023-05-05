# Muscade.jl

**Create and solve optimization-FEM models.**

!!! warning
    `Muscade` is work in progress.

    - Functionality is incomplete. 
    - Documentation is incomplete.
    - APIs are subject to change.
  
## Purpose

Optimization-FEM consists in optimising a target function (here called the cost function) under the constraint that a FEM model be exactly verified. This in turn implies that the problem has more unknowns than the model has equation (or at least: more unknowns than linearly independent equations, as would be the case with “insufficient” boundary conditions).

The unknowns that are dual to the FEM equations are noted X-dofs (the *response* of the model). The rest of the unknowns can be separated into U-dofs (varying with time, generaly *unknown loads*) and A-dofs (constant over time, generaly *unknown model parameters*). The conditions for such a constrained optimization problem to be well-posed are the Babushka-Brezzi conditions, which say, in essence “if you do not restrain, then at least measure”.

## Applications of optimization-FEM

Besides solving well-posed FEM problems (obtaining the response of a system, given adequate boundary conditions and known loading terms), many applications can be, or should be, possible within `Muscade`.

**Reliability analysis**: Finding a design point. What is the most probable combination of external loads `U` and strength of the structure `A` that may cause the response `X` to exceed, in one of many ways, the acceptable?

**Design optimization**: What is the cheapest way to engineer a system (for example a structure) that will survive a set of loading conditions?

**Load identification and monitoring**: Given incomplete and noisy measurements of the response of a system, what are the loads that are most likely to cause a response close to what has been measured?

**Optimal control**: how to steer a system with many dofs into a wanted behaviour?

**Model identification**: given enough measurements on the response of a system responding to at least partly unknown load, is it possible to adjust the model of the system (model calibration, damage detection)?

**Sensor array optimization**: how best to place sensors in a system in order to support the above applications?

## Documentation

See the [documentation](https://sintef.github.io/Muscade.jl/stable/index.html).

## Installing Muscade

In the REPL, type `]` (to go into package management mode), followed by 

- `add Muscade`.
- Press the `backspace` key to leave the package manager.