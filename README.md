# Muscade.jl

**Create and solve optimisation-FEM models.**
  
## Purpose

Solving FEM-constrained optimisation problems consists in optimising an objective function (here called the cost function) under the constraint that a FEM model be exactly verified. This in turn implies that the problem has more unknowns than the model has equation (or at least: more unknowns than linearly independent equations, as would be the case with “insufficient” boundary conditions).

The unknowns that are dual to the FEM equations are noted X-dofs (the *response* of the model). The rest of the unknowns can be separated into U-dofs (varying with time, generaly *unknown loads*) and A-dofs (constant over time, generaly *unknown model parameters*). The conditions for such a constrained optimisation problem to be well-posed are the Babushka-Brezzi conditions, which say, in essence “if you do not restrain, then at least measure”.

## Applications of optimisation-FEM

Besides solving well-posed FEM problems (obtaining the response of a system, given adequate boundary conditions and known loading terms), many applications can be, or should be, possible within `Muscade`.

**Reliability analysis**: Finding a design point. What is the most probable combination of external loads `U` and strength of the structure `A` that may cause the response `X` to exceed, in one of many ways, the acceptable?

**Design optimisation**: What is the cheapest way to engineer a system (for example a structure) that will survive a set of loading conditions?

**Load identification and monitoring**: Given incomplete and noisy measurements of the response of a system, what are the loads that are most likely to cause a response close to what has been measured?

**Optimal control**: how to steer a system with many dofs into a wanted behaviour?

**Model identification**: given enough measurements on the response of a system responding to at least partly unknown load, is it possible to adjust the model of the system (model calibration, damage detection)?

**Sensor array optimisation**: how best to place sensors in a system in order to support the above applications?

## Official documentation

See the [documentation](https://sintef.github.io/Muscade.jl/stable/index.html).

## Installing the latest release of Muscade.jl

In the REPL, type `]` (to go into package management mode), followed by 

- `add Muscade`.
- Press the `backspace` key to leave the package manager.

## Using the newest features of Muscade.jl

Some newer features, which have the vocation to be part of the next release, are generally pushed to the [dev branch](https://github.com/SINTEF/Muscade.jl/tree/dev). Pull that branch from the project repository, and in package mode, type `dev /path/to/local/copy/of/Muscade.jl`. The dev branch also has its online [documentation](https://sintef.github.io/Muscade.jl/dev/).

# Citation

If you use this package in published work, please cite (review on-going):

> P. Maincon and T. Sauder, "Muscade.jl : a Julia package for 
> FEM-constrained optimisation", ssrn:7268625
> (2026). https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7268625

```bibtex
@article{MainconSauder2026,
  author  = {Maincon, P. and Sauder, T.},
  title   = {Muscade.jl : a Julia package for FEM-constrained optimisation},
  journal = {SoftareX},
  year    = {2026},
  doi     = {10.2139/ssrn.7268625},
  url     = {https://papers.ssrn.com/sol3/papers.cfm?abstract_id=7268625},
}
```

Machine-readable metadata is in [`CITATION.cff`](CITATION.cff).