# Solvers

## `SweepX{ORDER}` solver

`SweepX{O}` is a non-linear solver for differential equations or order `O` in time.  `SweepX{1}` is a static solver, `SweepX{1}` an implicit-Euler solver and `SweepX{2}` a Newmark-β solver. 

!!! Info
    Muscade does not allow elements to have state variables, for example, plastic strain,
    or shear-free position for dry friction.  Where the element implements such physics, this 
    is implemented by introducing the state as a degree of freedom of the element, and solving
    for its evolution, *even in a static problem*, requires the use of `ORDER≥1`

`SweepX` solves forward FEM problems (not optimisation-FEM) (see [Theory](@ref)).  However, `SweepX` can be applied to models that have U and A-dofs. This is handled as follows: One input to `SweepX` is a `State`, which can come from [`initialize!`](@ref) or from the output of another solver. `SweepX` will keep the U- and A-dofs to the value in the input `State`. `initialize!` sets all dofs to zero, so when `SweepX` is given a `State` produced by `initialize!` the analysis starts with X-dofs equal to zero, and U- and A-dofs are kept zero throughout the analysis. 

`SweepX` handles inequality constraints (for example defined with the built-in [`DofConstraint`](@ref) element) using a simplified interior point method.

See [`SweepX`](@ref).

## `DirectXUA{OX,OU,IA}` solver

`DirectXUA` is a solver for small non-linear, static (`OX=0`) or dynamic (`OX=2`), optimisation-FEM solver, with load and model parameter identification. Given a model with costs (and possibly constraints) on U- and A-dofs, the solver will determine response (X-dofs) and unknown loads for each step (U-dofs). If (`IA=1`), the algorithm willalso estimate model parameters for the whole history (A-dofs).

Currently, the solver does not handle inequality constraints.

See [`DirectXUA`](@ref).