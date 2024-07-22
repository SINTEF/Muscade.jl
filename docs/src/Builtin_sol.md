# Built-in solvers

## `SweepX{O}` solver

`SweepX{O}` is a non-linear solver for differential equations or order `O` in time.  `SweepX{1}` is a static solver, `SweepX{1}` an implicit-Euler solver and `SweepX{2}` a Newmark-β solver. 

IMPORTANT NOTE: Muscade does not allow elements to have state variables, for example, plastic strain,
or shear-free position for dry friction.  Where the element implements such physics, this 
is implemented by introducing the state as a degree of freedom of the element, and solving
for its evolution, *even in a static problem*, requires the use of `ORDER≥1`

`SweepX` solves "normal" FEM (not optimisation-FEM) (see [Theory](@ref)).  However, `SweepX` can be applied to models that have U and A-dofs. This is handled as follows: One input to `SweepX` is a `State`, which can come from [`initialize!`](@ref) or from the output of another solver. `SweepX` will keep the U- and A-dofs to the value in the input `State`. `initialize!` sets all dofs to zero, so when `SweepX` is given a `State` produced by `initialize!` the analysis starts with X-dofs equal to zero, and U- and A-dofs are kept zero throughout the analysis. 

`SweepX` handles inequality constraints using a simplified interior point method. This works (reasonnably) well in concert with the built-in [`DofConstraint`](@ref) element.

See [`SweepX`](@ref).

## `StaticXUA` solver

!!! warning

    At this time, the handling of element-memory (see [Creating an element](@ref)) is not implemented.

`StaticXUA` is a non-linear, static, optimisation-FEM solver, with load and model parameter identification. Given a vector of static equilibrium configurations, obtained for example using `SweepX{0}`, on a model with costs (and possibly constraints) on U- and A-dofs, the solver will determine response (X-dofs) and unknown loads for each step (Udofs) as well as model parameters for the whole history (Adofs).

The solvers handles inequality constraints in the same way as `SweepX{0}`.

See [`StaticXUA`](@ref).