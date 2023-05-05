# Built-in solvers

## `StaticX` solver

!!! warning

    At this time, the handling of element-memory (see [Creating an element](@ref)) is not implemented.

`StaticX` is a non-linear, static, explicit solver for "normal" FEM (not optimisation-FEM):  At a succession of times ``t``, it will solve ``R(X(t),t)=0`` (see [Theory](@ref)).

`StaticX` can be applied to models that have U and A-dofs. This is handled as follows: One input to `StaticX` is a `State`, which can come from [`initialize!`](@ref) or from the output of another solver. `StaticX` will keep the U- and A-dofs to the value in the input `State`. `initialize!` sets all dofs to zero, so when `StaticX` is given a `State` produced by `initialize!` the analysis starts with X-dofs equal to zero, and U- and A-dofs are kept zero throughout the analysis. 

`StaticX` handles inequality constraints using a simplified interior point method, without a feasibility step. This works (reasonnably) well in concert with the built-in [`DofConstraint`](@ref) element.

See [`StaticX`](@ref).

## `StaticXUA` solver

!!! warning

    At this time, the handling of element-memory (see [Creating an element](@ref)) is not implemented.

`StaticXUA` is a non-linear, static, optimisation-FEM solver, with load and model parameter identification. Given a vector of static equilibrium configurations, obtained for example using `StaticX`, on a model with costs (and possibly constraints) on U- and A-dofs, the solver will determine response (X-dofs) and unknown loads for each step (Udofs) as well as model parameters for the whole history (Adofs).

The solvers handles inequality constraints in the same way as `StaticX`.

See [`StaticXUA`](@ref).