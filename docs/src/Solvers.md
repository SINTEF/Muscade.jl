# [Solvers](@id solvers)

## `EigX`: standard FEM modal analysis solver

The solver computes the eignmodes and oscillation frequencies of a model.

See the reference manual [`EigX`](@ref). 

## `SweepX`: standard FEM solver

`SweepX{O}` is a non-linear solver for differential equations of order `O` in time. This can be used for static and quasi static problems without hysterertic behaviour (plasticity, friction).

`SweepX{1}` is an implicit-Euler solver for differential equations of order `1` in time. This *must* be used for viscous problems, and for static and quasi static problem with hysteretic behaviour. The reason for this is that Muscade does not allow elements to have element-internal "state" variables (plastic strain, shear-free position for dry friction). Hence, where elements implement such physics, this is done by introducing the "state" as a degree of freedom of the element, and a corresponding equation.  This equation is the equation of evolution of the "state" variable, which involves the first order derivative of the variable in question *even in a static problem*.

`SweepX{2}` is a Newmark-β solver for differential equations of order `2` in time. A typical application is in structural dynamics. 

`SweepX` solves forward FEM problems (not optimisation-FEM) (see [Theory](@ref)).  However, `SweepX` can be applied to models that have ``U``- and ``A``-dofs. This is handled as follows: One input to `SweepX` is a `State`, which can come from [`initialize!`](@ref) or from the output of another solver. `SweepX` will keep the ``U``- and ``A``-dofs to the value in the input `State`. `initialize!` sets all dofs to zero, so when `SweepX` is given a `State` produced by `initialize!` the analysis starts with ``X``-dofs equal to zero, and ``U``- and ``A``-dofs are kept zero throughout the analysis. 

`SweepX` handles inequality constraints (for example defined with the built-in [`DofConstraint`](@ref) element) using a simplified interior point method.

See the reference manual [`SweepX`](@ref).   

## `DirectXUA`: non-linear inverse solver

`DirectXUA` is a solver for non-linear, static (`OX=0`), first order (`OX=1`) or dynamic (`OX=2`), optimisation-FEM problems. The same remarks on "state" variables and the choice of `OX` as for `SweepX` apply here. 

`DirectXUA` is designed for load and model parameter identification. Given a model with costs (and possibly constraints) on U- and ``A``-dofs, the solver will determine response (``X``-dofs) and unknown loads for each step (``U``-dofs). If (`IA=1`), the algorithm will also estimate model parameters for the whole history (``A``-dofs).

`OU` specifies the order of the derivatives of ``U``-dofs that appear in the target function.  For example, if `OU=0`, then costs should only be associated to the *value* of unknown external loads: the prior information on the unknown load process is that it is a white noise process. `OU≥1` allows to provide prior information in the form of colored processes.  

!!! info
    The solver is "direct" in that it solves all the degrees of freedom *at all the steps* at the same time. This introduces a limitation on the number of degrees of freedom and time steps that can be handled.  Improving performance for large problems is on-going work.

!!! info
    Currently, the solver does not handle inequality constraints.

See the reference manual [`DirectXUA`](@ref).

## `FreqXU`: linear inverse solver

`FreqXU` is a solver for *linear*, static (`OX=0`), first order (`OX=1`) or dynamic (`OX=2`), optimisation-FEM problems. The same remarks on "state" variables and the choice of `OX` as for `SweepX` apply here. The same remark on the choice of `OU` as for `DirectXUA` applies here.

"Linear" is to be understood as follows: the solver computes the Hessian (second order derivative) of the Lagrangian at a reference time, and assumes that this Hessian does not vary over time.  This allows to use the Fourier transform to transform the set of differential equations into algegraic equations int he frequency domain: the amplitudes of the ``X``- and ``U``-dofs can be solved for each frequency separately, so that the algorithm is linear in the number of steps.

The Hessian of the Lagrangian include the tangential matrices of the finite element model.  Where costs are interpreted as suprisals, this amounts to approximating probability distributions by a Gauss distribution.

See the reference manual [`FreqXU`](@ref).
