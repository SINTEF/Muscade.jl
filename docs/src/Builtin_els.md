# Built-in elements

`Muscade` does not physical elements, implementing the solution of some differential equation: this is left to the apps built on `Muscade`.  Rather, it provide provides some general purpose elements which are practical in some setting in order to introduce loads, costs or  constraints.

[`DofLoad`](@ref) adds a time-varying load on a single X-dof.  Elements for more general loads, in particular, consistent loads on element boundaries or domain, or follower loads, need to be provided by apps if required.

[`DofCost`](@ref) adds a cost as a function of either X-dofs ,U-dofs (and/or their derivatives), A-dofs and time, or as a function of A-dofs alone. Elements for costs on unknwn distributed load *fields* (over boundary or domain) must be provided by apps if required.

[`SingleDofCost`](@ref) provides a simplified syntax for costs on a single dof.

[`ElementCost`](@ref) adds a cost on a combination of one element's dofs and element-results.

[`DofConstraint`](@ref) adds a constraint to a combination of *values* (no time derivatives) of dofs. The constraints can switch over time between equality, inequality and "off". Inequality constraints are handled using an interior point method without feasibility step.

[`ElementConstraint`](@ref) adds a constraint to a function of internal results from one element. The constraints can switch over time between equality, inequality and "off". Inequality constraints are handled using an interior point method without feasibility step.

[`Hold`](@ref) provides a simplified syntax to set a single X-dof to zero.

[`QuickFix`](@ref) allows to rapidly create a simple element. Apps should not use this.
