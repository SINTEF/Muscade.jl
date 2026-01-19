# Diagnostic

The recommended approach is to always do things right the first time.  For those of us that have not yet adopted this wise strategy, `Muscade.jl` provides an ample supply of tools for diagnostic and testing.

## Creating models

[`Muscade.describe`](@ref) allows to study various aspects of the model being constructed.  It can also be used to display a `State` (for models of small size).

[`Muscade.study_scale`](@ref) looks at the order of magnitudes of the gradient and Jacobian of a model, by doftype. This can help to chose an adequate scaleing.

[`Muscade.study_singular`](@ref) looks for vectors than span the nullspace (aka. kernel) of the Hessian of the model.  Inother words, it looks for combinations of dofs which are undefined because there is neither a constraint or a cost associated to it.

[`Muscade.Monitor`](@ref) and element that can be "inserted between another element and the model, to monitor traffic to and from `residual` or `lagrangian`.

## Testing elements

When developing a new element, it is advisable to test the constructor, and `residual` or `lagrangian` in a direct call (outside of any Muscade solver), and examine the returned outputs.

[`Muscade.diffed_residual`](@ref) and [`Muscade.diffed_lagrangian`](@ref) can be used to test the element's automatic differentiation.
Generaly, automatic differentiation is unproblematic, but when advanced tools are used (e.g. [`revariate`](@ref) and [`chainrule`](@ref)), then the derivatives should be inspected.  

[`Muscade.@typeof`](@ref) allows to determine the type of the return variable of a function, as inferred by the compiler.  Usefull to study type stability.

[`Muscade.print_element_array`](@ref) allows to display vectors of matrices, annotating them with the descriptionsof corresponding dofs.

[`Muscade.SpyAxis`](@ref) allows to test calls to graphic generating functions in the element's [`Muscade.display_drawing!`](@ref) method.

## Testing solvers

When developing a new solver, matrix sparsity structures are of interest.

[`Muscade.plot_matrix_sparsity`](@ref) allows to plot the sparsity structure of a `SparseMatrixCSC`.

[`Muscade.plot_block_matrix_sparsity`](@ref) allows to plot the sparsity structure of a "matrix of matrices".

