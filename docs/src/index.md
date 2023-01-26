```@meta
CurrentModule = Muscade
```

# Muscade.jl

## [Purpose](@id purpose)
[Muscade.jl](https://github.com/SINTEF/Muscade.jl) is a framework to support the rapid development of domain-specific, multiphysics optimization-FEM models.

"Framework" refers tot he fact that Muscade does not provide finite element types modeling any specific domain of physics: rather, it facilitates their development.  Using automatic differentiation and automatic result extraction, the code needed to implement new elements is minimized.  The API to create new elements is quite classic (given degrees of freedom, compute their duals, aka. residuals): Muscade is not a *modeling language* like for example UFL which provides a formalism in which a "well posed problem" can be defined with domain, differential equations and boundary conditions.

Degrees of freedom (dofs) in Muscade have a `field`: application developers can introduce new fields like (`:tx` for translation in `x` direction, `c` for solute concentration, etc...), and associate a scale to them, to aid solvers handle unknowns of differents orders of magnitudes.  This facilitates the implementation of multiphysical models.

Muscade handles optimization-FEM problems, that is, optimization problems constrained by equilibrium of the FEM model.
```math
\forall t, \overline{R}\left(\overline{X},t) = \overline{0}    
```

## [Reference manual](@id reference)
```@index
```

```@autodocs
Modules = [Muscade]
Order   = [:module,:constant,:type,:function,:macro, :type]
```
