# FEM-optimization

## Classes of degrees of freedom

`Muscade` introduces 3 *classes* of degrees of freedom (dofs). 

**Xdofs** are the dofs normaly encountered in normal "forward" FEM analysis.  They provide a discrete representation of the response of the system. There is a one-to-one relation (a "duality") between the Xdofs ``X`` and the *residuals* ``R``, which are the discretized form of he differential equations we seek to solve. In forward FEM, we seek to solve a problem of the form

```math
\forall t, R(X(t),t)=0
```

**Udofs** are additional dofs that can be used to represent additional loads on the system. Like Xdofs, Udofs are time-dependant. Unlike Xdofs, there str no residuald (no new equations) corresponding to them.

**Adofs** are additional dofs that can be used to parameterise the model. Adofs are *not* time-dependant,and there are no residuals corresponing to them.

Equilibrium now requires

```math
\forall t, R(X(t),U(t),A,t)=0
```

which is ill-posed (unknowns ``U`` and ``A`` have been added, but there are no new equations).

## Lagrangian

To simplify the presentations in the following, we drop time ``t`` from the notations.  This can be interpreted in two ways:

1. We are looking at a problem at a single instant (static).
2. We are solving an evolution problem (dynamic), and ``X`` and ``U`` are now vectors-valued functions of time, and ``R``represents a vector of ordinary differential equations (over time). 

To make the problem well-posed again, we introduce a cost (or target, or surprisal) function ``Q(X,U,A)``.  We then seek to make ``Q`` extremal (finding a local minimum, maximum), under the constraint that the above equations be verified.  This is equivalent to finding an extremal (a saddle point) of the *Lagrangian* (not to be confused with the potential in Lagrangian mechanics)

```math
L(\Lambda,X,U,A)= Q(X,U,A) + Λ\cdot R(X,U,A)
```

where ``Λ`` are Lagrange multipliers.  There is again a one-to-one correspondance between Λdofs and Xdofs.

For evolution problems, the above dot product includes an integral over time.  The gradients of ``L`` are ordinary differential equations (in time), found using Frechet derivatives.

## Solvers

Solvers in `Muscade` solve problems of different forms including forward `X` problems:

```math
R(X,0,0)=0
```

optimisation `XUA` problems:

```math
\begin{aligned}
\nabla_\Lambda L(\Lambda,X,U,A) &= 0 \\
\nabla_X       L(\Lambda,X,U,A) &= 0 \\
\nabla_U       L(\Lambda,X,U,A) &= 0 \\
\nabla_A       L(\Lambda,X,U,A) &= 0
\end{aligned}
```

as well as `XU` problems:

```math
\begin{aligned}
\nabla_\Lambda L(\Lambda,X,U,0) &= 0 \\
\nabla_X       L(\Lambda,X,U,0) &= 0 \\
\nabla_U       L(\Lambda,X,U,0) &= 0 
\end{aligned}
```

and `XA` problems:

```math
\begin{aligned}
\nabla_\Lambda L(\Lambda,X,0,A) &= 0 \\
\nabla_X       L(\Lambda,X,0,A) &= 0 \\
\nabla_A       L(\Lambda,X,0,A) &= 0 
\end{aligned}
```

## New elements

When creating a new element in `Muscade` to create an application, one can either implement a contribution to ``R(X,U,A)`` by implementing a method `Muscade.residual` for the element, or a contribution to ``L(\Lambda,X,U,A)`` by implementing `Muscade.lagrangian`.  For performance, when possible (when implementing an "physical" element discretizing the differential equations), use `Muscade.residual`. 