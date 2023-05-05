# Theory

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

To simplify the presentations in the following, we drop time ``t`` from the notations.  This can be interpreted in two ways:

1. We are looking at a problem at a single instant (static).
2. We are solving an evolution problem (dynamic), and ``X`` and ``U`` are now vectors-valued functions of time, and ``R``represents a vector of ordinary differential equations (over time). 

## Target function

To make the problem well-posed again, we introduce a target function ``Q(X,U,A)``.  We then seek to make ``Q`` extremal (finding a local minimum, maximum), under the constraint ``R(X,U,A)=0``. Depending on the relevant application, ``Q`` can represent different concepts.

In design optimisation, one can associate a monetary value to ``A`` (building stronger costs more money), and to ``X`` (some system responses, including failure, would cost money). The objective is to find the design with the lowest associated cost.

In an load estimation problem where part of the response is measured, we wish to find the most probable unknown load ``U`` (large loads are not likely) and response (a computed response ``X`` that drasticaly disagrees with the actual measurements is not likely) - under the constraint that load and response verify equilibrium.  In this case, we minimize the *surprisal* ``s = -\log(P)`` instead of maximizing ``P``.  Making probabilities or surprisal extremal is equivalent.   However, since probabilities multiply (the joint probability density of independant variables is equal to the product of the marginal probability densities), the corresponding surprisals add, which fits into the system of having elements *add* contribution to the system of equations to be solved. Further, it is typicaly easier to find the e surprisals are typicaly the former is numericaly "kinder".

## Lagrangian

Making ``Q`` extremal under the constraint ``R(X,U,A)=0`` is equivalent to finding an extremal (a saddle point) of the *Lagrangian* ``L(\Lambda,X,U,A)`` (not to be confused with the potential in Lagrangian mechanics)

```math
L(\Lambda,X,U,A) = Q(X,U,A) + Λ\cdot R(X,U,A)
```

where ``Λ`` are Lagrange multipliers.  There is again a one-to-one correspondance between Λdofs and residuals, and this between Λdofs and Xdofs.  One result of this correspondance is that when implmenting a new element, the method `doflist` does not list the Λdofs (this would otherwise just have been a compulsory repetition of the list of Xdofs). 

For evolution problems, the dot product ``Λ\cdot R(X,U,A)`` includes an integral over time: the Lagrangian is a *functional*, and the gradients of ``L`` are ordinary (or integral) differential equations (in time), found using Frechet derivatives.

## Elements

In the original finite element analysis, the *elements* in a model form a partition of a domain over which differential equations are to be solved. We will call these "physical elements".  These elements provide additive contributions to the system of equations ``R(X,U,A)=0``. Further contributions can come from external loads.

When doing optimization-FEM in `Muscade`, the elements provide additive contributions to the Lagrangian scalar ``L``, instead of to the residual vector ``R``. Actualy when creating a new element in `Muscade` to create an application, one can either implement a contribution to ``R`` by implementing a method `Muscade.residual` for the element, or a contribution to ``L`` by implementing `Muscade.lagrangian`.  For performance, whenever possible (for example when implementing a physical element), prefer `Muscade.residual`. 

In `Muscade`, the broad view is taken that anything that contributes to ``L`` is an element (even if no part of the domain is actualy associated to the element).  This more general definition of "elements" includes a variety of types:

- Physical element, disretizing differential equations over a part of the domain
- Known external loads on the boundary (non-essential boundary conditions)
- Known external loads in the domain
- Constrained dofs (essential boundary conditions)
- Holonomic equality and inequality constraints (contact)
- Optimisation constraints (e.g. stresses shal not exceed some limit at any point within part of the domain)
- Response measurements (surprisal on Xdofs)
- Unknown external loads (surprisal on Udofs)
- Damage (surprisal on Adofs)
- Cost of unfavorable response (cost on Xdofs)
- Cost of actuators (cost on Udofs)
- Cost of building a system (cost on Adofs)

Because "everything" is an element, app developers can express a wide range of ideas through `Muscade`'s element API.

## Solvers

Solvers in `Muscade` solve problems of different forms including forward `X` problems

```math
R(X,0,0)=0
```

optimisation `XUA` problems

```math
\begin{aligned}
\nabla_\Lambda L(\Lambda,X,U,A) &= 0 \\
\nabla_X       L(\Lambda,X,U,A) &= 0 \\
\nabla_U       L(\Lambda,X,U,A) &= 0 \\
\nabla_A       L(\Lambda,X,U,A) &= 0
\end{aligned}
```

as well as corresponding `XU` and `XA` problems.


