# Theory

## Classes of degrees of freedom

`Muscade` introduces 3 *classes* of degrees of freedom (dofs). 

**Xdofs** are the dofs normaly encountered in normal "forward" FEM analysis.  They provide a discrete representation of the response of the system. There is a one-to-one relation (a "duality") between the ``X``dofs and the *residuals* ``R``, which are the discretized form of he differential equations we seek to solve. In forward FEM, we seek to solve a problem of the form

```math
\forall t, R(X(t),t)=0
```

**Udofs** are additional dofs that can be used to represent additional loads on the system. Like ``X``dofs, ``U``dofs are time-dependent. Unlike ``X``dofs, there is no residual (no new equations) corresponding to them.

**Adofs** are additional dofs that can be used to parameterise the model. Adofs are *not* time-dependant,and there are no residuals corresponing to them.

Equilibrium now requires

```math
\forall t, R(X(t),U(t),A,t)=0
```

which is ill-posed (unknowns ``U`` and ``A`` have been added, but without adding any new equation).

To simplify the presentations in the following, we drop time ``t`` from the notations.  This can be interpreted in two ways. Either we are looking at a problem at a single instant (static). Or we are solving an evolution problem (dynamic), and ``X`` and ``U`` are now vectors-valued functions of time, and ``R`` represents a vector of ordinary differential equations (over time). 

## Target function

To make the problem well-posed again, we introduce a target function ``Q(X,U,A)``.  We then seek to make ``Q`` stationary (finding a local minimum or maximum), under the constraint ``R(X,U,A)=0``. Depending on the relevant application, the target function ``Q`` can represent different concepts.

In design optimisation, one can associate a monetary value to ``A`` (building stronger costs more money), and to ``X`` (some system responses, including failure, would cost money). The objective is to find the design with the lowest associated cost.

In a load estimation problem where part of the response is measured, we wish to find the most probable unknown load ``U`` (large loads are not likely) and response (a computed response ``X`` that drasticaly disagrees with the actual measurements is not likely) - under the constraint that load and response verify equilibrium.  In this case, we minimize the *surprisal* ``s = -\log(P)`` instead of maximizing ``P``.  Maximizing probabilities or making surprisal stationary is equivalent. However, since probabilities of independent variables multiply, the corresponding surprisals add, which fits into the system of elements that *add* contribution to the system of equations to be solved. Further, the surprisals are typicaly numericaly "kinder".

## Lagrangian

Making ``Q`` stationary under the equilibrium constraints ``R(X,U,A)=0`` is equivalent to finding a stationary point (a saddle point) of the *Lagrangian*

```math
L(\Lambda,X,U,A) = Q(X,U,A) + R(X,U,A) \cdot \Lambda
```

where ``Λ`` are Lagrange multipliers.  

There is a one-to-one correspondance between Lagrange multipliers ``Λ``dofs and residuals, and hence between ``Λ``dofs and ``X``dofs.  One result of this correspondance is that when implmenting a new element, the method `doflist` does not list the ``Λ``dofs (this would otherwise just have been a compulsory repetition of the list of ``X``dofs). 

For evolution problems (involving a time-dependency), the dot product ``R(X,U,A) \cdot \Lambda`` includes an integral over time: the Lagrangian is a *functional*, and the gradients of ``L`` are ordinary (or integral) differential equations (in time), found using Frechet derivatives.


## Elements

In the original finite element analysis, the *elements* in a model form a partition of a domain over which differential equations are to be solved. We will call these "physical elements".  These elements provide additive contributions to the system of equations ``R(X,U,A)=0``. Further contributions can come from external loads.

When doing optimization-FEM in `Muscade`, the elements provide additive contributions to the Lagrangian scalar ``L``, instead of to the residual vector ``R``. Actually when creating a new element in `Muscade` to create an application, one can either implement a contribution to ``R`` by implementing a method `Muscade.residual` for the element, or a contribution to ``L`` by implementing `Muscade.lagrangian`.  For performance, whenever possible (for example when implementing a physical element), prefer `Muscade.residual`. 

In `Muscade`, the broad view is taken that anything that contributes to ``L`` is an element (even if no part of the domain is actually associated to the element).  This more general definition of "elements" includes a variety of types:

- Physical element, discretizing differential equations over a part of the domain
- Known external loads on the boundary (non-essential boundary conditions)
- Known external loads in the domain
- Constrained dofs (essential boundary conditions)
- Holonomic equality and inequality constraints (contact)
- Optimisation constraints (e.g. stresses shal not exceed some limit at any point within part of the domain)
- Response measurements (surprisal on Xdofs)
- Unknown external loads (surprisal on Udofs)
- Observed damage (surprisal on Adofs)
- Cost of unfavorable response (cost on Xdofs)
- Cost of actuators (cost on Udofs)
- Cost of building a system (cost on Adofs)

Because "everything" is an element, app developers can express a wide range of ideas through `Muscade`'s element API.

## Physical and optimisation constraints

*Physical constraints* (including contact or Dirichlet/essential boundary conditions) are added to the model by using an element that adds adds a dof of class ``X`` to the model.  This new dof ``X_λ`` is a Lagrange multiplier for the constraint. If we note ``R``, ``X`` and ``Λ`` the list of dofs before adding the constraint element and ``R^*``, ``X^*`` and ``Λ^*`` after, then  

```math
L^*(\Lambda^*,X^*,U,A) = Q(X,U,A)  + R^*(X^*,U,A) \cdot Λ^* 
```

with

```math
\begin{aligned}
Λ^* &= [Λ,Λ_λ]\\
X^* &= [X,X_λ]\\
R^*(X^*,U,A) &= \left[R(X,U,A) - ∇g_x(X,U,A) \cdot X_λ \; , \; g_x(X,U,A)\right]
\end{aligned}
```

*Optimisation constraints* allow to define that some situations are impossible, or inacceptable. For example, in a design optimisation analysis, excessive stresses would lead to failure, so would be constrained to remain under a given threshold. In a tracking or optimal control problem, an actuator force may not exceed some limit for the actuator's capacity.

Optimisation constraints that have to be verified at every step (for example stresses that must remain below a critical level, at any time) require a Lagrange multiplier that changes over time, and that is thus of class ``U``. Optimisation constraints that act only on ``A``dofs (for example, there is a limit to the strength of steel we can order)  require a Lagrange multiplier of class ``A``.  With optimisation constraints, the Lagrangian is of the form

```math
L^*(\Lambda,X,U^*,A^*) = Q^*(X,U^*,A^*)  + R(X,U,A) \cdot Λ 
```

with

```math
\begin{aligned}
U^* &= [U,U_λ]\\
A^* &= [A,A_λ]\\
Q^*(X,U^*,A^*) &= Q(X,U,A) +  g_u(X,U,A) \cdot U_\lambda + g_a(A) \cdot A_λ \\
\end{aligned}
```




