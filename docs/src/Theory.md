# Theory

## Classes of degrees of freedom

`Muscade` introduces 3 *classes* of degrees of freedom (dofs). 

**``X``-dofs** are the dofs normaly encountered in normal "forward" FEM analysis.  They provide a discrete representation of the response of the system. There is a one-to-one relation (a "duality") between the ``X``-dofs and the *residuals* ``R``, which are the discretized form of the differential equations we seek to solve. In forward FEM, we formulate a discreet problem of the form

```math
\forall t, R(X(t),t)=0
```

**``U``-dofs** are additional dofs that can be used to represent additional *unknown* loads on the system. Like ``X``-dofs, ``U``-dofs are time-dependent. Unlike ``X``-dofs, there is no residual (no new equations) corresponding to them.

**``A``-dofs** are additional dofs that can be used to represent *unknown* model parameters. ``X``-dofs are *not* time-dependent, and there are no residuals corresponing to them.

We now formulate the discrete finite element model of the form

```math
\forall t, R(X(t),U(t),A,t)=0
```

which is ill-posed (unknowns ``U`` and ``A`` have been added, but without adding any new equation).

To simplify the presentations in the following, we drop time ``t`` from the notations.  This can be interpreted in two ways:

1) We are solving a problem at a single instant (static). 
2) We are solving an evolution problem (dynamic).  ``X``, ``U`` and ``R`` are now vectors-valued functions of time, and ``R=0`` denotes that the discretied differential equations must be verified at all times. 

## Target function

To make the problem well-posed again, we introduce a target function ``Q(X,U,A)``.  We then seek to make ``Q`` stationary (finding a local minimum or maximum), under the constraint ``R(X,U,A)=0``. Depending on the relevant application, the target function ``Q`` can represent different concepts.

### Target function as financial value
In design optimisation, one can associate a monetary value to ``A`` (building stronger costs more money), and to ``X`` (some system responses, including failure, would cost money). The objective is to find the design with the lowest total cost,- under the constraint that load and response verify equilibrium.

### Target function as surprisal
In a load estimation problem where part of the response is measured, we wish to find the most probable unknown load ``U`` (very large loads are not likely) and response (a computed response ``X`` that drasticaly disagrees with the actual measurements is not likely) - under the constraint that load and response verify equilibrium.  

In Muscade, this *must* be handled by minimizing the *surprisal* ``Q = -\log(P)`` instead of maximizing the probability density ``P``. 

1) ``Q`` and ``P`` have the same extrema, so finding a minimum of ``Q`` does give a maximum of ``P``.  
2) The joint probability density of independant random variables is equal to the *product* of each variable's probability density. Because ``\log(a+b)=\log(a)+\log(b)``, the joint surprisal of independant random variables is equal to the *sum* of each variable's surprisal. In forward finite element analysis, the load vectors and incremental matrices of elements are *added* into a system vector and system matrix.  Muscade extends this logic: contributions to the Lagrangian from various elements are *added* together.
3) Generaly speaking, surprisals are numericaly better behaved that probability densities.  
 
To illustrate the last point, consider a multinormal Gaussian probability density distribution and its surprisal:
   
```math
\begin{aligned}
P(X) &= (2\pi)^{-k/2} \det(Σ)^{-1/2} \exp \left(-\frac{1}{2} \left( X-\mu \right)^T \cdot \Sigma^{-1} \cdot \left( X-\mu \right) \right) \\
Q(X) &= -\log(P(X)) \\
     &= k + \frac{1}{2} \left( X-\mu \right)^T \cdot Q_{XX} \cdot \left( X-\mu \right)
\end{aligned}
```

where ``k`` is a constant, whose value does not affect the extremal ``X``, so ``k`` will be ignored in contributions to the target function. The constant matrix ``Q_{XX}`` is the inverse of the covariance matrix ``\Sigma``. In this example, ``s`` is a quadratic function of ``X``: its Hessian is constant and has value ``Q_{XX}``.  A Newton-Rapshon algorithm solves the minimization of a quadratic function exactly, in one iteration.  Further, a Newton-Raphson algorithm solves the minimization of a quadratic function, constrained by linear constraints, exaclty, in one iteration.  By contrast, a Newton-Raphson algorithm applied to finding the zero of the derivative of a Gaussian distribution, is prone to diverge.

Hence, the code for an element defining a "cost" on a single dof ``X``, measured to a value ``v(t)`` with a Gaussian measurement error of standard deviation ``\sigma`` would implemented the function 

```math
Q(X,t) = \frac{1}{2} \frac{1}{\sigma^2} \left( X-v(t) \right)^2
```

### Contributions to the target function

In Muscade, additive contributions to the target function are called "costs", although this is admitedly an inadequate name for an additive contribution to a surprisal. The costs are implemented as elements.  

One example of built-in element defining a cost is [`DofCost`](@ref), an element to add a cost which is a function of the value of a dof. Examples would include:

- A cost on a ``U``-dof, a surprisal expressing prior knowledge of the magnitude of forces that may be explaining measured response.
- A cost on a ``X``-dof, a financial cost incured if a point in the model has a too high value, causing a failure.

Another built-in element is [`ElementCost`](@ref), an element to add a cost which is a function of an element-result. For example, a bar element can have "axial strain" as an internal value, that is, a value it computes that is neither a dof nor a residual. 

## Constrained optimization

Making ``Q`` stationary under the equilibrium constraints ``R(X,U,A)=0`` is equivalent to finding a stationary point (a saddle point) of the *Lagrangian* ``L``

```math
L(\Lambda,X,U,A) = Q(X,U,A) + R(X,U,A) \cdot \Lambda
```

where ``Λ`` are [Lagrange multipliers](https://en.wikipedia.org/wiki/Lagrange_multiplier) (also known as *adjoint state variables*).  

There is a one-to-one correspondance between Lagrange multipliers ``Λ``dofs and residuals  ``R``, and hence between ``Λ``-dofs and ``X``-dofs.  One result of this correspondance is that when implementing a new element, the method that must be provided `doflist` does not list the ``Λ``-dofs (this would otherwise just have been a compulsory repetition of the list of ``X``-dofs). 

For evolution problems (involving a time-dependency), the dot product ``R(X,U,A) \cdot \Lambda`` includes an integral over time: the Lagrangian is a *functional*, and the gradients of ``L`` are ordinary differential equations in time, found using [functional derivatives](https://en.wikipedia.org/wiki/Functional_derivative).  

## Elements

In the original finite element analysis, the *elements* in a model form a partition of a domain over which differential equations are to be solved. We will call these "physical elements".  These elements provide additive contributions to the system of equations ``R(X,U,A)=0``. Further contributions can come from external loads.

When doing optimization-FEM in `Muscade`, the elements provide additive contributions to the Lagrangian scalar ``L``, instead of to the residual vector ``R``. Actually when creating a new element in `Muscade` to create an application, one can either implement a contribution to ``R`` by implementing a method `Muscade.residual` for the element, or a contribution to ``L`` by implementing `Muscade.lagrangian`.  For performance, whenever possible (for example when implementing a physical element), prefer `Muscade.residual`. 

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

Optimisation constraints that have to be verified at every step (for example stresses that must remain below a critical level, at any time) require a Lagrange multiplier that changes over time, and that is thus of class ``U``. Optimisation constraints that act only on ``A``-dofs (for example, there is a limit to the strength of steel we can order)  require a Lagrange multiplier of class ``A``.  With optimisation constraints, the Lagrangian is of the form

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

## Literature

For more details on the theory and example of applications, see [maincon04a](@cite), [maincon04b](@cite), [maree04](@cite), [barnardo04](@cite), [hauser06](@cite), [wu07](@cite), [maincon08](@cite), [hauser08](@cite), [wu08](@cite), [wu09](@cite), [maincon13](@cite). 



