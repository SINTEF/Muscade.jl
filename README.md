# Muscade.jl

Muscade.jl is work in progress, and the ambition is to have a documented version with enough functionality and performance, fit for public use by the end of 2023.

Muscade.jl is a framework for "FEM optimisation", initiated by SINTEF, and made available under an MIT license.  

Muscade.jl solves constrained optimisation problems, where the constraints include that a FEM model must be at equilibrium.  The target function to be optimised could be a cost ("what is the cheapest design that survives all the load cases?") or the logarithm of a probability density ("what loads and/or model adjustements would explain the measured response?").

To build applications using Muscade one defines the necessary finite-element types, including definitions of costs or surprisal to handle specific problems.  It is not foreseen to create a comprehensive library of elements within Muscade.jl, this is left to other packages.

Emphasis is on performance, by encouraging element code that is typestable and uses memory on the stack.  Parallelisation will be supported, first by using multithreading.  Emphasis is also on making it easy to create applications. Automatic differentiation makes it unecessary to formulate and implement incremental forms, or gradients of target functions, and EspyInsideFunction.jl is used to make element-internal results available to the user, without taxing the developer with trivial work.  

