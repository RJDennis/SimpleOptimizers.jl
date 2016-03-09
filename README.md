# SimpleOptimizers.jl
A small collection of simple function optimizers.

Julia has a number of very good optimization packages.  Unlike those other packages, SimpleOptimizers.jl does not use JuMP.jl and it focuses on optimization methods that do not require derivatives.  Yes, it contains some of the usual optimizers: Newton-Raphson, Conjugate-gradient, etc, also with versions that impose box-constraints, but the main focus is on methods that are derivative free.  Naturally such methods are less efficient, but they are sometimes better behaved.

The derivative-free optimization methods included in this package are:

- Hooke-Jeeves
- Luus-Jaakola
- Nelder-Mead
- Genetic-search
- Mesh-search
- Pattern-search
- Quadratic interpolation

The derivative-based methods are:

- Newton-Raphson
- BFGS
- Conjugate gradient
