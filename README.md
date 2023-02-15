# NumericalMethods
This repository contains files to complete a variety of numerical algorithms relevant to computational physics and chemistry.

## Integrals
- integrals.py --> functions to calculate numerical integrals using the trapezoidal rule, Simpson's rule, Romberg, or Gaussian quadrature method
- gaussxw.py --> calculates the integration points and Gaussian weights for performing the Gaussian quadrature method
- gaussquad_ex.py --> an example problem using Gaussian quadrature

## Linear Systems of Equations
- Substitutions.py --> functions to perform backward or forward substitution
- GaussElim.py --> performs the Gaussian elimination algorithm with partial pivoting
- Decompositions.py --> functions to perform the LU (with or without partial pivoting) and QR decompositions
- LUsolve.py --> performs the LU algorithm with or without partial pivoting
- eigQR.py --> performs the QR algorithm for calculating the eigenvalues and eigenvectors of a matrix

## Ordinary and partial differential equations
- RungeKutta_ODE.py --> performs the fourth-order Runge-Kutta method for first-order ODEs
- SecOrder_ODEs.py --> performs the fourth-order Runge-Kutta method for a second-order ODE rewritten as two coupled first-order ODEs
- Verlet_ODE.py --> performs the Verlet method for ODEs
- BulirschStoer_ODE.py --> performs the Bulirsch-Stoer method for ODEs
- JacobiORGS_PDE.py --> performs the Jacobi method using overrelaxation/Gauss-Seidel speed up for PDEs
- FTCS_PDE.py --> performs the FTCS method for PDEs

## Monte Carlo
- MonteCarlo.py --> performs as Monte Carlo simulation for the decay of (213)Bi as an example
- MCintegration.py --> performs Monte Carlo integration for the volume of a 10D unit hypersphere as an example
- MCIimportancesampling.py --> performs Monte Carlo integration with importance sampling on an example integral
- SimulatedAnnealing.py --> performs the simulated annealing method to find the global minimum of a function with two examples provided