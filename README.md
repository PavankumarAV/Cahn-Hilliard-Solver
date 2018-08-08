# Cahn-Hilliard-Solver
This can be used to explore solution of 1-D Cahn-Hilliard Equation using PDE and explore metastable states using Continuously adaptive switching numerical scheme.
What solver can do:
  1. Compute numerical solution of Cahn-Hilliard equation for any initial condition and display as a animation/movie. The solution is       constructed/approximated using spectral methods .
  2. Compute numerical solution using ODE obtained in [1] for fast transitions and display as
  animation.
  3. Explore metastable states ( use switching scheme). It used both PDE and ODE with continuous
  adaptation of time, and display as animation.
What it can not do:
  1. Numerical blow up if the initial condition is not monotone stationary solution ( for example
  h = [0.2, 0.8] h = [0.3, 0.7].
  2. Switching scheme cannot take more than 4 transition layers (But easily extendable to N−
  transition layers).
