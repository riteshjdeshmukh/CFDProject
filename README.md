# 2D Navier-Stokes CFD Solver

A custom Computational Fluid Dynamics (CFD) engine written from scratch in standard C++. This project simulates the classic **Lid-Driven Cavity Flow** problem using Chorin's Projection Method.

## Features
* **Zero External Dependencies:** Built entirely using standard C++ library components (`<cmath>`, `<cstdio>`, `<vector>`).
* **Custom Memory Management:** Uses a flat 1D array mapped to 2D space for optimal CPU cache locality.
* **Algorithm:** Implements a staggered MAC grid, Forward Euler time integration, and a Conjugate Gradient linear solver for the pressure Poisson equation.
* **Visualization:** Exports raw velocity and pressure field data sequentially into `.vtk` formats, ready for rendering in ParaView.

## How to Compile and Run
Ensure you have a C++ compiler (`g++`, `clang++`, or MSVC) installed.

1. Clone the repository:
   `git clone https://github.com/YourUsername/CFDProject.git`
2. Compile the code with O3 optimization:
   `g++ -O3 cfd.cpp -o cfd_solver`
3. Run the executable:
   `./cfd_solver` (or `cfd_solver.exe` on Windows)
4. Open the resulting `fluid_output_*.vtk` files in ParaView to visualize the vortex.
