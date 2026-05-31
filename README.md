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

Initial Pressure Distribution : 
<img width="1565" height="792" alt="InititalPressureDistribution" src="https://github.com/user-attachments/assets/b5b6fd34-aae6-44ad-99ab-106165e1d3c3" />

Final Pressure Distribution : 
<img width="1561" height="792" alt="FinalPressureDistribution" src="https://github.com/user-attachments/assets/2a23d4f8-dd42-4111-9342-cac2e7d55086" />
