#include <cstdio>
#include <cmath>

// Custom contiguous 2D Matrix structure for cache locality
class Matrix2D {
public:
    int cols;
    int rows;
    double* data;

    Matrix2D(int c, int r) : cols(c), rows(r) {
        data = new double[cols * rows];
        setZero();
    }

    ~Matrix2D() {
        delete[] data;
    }

    inline double& operator()(int x, int y) {
        return data[y * cols + x];
    }

    inline const double& operator()(int x, int y) const {
        return data[y * cols + x];
    }

    void setZero() {
        for (int i = 0; i < cols * rows; ++i) {
            data[i] = 0.0;
        }
    }
};

// Main 2D CFD Solver Class
class CFDSolver2D {
private:
    int Nx, Ny;          // Number of fluid cells
    double dx, dy;       // Spatial steps
    double dt;           // Time step
    double nu;           // Kinematic viscosity
    double rho;          // Fluid density

    // Staggered MAC Grid fields
    Matrix2D u;          // Horizontal velocity: (Nx + 1) x Ny
    Matrix2D v;          // Vertical velocity:   Nx x (Ny + 1)
    Matrix2D p;          // Pressure field:       Nx x Ny

    // Intermediate fields
    Matrix2D u_star;
    Matrix2D v_star;
    Matrix2D rhs;        // Right-hand side of Poisson Eq

    // Conjugate Gradient workspace matrices
    Matrix2D r;
    Matrix2D d;
    Matrix2D q;

    // Enforce physical boundary conditions (Lid-driven cavity)
    void applyBoundaryConditions() {
        // 1. Left and Right Walls (No-slip, solid)
        for (int j = 0; j < Ny; ++j) {
            u(0, j) = 0.0;       // No flow through left wall
            u(Nx, j) = 0.0;      // No flow through right wall
        }
        for (int j = 0; j <= Ny; ++j) {
            // Ghost cell adjustments for tangential velocities
            // v at walls averaged with neighbors must equal 0
            if (j < Ny) {
                // Approximate v on left/right boundary cells
            }
        }

        // 2. Bottom Wall (No-slip, solid, stationary)
        for (int i = 0; i < Nx; ++i) {
            v(i, 0) = 0.0;       // No flow through bottom wall
        }

        // 3. Top Wall (Moving Lid: U = 1.0)
        for (int i = 0; i < Nx; ++i) {
            v(i, Ny) = 0.0;      // No flow through top wall
        }
        // Ghost cell interpolation for u at top/bottom boundaries
        // Handled directly inside advection/diffusion loops to ensure strict safety limits.
    }

    // Compute the matrix-vector product Ap for the Conjugate Gradient solver
    void computeLaplacian(const Matrix2D& pressureField, Matrix2D& laplacianOut) {
        double idx2 = 1.0 / (dx * dx);
        double idy2 = 1.0 / (dy * dy);

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                // Neumann boundary condition tracking: dP/dn = 0 at walls
                double p_left  = (i > 0)      ? pressureField(i - 1, j) : pressureField(i, j);
                double p_right = (i < Nx - 1) ? pressureField(i + 1, j) : pressureField(i, j);
                double p_down  = (j > 0)      ? pressureField(i, j - 1) : pressureField(i, j);
                double p_up    = (j < Ny - 1) ? pressureField(i, j + 1) : pressureField(i, j);

                laplacianOut(i, j) = (p_left - 2.0 * pressureField(i, j) + p_right) * idx2 +
                                     (p_down - 2.0 * pressureField(i, j) + p_up)    * idy2;
            }
        }
    }

    // Basic vector operations for the Conjugate Gradient Solver
    double dotProduct(const Matrix2D& m1, const Matrix2D& m2) {
        double result = 0.0;
        for (int i = 0; i < Nx * Ny; ++i) {
            result += m1.data[i] * m2.data[i];
        }
        return result;
    }

public:
    CFDSolver2D(int nx, int ny, double w, double h, double viscosity, double density)
        : Nx(nx), Ny(ny), nu(viscosity), rho(density),
          u(nx + 1, ny), v(nx, ny + 1), p(nx, ny),
          u_star(nx + 1, ny), v_star(nx, ny + 1), rhs(nx, ny),
          r(nx, ny), d(nx, ny), q(nx, ny) 
    {
        dx = w / Nx;
        dy = h / Ny;
        dt = 0.001; // Safe initial time-step
    }

    // Step 1: Compute intermediate fields using explicit updates
    void computeIntermediateVelocity() {
        // Horizontal intermediate tracking loop
        for (int j = 0; j < Ny; ++j) {
            for (int i = 1; i < Nx; ++i) {
                // Handle top/bottom boundary ghost values for u
                double u_top  = (j == Ny - 1) ? 2.0 * 1.0 - u(i, j)   : u(i, j + 1); // Lid speed = 1.0
                double u_bot  = (j == 0)      ? -u(i, j)              : u(i, j - 1);
                
                double diff_u = (u(i + 1, j) - 2.0 * u(i, j) + u(i - 1, j)) / (dx * dx) +
                                (u_top - 2.0 * u(i, j) + u_bot) / (dy * dy);

                double u_avg = u(i, j);
                double du_dx = (u(i + 1, j) - u(i - 1, j)) / (2.0 * dx);

                // Safe 4-point staggered cell interpolation for v components
                double v_avg = 0.25 * (v(i, j) + v(i, j + 1) + v(i - 1, j) + v(i - 1, j + 1));
                double du_dy = (u_top - u_bot) / (2.0 * dy);

                double advect_u = u_avg * du_dx + v_avg * du_dy;
                u_star(i, j) = u(i, j) + dt * (nu * diff_u - advect_u);
            }
        }

        // Vertical intermediate tracking loop
        for (int j = 1; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double v_left  = (i == 0)      ? -v(i, j)            : v(i - 1, j);
                double v_right = (i == Nx - 1) ? -v(i, j)            : v(i + 1, j);

                double diff_v = (v_right - 2.0 * v(i, j) + v_left) / (dx * dx) +
                                (v(i, j + 1) - 2.0 * v(i, j) + v(i, j - 1)) / (dy * dy);

                double u_avg = 0.25 * (u(i, j) + u(i, j - 1) + u(i + 1, j) + u(i + 1, j - 1));
                double dv_dx = (v_right - v_left) / (2.0 * dx);
                
                double v_avg = v(i, j);
                double dv_dy = (v(i, j + 1) - v(i, j - 1)) / (2.0 * dy);

                double advect_v = u_avg * dv_dx + v_avg * dv_dy;
                v_star(i, j) = v(i, j) + dt * (nu * diff_v - advect_v);
            }
        }
    }

    // Step 2: Compute Divergence of Intermediate Velocity Field
    void computePoissonRHS() {
        double factor = rho / dt;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double div = (u_star(i + 1, j) - u_star(i, j)) / dx +
                             (v_star(i, j + 1) - v_star(i, j)) / dy;
                rhs(i, j) = factor * div;
            }
        }
    }

    // Step 3: Pure Low-Level Conjugate Gradient (Linear System Solver)
    // Step 3: Pure Low-Level Conjugate Gradient (Linear System Solver)
    void solvePressurePoisson() {
        p.setZero(); // Initial guess
        
        // FIX: Removed ", Nx, Ny" from this call
        computeLaplacian(p, r); 

        // r = rhs - Ap
        for (int i = 0; i < Nx * Ny; ++i) {
            r.data[i] = rhs.data[i] - r.data[i];
            d.data[i] = r.data[i];
        }

        double delta_old = dotProduct(r, r);
        double delta_reason = delta_old;
        int max_iter = 500;
        double tol = 1e-6;

        for (int iter = 0; iter < max_iter; ++iter) {
            if (delta_old < tol * delta_reason) break;

            // FIX: Removed ", Nx, Ny" from this call as well
            computeLaplacian(d, q);
            
            double alpha = delta_old / dotProduct(d, q);

            for (int i = 0; i < Nx * Ny; ++i) {
                p.data[i] += alpha * d.data[i];
                r.data[i] -= alpha * q.data[i];
            }

            double delta_new = dotProduct(r, r);
            double beta = delta_new / delta_old;

            for (int i = 0; i < Nx * Ny; ++i) {
                d.data[i] = r.data[i] + beta * d.data[i];
            }
            delta_old = delta_new;
        }
    }

    // Step 4: Velocity Correction Mapping
    void correctVelocities() {
        double factor = dt / rho;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 1; i < Nx; ++i) {
                u(i, j) = u_star(i, j) - factor * (p(i, j) - p(i - 1, j)) / dx;
            }
        }
        for (int j = 1; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                v(i, j) = v_star(i, j) - factor * (p(i, j) - p(i, j - 1)) / dy;
            }
        }
    }

    // Execute one full temporal iteration
    void stepSimulation() {
        applyBoundaryConditions();
        computeIntermediateVelocity();
        computePoissonRHS();
        solvePressurePoisson();
        correctVelocities();
    }

    // Write internal engine metrics to structured standard ASCII VTK
    void writeVTK(int fileIndex) {
        char filename[64];
        std::sprintf(filename, "fluid_output_%04d.vtk", fileIndex);
        std::FILE* fp = std::fopen(filename, "w");
        if (!fp) return;

        std::fprintf(fp, "# vtk DataFile Version 3.0\n");
        std::fprintf(fp, "2D Navier-Stokes Fluid Cavity Simulation\n");
        std::fprintf(fp, "ASCII\n");
        std::fprintf(fp, "DATASET STRUCTURED_GRID\n");
        std::fprintf(fp, "DIMENSIONS %d %d 1\n", Nx, Ny);
        std::fprintf(fp, "POINTS %d float\n", Nx * Ny);

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                std::fprintf(fp, "%f %f 0.0\n", i * dx, j * dy);
            }
        }

        std::fprintf(fp, "POINT_DATA %d\n", Nx * Ny);
        std::fprintf(fp, "VECTORS velocity float\n");
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                // Interpolate staggered values to cell center point
                double u_center = 0.5 * (u(i, j) + u(i + 1, j));
                double v_center = 0.5 * (v(i, j) + v(i, j + 1));
                std::fprintf(fp, "%f %f 0.0\n", u_center, v_center);
            }
        }

        std::fprintf(fp, "SCALARS pressure float 1\n");
        std::fprintf(fp, "LOOKUP_TABLE default\n");
        for (int i = 0; i < Nx * Ny; ++i) {
            std::fprintf(fp, "%f\n", p.data[i]);
        }

        std::fclose(fp);
    }
};

int main() {
    // Keeping the grid size at 50x50 for speed
    CFDSolver2D solver(50, 50, 1.0, 1.0, 0.01, 1.0);

    std::printf("[CFD Engine Launched] Processing steps...\n");

    // Total steps are 20,000 (20.0 seconds of physical fluid time)
    int totalSteps = 20000; 
    
    // Saving a frame every 1000 steps so we generate 20 files total
    int saveInterval = 1000; 

    for (int step = 0; step <= totalSteps; ++step) {
        solver.stepSimulation();

        if (step % saveInterval == 0) {
            solver.writeVTK(step / saveInterval);
            std::printf(" Saved snapshot frame: %d / %d\n", step, totalSteps);
        }
    }

    std::printf("[Simulation Finished] Open the output .vtk files directly inside ParaView.\n");
    return 0;
}