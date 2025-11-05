/*
 * lambert_solver.h
 *
 * A high-speed C++ implementation of a Universal Variables Lambert Solver.
 * Based on the algorithm by Dario Izzo, "Revisiting Lambert's Problem" (2014).
 *
 * This is a single-header-file implementation. To use it:
 * 1. #include "lambert_solver.h"
 * 2. Call the `solve_lambert` function.
 */

#ifndef LAMBERT_SOLVER_H
#define LAMBERT_SOLVER_H

#include <cmath>
#include <stdexcept>
#include <iostream>

namespace LambertSolver {

    // --- Vector Operations ---
    // We'll use simple C-style arrays as inputs, matching SPICE.

    inline double vnorm_3d(const double v[3]) {
        return std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    }

    inline double vdot_3d(const double v1[3], const double v2[3]) {
        return v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
    }

    inline void vcross_3d(const double v1[3], const double v2[3], double out[3]) {
        out[0] = v1[1] * v2[2] - v1[2] * v2[1];
        out[1] = v1[2] * v2[0] - v1[0] * v2[2];
        out[2] = v1[0] * v2[1] - v1[1] * v2[0];
    }

    // --- Stumpff Functions ---
    // These are essential for the universal variable formulation.

    inline double stumpff_C2(double z) {
        if (z > 1e-6) {
            return (1.0 - std::cos(std::sqrt(z))) / z;
        } else if (z < -1e-6) {
            return (1.0 - std::cosh(std::sqrt(-z))) / z;
        } else {
            // Use Taylor series for small z
            return 0.5 - z / 24.0 + z * z / 720.0 - z * z * z / 40320.0;
        }
    }

    inline double stumpff_C3(double z) {
        if (z > 1e-6) {
            return (std::sqrt(z) - std::sin(std::sqrt(z))) / (z * std::sqrt(z));
        } else if (z < -1e-6) {
            return (std::sinh(std::sqrt(-z)) - std::sqrt(-z)) / (-z * std::sqrt(-z));
        } else {
            // Use Taylor series for small z
            return 1.0 / 6.0 - z / 120.0 + z * z / 5040.0 - z * z * z / 362880.0;
        }
    }

    /*
     * Solves Lambert's Problem using a universal variable (Izzo's) method.
     *
     * @param mu      Gravitational parameter of the central body (e.g., sun_gm).
     * @param r1_vec  Position vector of the departure body [3].
     * @param r2_vec  Position vector of the arrival body [3].
     * @param tof_sec Time of flight in seconds.
     * @param prograde If true, computes the prograde trajectory (default).
     * @param v1_out  [OUTPUT] The resulting departure velocity vector [3].
     * @param v2_out  [OUTPUT] The resulting arrival velocity vector [3].
     * @return        True on success, false on failure (e.g., no solution).
     */
    bool solve_lambert(const double mu, const double r1_vec[3], const double r2_vec[3],
                       double tof_sec, bool prograde,
                       double v1_out[3], double v2_out[3]) {
        
        double r1_mag = vnorm_3d(r1_vec);
        double r2_mag = vnorm_3d(r2_vec);

        double r1_dot_r2 = vdot_3d(r1_vec, r2_vec);

        double cos_delta_nu = r1_dot_r2 / (r1_mag * r2_mag);
        if (cos_delta_nu > 1.0) cos_delta_nu = 1.0;
        if (cos_delta_nu < -1.0) cos_delta_nu = -1.0;
        
        double delta_nu = std::acos(cos_delta_nu);

        // Determine prograde/retrograde transfer direction
        double h_vec[3];
        vcross_3d(r1_vec, r2_vec, h_vec);
        
        if (prograde) {
            if (h_vec[2] < 0) { // Z-component of angular momentum
                delta_nu = 2.0 * M_PI - delta_nu;
            }
        } else {
            if (h_vec[2] >= 0) {
                delta_nu = 2.0 * M_PI - delta_nu;
            }
        }
        
        double A = std::sqrt(r1_mag * r2_mag * (1.0 + cos_delta_nu));
        if (A == 0.0) {
            // Cannot solve (e.g., r1 or r2 is zero vector)
            return false;
        }

        // --- Root-Finding (Newton's Method) ---
        // We are solving for 'z' (the universal variable squared)
        
        double z = 0.0; // Initial guess
        double y, C2, C3, y_prime;
        double tof_error = 1.0;
        int max_iter = 100;
        int iter = 0;

        double sqrt_mu = std::sqrt(mu);

        for (iter = 0; iter < max_iter; ++iter) {
            C2 = stumpff_C2(z);
            C3 = stumpff_C3(z);

            double C3z = C3 * z;
            double C2z = C2 * z;

            y = r1_mag + r2_mag + A * (z * C3 - 1.0) / std::sqrt(C2);
            if (y < 0.0) {
                // No physical solution for this z, try a different guess
                z = z * 0.5; // Simple bisection step
                continue;
            }

            double x = std::sqrt(y / C2);
            double tof_calc = (x * x * x * C3 + A * std::sqrt(y)) / sqrt_mu;

            tof_error = tof_calc - tof_sec;

            if (std::abs(tof_error) < 1e-7) {
                // Converged!
                break;
            }

            // Calculate derivative for Newton's step
            if (z == 0.0) {
                y_prime = (std::sqrt(2.0)/40.0) * y + (A/8.0)*(std::sqrt(y) + A*std::sqrt(1.0/(2.0*y)));
            } else {
                y_prime = (y / (4.0 * z)) * (C2 - 3.0 * C3 / (2.0 * C2)) + 
                          (A / 8.0) * (3.0 * C3 * std::sqrt(y) / C2 + A * std::sqrt(C2 / y));
            }
            
            double tof_prime = y_prime / sqrt_mu;
            
            // Newton's step
            z = z - tof_error / tof_prime;
        }

        if (iter == max_iter) {
            // Failed to converge
            return false;
        }

        // --- Calculate Lagrange Multipliers ---
        double f = 1.0 - y / r1_mag;
        double g = A * std::sqrt(y / mu);
        double g_dot = 1.0 - y / r2_mag;
        double f_dot = (f * g_dot - 1.0) / g;

        // --- Calculate Velocity Vectors ---
        for (int i = 0; i < 3; ++i) {
            v1_out[i] = (r2_vec[i] - f * r1_vec[i]) / g;
            v2_out[i] = (g_dot * r2_vec[i] - r1_vec[i]) / g;
        }

        return true;
    }

} // namespace LambertSolver

#endif // LAMBERT_SOLVER_H
