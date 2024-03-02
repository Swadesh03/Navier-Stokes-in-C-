#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>
#include <numeric>

class ductFlow
{
    public:
        ductFlow(); // constructor to inititate the object
        
        // call this to start the simulation
        void init();

    private:

        // to calculate b vector, later used in pressure_poisson function
        std::vector<std::vector<double>> build_up_b(std::vector<std::vector<double>> b,
                                           double rho, double dt,
                                           const std::vector<std::vector<double>>& u,
                                           const std::vector<std::vector<double>>& v,
                                           double dx, double dy);

        // to calculate the pressure term
        std::vector<std::vector<double>> pressure_poisson(std::vector<std::vector<double>>& p,
                                                  double dx, double dy,
                                                  const std::vector<std::vector<double>>& b,
                                                  int nit);

        // to calculate the u and v at each time step, and later update the boundary conditions u & v
        std::pair<std::vector<std::vector<double>>, std::vector<std::vector<double>>> update_uv_boundary_conditions(std::vector<std::vector<double>>& un,
                              std::vector<std::vector<double>>& vn,
                              const std::vector<std::vector<double>>& p,
                              double dt, double dx, double dy, double rho, double nu, double F);

        // to calculate the error between previous and the current time step
        double calculate_udiff(const std::vector<std::vector<double>>& u, const std::vector<std::vector<double>>& un);

        // to return the vector with defined number of points
        std::vector<double> linspace(double start, double end, int numPoints);

        // to flattern 2d vector in 1d vector, used for plotting the results
        std::vector<double> flatten(const std::vector<std::vector<double>>& input);

        // to plot the variable
        void plot_fun(const std::vector<double>& x,
                 const std::vector<double>& y,
                 const std::vector<double>& u,
                 const std::vector<double>& v);

        // to initiate the velocity at the inlet
        std::vector<std::vector<double>> vel_inlet(std::vector<std::vector<double>>& u, double ux, int ny);

        // to initiate the pressure at the inlet
        std::vector<std::vector<double>> pressure_inlet(std::vector<std::vector<double>>& p, double px, int ny);

        // to calculate the pressure at the outlet
        double pressure_outlet(std::vector<std::vector<double>>& p, int ny, int nx);

        // to calculate and plot the wall shear stress
        double wall_shear(const std::vector<double>& x, const std::vector<std::vector<double>>& u, double dy, int nx, double nu);

        // Constant variable declaration ---
        const int nx = 41;  // no of points in x
        const int ny = 41;  // no of points in y
        const int nt = 1000; // no of time steps
        const int nit = 50; // no of iteration to find the optimum value of pressure
        const double dx = 2.0 / (nx - 1);   // calculation of dx
        const double dy = 2.0 / (ny - 1);   // calculation of dy
        const double rho = 1; // density of fluid
        const double nu = 0.1;  // viscosity of fluid
        const double F = 1.0;   // source term in the x-momentum eqn to mimic the effect of pressure-driven flow
        const double dt = 0.01; // each time step length - 0.001s
        const double start = 0.0; // length of the channel
        const double end_x = 2.0; // width of the channel
        const double end_y = 2.0; // width of the channel

        double udiff = 1.0;
        int stepcount = 1;

};