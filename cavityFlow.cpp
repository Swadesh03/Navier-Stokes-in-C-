#include <iostream>
#include <fstream>
#include <vector>
#include <tuple>

// constant variable decalaration
const int nx = 41;  // no of points in x
const int ny = 41;  // no of points in y
const int nt = 100; // no of time steps
const int nit = 50; // no of iteration to find the optimum value of pressure
const double dx = 2.0 / (nx - 1);   // calculation of dx
const double dy = 2.0 / (ny - 1);   // calculation of dy
const double rho = 1.0; // density of fluid
const double nu = 0.1;  // viscosity of fluid
const double dt = 0.001;   // each time step length - 0.001s

/*  function -- "build_up_b"  #########################################################
    return b vector -- which can be later used in poisson equation
    b is just a part of poisson equation, seperated from poisson eqn
    to keep it short.
    input - b, rho, dt, u, v, dx, dy
    output - b vector
*/
std::vector<std::vector<double>> build_up_b(std::vector<std::vector<double>> b,
                                           double rho, double dt,
                                           const std::vector<std::vector<double>>& u,
                                           const std::vector<std::vector<double>>& v,
                                           double dx, double dy) {
    std::vector<std::vector<double>> result = b;

    for (int i = 1; i < ny - 1; ++i) {
        for (int j = 1; j < nx - 1; ++j) {
            result[i][j] = (rho * (1.0 / dt *
                ((u[i][j + 1] - u[i][j - 1]) / (2.0 * dx) + (v[i + 1][j] - v[i - 1][j]) / (2.0 * dy)) -
                ((u[i][j + 1] - u[i][j - 1]) / (2.0 * dx)) * ((u[i][j + 1] - u[i][j - 1]) / (2.0 * dx)) -
                2.0 * ((u[i + 1][j] - u[i - 1][j]) / (2.0 * dy) * (v[i][j + 1] - v[i][j - 1]) / (2.0 * dx)) -
                ((v[i + 1][j] - v[i - 1][j]) / (2.0 * dy)) * ((v[i + 1][j] - v[i - 1][j]) / (2.0 * dy))));
        }
    }

    return result;
}

/*  function -- "pressure poisson"  ###########################################################
    solves the poisson eqn to find out the pressure
    input - p, dx, dy, b, nit
    output - pressure vector
*/
std::vector<std::vector<double>> pressure_poisson(std::vector<std::vector<double>>& p,
                                                  double dx, double dy,
                                                  const std::vector<std::vector<double>>& b,
                                                  int nit) {
    std::vector<std::vector<double>> pn = p;

    for (int q = 0; q < nit; ++q) {
        pn = p;

        for (int i = 1; i < ny - 1; ++i) {
            for (int j = 1; j < nx - 1; ++j) {
                p[i][j] = ((pn[i][j + 1] + pn[i][j - 1]) * dy * dy +
                           (pn[i + 1][j] + pn[i - 1][j]) * dx * dx) /
                          (2.0 * (dx * dx + dy * dy)) -
                          (dx * dx * dy * dy) / (2.0 * (dx * dx + dy * dy)) * b[i][j];
            }
        }

        // Boundary conditions
        for (int i = 0; i < ny; ++i) {
            p[i][nx - 1] = p[i][nx - 2]; // dp/dx = 0 at x = 2
            p[i][0] = p[i][1];            // dp/dx = 0 at x = 0
        }
        for (int j = 0; j < nx; ++j) {
            p[0][j] = p[1][j];            // dp/dy = 0 at y = 0
            p[ny - 1][j] = 0.0;           // p = 0 at y = 2
        }
    }

    return p;
}


/*  function -- "cavity_flow"  ###########################################################
    solves the flow eqn in cavity
    input - u, v, dt, dx, dy, p, nu
    output - u, v, p
*/
std::tuple<std::vector<std::vector<double>>, std::vector<std::vector<double>>, std::vector<std::vector<double>>>
cavity_flow(int nt, std::vector<std::vector<double>> u, std::vector<std::vector<double>> v,
            double dt, double dx, double dy, std::vector<std::vector<double>>& p, double rho, double nu) {
    std::vector<std::vector<double>> un = u;
    std::vector<std::vector<double>> vn = v;
    std::vector<std::vector<double>> b(ny, std::vector<double>(nx, 0.0));

    for (int n = 0; n < nt; ++n) {
        un = u;
        vn = v;

        b = build_up_b(b,rho, dt, u, v, dx, dy);
        p = pressure_poisson(p, dx, dy, b, nit);

        for (int i = 1; i < ny - 1; ++i) {
            for (int j = 1; j < nx - 1; ++j) {
                u[i][j] = (un[i][j] - un[i][j] * dt / dx * (un[i][j] - un[i][j - 1]) -
                           vn[i][j] * dt / dy * (un[i][j] - un[i - 1][j]) -
                           dt / (2 * rho * dx) * (p[i][j + 1] - p[i][j - 1]) +
                           nu * (dt / (dx * dx) * (un[i][j + 1] - 2 * un[i][j] + un[i][j - 1]) +
                                 dt / (dy * dy) * (un[i + 1][j] - 2 * un[i][j] + un[i - 1][j])));

                v[i][j] = (vn[i][j] - un[i][j] * dt / dx * (vn[i][j] - vn[i][j - 1]) -
                           vn[i][j] * dt / dy * (vn[i][j] - vn[i - 1][j]) -
                           dt / (2 * rho * dy) * (p[i + 1][j] - p[i - 1][j]) +
                           nu * (dt / (dx * dx) * (vn[i][j + 1] - 2 * vn[i][j] + vn[i][j - 1]) +
                                 dt / (dy * dy) * (vn[i + 1][j] - 2 * vn[i][j] + vn[i - 1][j])));
            }
        }

        // Boundary conditions
        for (int i = 0; i < ny; ++i) {
            u[i][0] = 0;
            u[i][nx - 1] = 0;
            v[i][0] = 0;
            v[i][nx - 1] = 0;
        }
        for (int j = 0; j < nx; ++j) {
            u[0][j] = 0;
            u[ny - 1][j] = 1;  // set velocity on cavity lid equal to 1
            v[0][j] = 0;
            v[ny - 1][j] = 0;
        }
    }

    return std::make_tuple(u, v, p);
}

/*  function -- "linspace" ###################################################################
    a replication of linspace command in python
    input - start point, end points, number of points in between
    output - a vector
*/
std::vector<double> linspace(double start, double end, int numPoints) {
    std::vector<double> result(numPoints);
    double step = (end - start) / (numPoints - 1);
    for (int i = 0; i < numPoints; ++i) {
        result[i] = start + i * step;
    }
    return result;
}

int main() {

    int nx = 41;
    int ny = 41;

    // x & y variable declaration
    std::vector<double> x = linspace(0.0, 2.0, nx);
    std::vector<double> y = linspace(0.0, 2.0, ny);

    // Initialize variables
    std::vector<std::vector<double>> u(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> v(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> p(ny, std::vector<double>(nx, 0.0));
    std::vector<std::vector<double>> b(ny, std::vector<double>(nx, 0.0));
    int nt = 100;

    // Call cavity_flow function
    std::tie(u, v, p) = cavity_flow(nt, u, v, dt, dx, dy, p, rho, nu);

    std::ofstream file;

    file.open("x.txt");
    for (const double &val : x) {
        file << val << "\n";
    }
    file.close();

    file.open("y.txt");
    for (const double &val : y) {
        file << val << "\n";
    }
    file.close();

    // Save data to text files
    std::ofstream u_file("u.txt");
    std::ofstream v_file("v.txt");
    std::ofstream p_file("p.txt");

    for (int i = 0; i < ny; ++i) {
        for (int j = 0; j < nx; ++j) {
            u_file << u[i][j] << " ";
            v_file << v[i][j] << " ";
            p_file << p[i][j] << " ";
        }
        u_file << "\n";
        v_file << "\n";
        p_file << "\n";
    }

    u_file.close();
    v_file.close();
    p_file.close();

    return 0;
}
