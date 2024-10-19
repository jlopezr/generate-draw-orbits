#include <stdio.h>
#include <math.h>
#include <time.h>

// Constants
const double G = 6.67430e-11;  // Gravitational constant (m^3 kg^-1 s^-2)
const double M = 5.97219e24;   // Mass of Earth (kg)
const double PI = 3.141592653589793;
const double R_EARTH = 6371000;  // Radius of Earth (meters)

// Function to compute the mean motion (n)
double mean_motion(double semi_major_axis) {
    return sqrt(G * M / pow(semi_major_axis, 3));  // radians per second
}

// Function to solve Kepler's equation using Newton's method
double solve_keplers_equation(double mean_anomaly, double eccentricity) {
    double E = mean_anomaly;  // Initial guess: E = M
    double tolerance = 1e-6;  // Tolerance for the solution
    double delta;

    do {
        delta = E - eccentricity * sin(E) - mean_anomaly;
        E = E - delta / (1 - eccentricity * cos(E));
    } while (fabs(delta) > tolerance);

    return E;
}

// Function to simulate and generate coordinates from Keplerian elements
void simulate_orbit_kepler(double a, double e, double i, double omega, double Omega, double M0, double total_time, int n_steps) {
    double n = mean_motion(a);  // Mean motion
    double real_time_step = total_time / n_steps;  // Time step in real orbit time

    printf("Simulating satellite orbit using Keplerian elements:\n");
    printf("Semi-major axis: %.2f meters\n", a);
    printf("Eccentricity: %.2f\n", e);
    printf("Inclination: %.2f degrees\n", i);
    printf("Argument of periapsis (omega): %.2f degrees\n", omega);
    printf("Longitude of ascending node (Omega): %.2f degrees\n", Omega);
    printf("Mean anomaly at epoch (M0): %.2f degrees\n", M0);
    printf("Total simulated time: %.2f seconds\n", total_time);
    printf("Time step: %.2f seconds\n\n", real_time_step);

    // Loop over time steps to compute the satellite's position
    for (int step = 0; step < n_steps; step++) {
        double t = step * real_time_step;  // Time
        double M = M0 + n * t;  // Mean anomaly at time t (radians)
        double E = solve_keplers_equation(M, e);  // Solve for eccentric anomaly (radians)

        // Compute true anomaly (nu) from eccentric anomaly (E)
        double nu = 2 * atan2(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2));

        // Orbital radius
        double r = a * (1 - e * cos(E));

        // Position in the orbital plane (PQW coordinates)
        double xp = r * cos(nu);  // X position in orbital plane
        double yp = r * sin(nu);  // Y position in orbital plane
        //double zp = 0.0;  // Z position is always 0 in the orbital plane

        // Convert PQW coordinates to ECI (Earth-Centered Inertial) coordinates
        double cos_Omega = cos(Omega);
        double sin_Omega = sin(Omega);
        double cos_omega = cos(omega);
        double sin_omega = sin(omega);
        double cos_i = cos(i);
        double sin_i = sin(i);

        // Rotation matrix to convert from orbital plane to ECI
        double X = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * xp +
                   (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * yp;

        double Y = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * xp +
                   (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * yp;

        double Z = (sin_omega * sin_i) * xp + (cos_omega * sin_i) * yp;

        // Output the current position
        printf("Step %d | Time: %.2f s | Position: (X: %.2f m, Y: %.2f m, Z: %.2f m)\n",
               step + 1, t, X, Y, Z);

        // Sleep for real-time simulation
        struct timespec ts;
        ts.tv_sec = 0;
        ts.tv_nsec = (long)(real_time_step * 1e9);  // Convert seconds to nanoseconds
        nanosleep(&ts, NULL);
    }
}

int main() {
    // Example Keplerian elements
    double a = R_EARTH + 400000;  // Semi-major axis in meters (Earth radius + altitude)
    double e = 0.01;  // Eccentricity (0 for circular, >0 for elliptical)
    double i = 0.0 * PI / 180.0;  // Inclination in radians (0 for equatorial orbit)
    double omega = 0.0 * PI / 180.0;  // Argument of periapsis in radians
    double Omega = 0.0 * PI / 180.0;  // Longitude of ascending node in radians
    double M0 = 0.0 * PI / 180.0;  // Mean anomaly at epoch in radians
    double total_simulation_time = 5400.0;  // Simulated time in seconds (e.g., 1.5 hours)
    int n_steps = 360;  // Number of time steps

    simulate_orbit_kepler(a, e, i, omega, Omega, M0, total_simulation_time, n_steps);

    return 0;
}
