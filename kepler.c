#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <string.h>

// Constants
const double G = 6.67430e-11;  // Gravitational constant (m^3 kg^-1 s^-2)
const double M = 5.97219e24;   // Mass of Earth (kg)
const double PI = 3.141592653589793;
const double R_EARTH = 6371000;  // Radius of Earth (meters)
const double OMEGA_EARTH = 7.2921159e-5;  // Earth's rotation rate (radians per second)

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

// Function to compute Greenwich Sidereal Time (GST) in radians
double compute_gst(double t) {
    return fmod(OMEGA_EARTH * t, 2 * PI);
}

// Function to simulate and generate coordinates from Keplerian elements
void simulate_orbit_kepler(double a, double e, double i, double omega, double Omega, double M0, double total_time, int n_steps, int wait) {
    double n = mean_motion(a);  // Mean motion
    double real_time_step = total_time / n_steps;  // Time step in real orbit time
    double desired_real_time_duration = 60.0;  // Desired real-time duration for the entire simulation (in seconds)
    double sleep_duration = desired_real_time_duration / n_steps;  // Sleep duration for each step (in seconds)

    printf("Simulating satellite orbit using Keplerian elements:\n");
    printf("Semi-major axis: %.2f meters\n", a);
    printf("Eccentricity: %.2f\n", e);
    printf("Inclination: %.2f degrees\n", i * 180.0 / PI);
    printf("Argument of periapsis (omega): %.2f degrees\n", omega * 180.0 / PI);
    printf("Longitude of ascending node (Omega): %.2f degrees\n", Omega * 180.0 / PI);
    printf("Mean anomaly at epoch (M0): %.2f degrees\n", M0 * 180.0 / PI);
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

        // Convert PQW coordinates to ECI (Earth-Centered Inertial) coordinates
        double cos_Omega = cos(Omega);
        double sin_Omega = sin(Omega);
        double cos_omega = cos(omega);
        double sin_omega = sin(omega);
        double cos_i = cos(i);
        double sin_i = sin(i);

        // Rotation matrix to convert from orbital plane to ECI
        double X_eci = (cos_Omega * cos_omega - sin_Omega * sin_omega * cos_i) * xp +
                       (-cos_Omega * sin_omega - sin_Omega * cos_omega * cos_i) * yp;

        double Y_eci = (sin_Omega * cos_omega + cos_Omega * sin_omega * cos_i) * xp +
                       (-sin_Omega * sin_omega + cos_Omega * cos_omega * cos_i) * yp;

        double Z_eci = (sin_omega * sin_i) * xp + (cos_omega * sin_i) * yp;

        // Compute Greenwich Sidereal Time (GST)
        double gst = compute_gst(t);

        // Convert ECI coordinates to ECEF coordinates
        double cos_gst = cos(gst);
        double sin_gst = sin(gst);

        double X_ecef = cos_gst * X_eci + sin_gst * Y_eci;
        double Y_ecef = -sin_gst * X_eci + cos_gst * Y_eci;
        double Z_ecef = Z_eci;

        // Output the current position
        printf("Step %d | Time: %.2f s | ECI Position: (x: %.2f m, y: %.2f m, z: %.2f m) | ECEF Position: (X: %.2f m, Y: %.2f m, Z: %.2f m)\n",
               step + 1, t, X_eci, Y_eci, Z_eci, X_ecef, Y_ecef, Z_ecef);

        fflush(stdout);

        // Sleep for real-time simulation if wait is enabled
        if (wait) {
            struct timespec ts;
            ts.tv_sec = (time_t)sleep_duration;
            ts.tv_nsec = (long)((sleep_duration - ts.tv_sec) * 1e9);  // Convert seconds to nanoseconds
            nanosleep(&ts, NULL);
        }
    }
}

void set_example_orbit(const char *example, double *a, double *e, double *i, double *omega, double *Omega, double *M0) {
    if (strcmp(example, "LEO") == 0) {
        *a = R_EARTH + 400000;
        *e = 0.01;
        *i = 51.6 * PI / 180.0;
        *omega = 0.0;
        *Omega = 0.0;
        *M0 = 0.0;
    } else if (strcmp(example, "GEO") == 0) {
        *a = R_EARTH + 35786000;
        *e = 0.0;
        *i = 0.0;
        *omega = 0.0;
        *Omega = 0.0;
        *M0 = 0.0;
    } else if (strcmp(example, "Molniya") == 0) {
        *a = 26560000;
        *e = 0.74;
        *i = 63.4 * PI / 180.0;
        *omega = 270.0 * PI / 180.0;
        *Omega = 0.0;
        *M0 = 0.0;
    } else if (strcmp(example, "SSO") == 0) {
        *a = R_EARTH + 600000;
        *e = 0.001;
        *i = 98.0 * PI / 180.0;
        *omega = 0.0;
        *Omega = 0.0;
        *M0 = 0.0;
    } else if (strcmp(example, "EquatorialCircular") == 0) {
        *a = R_EARTH + 400000;
        *e = 0.0;
        *i = 0.0;
        *omega = 0.0;
        *Omega = 0.0;
        *M0 = 0.0;
    } else if (strcmp(example, "EquatorialElliptic") == 0) {        
        double rp = R_EARTH + 400000;  // Perigee
        double ra = R_EARTH + 2000000; // Apogee

        *a = (rp + ra) / 2;
        *e = (ra - rp) / (ra + rp);
        *i = 0.0;
        *omega = 0.0;
        *Omega = 0.0;
        *M0 = 0.0;
    } else if (strcmp(example, "Inclined") == 0) {
        *a = R_EARTH + 4000000;
        *e = 0.5;
        *i = 45.0 * PI / 180.0;
        *omega = 0.0;
        *Omega = 0.0;
        *M0 = 0.0;
    } else {
        fprintf(stderr, "Unknown example orbit: %s\n", example);
        //Print available example orbits
        fprintf(stderr, "Available example orbits: LEO, GEO, Molniya, SSO, EquatorialCircular, EquatorialElliptic, Inclined\n");
        exit(EXIT_FAILURE);
    }
}

int main(int argc, char *argv[]) {
    // Default Keplerian elements
    double a = R_EARTH + 400000;  // Semi-major axis in meters (Earth radius + altitude)
    double e = 0.01;  // Eccentricity (0 for circular, >0 for elliptical)
    double i = 0.0 * PI / 180.0;  // Inclination in radians (0 for equatorial orbit)
    double omega = 0.0 * PI / 180.0;  // Argument of periapsis in radians
    double Omega = 0.0 * PI / 180.0;  // Longitude of ascending node in radians
    double M0 = 0.0 * PI / 180.0;  // Mean anomaly at epoch in radians
    double total_simulation_time = 5400.0;  // Simulated time in seconds (e.g., 1.5 hours)
    int n_steps = 360;  // Number of time steps
    int wait = 0;  // Default to no wait
    char *example = NULL;  // Example orbit

    // Define long options
    static struct option long_options[] = {
        {"semi-major-axis", required_argument, 0, 'a'},
        {"eccentricity", required_argument, 0, 'e'},
        {"inclination", required_argument, 0, 'i'},
        {"argument-of-periapsis", required_argument, 0, 'w'},
        {"longitude-of-ascending-node", required_argument, 0, 'O'},
        {"mean-anomaly", required_argument, 0, 'M'},
        {"total-time", required_argument, 0, 't'},
        {"steps", required_argument, 0, 'n'},
        {"wait", no_argument, 0, 'W'},
        {"example", required_argument, 0, 'x'},
        {0, 0, 0, 0}
    };

    // Parse command-line arguments
    int opt;
    int option_index = 0;
    while ((opt = getopt_long(argc, argv, "a:e:i:w:O:M:t:n:Wx:", long_options, &option_index)) != -1) {
        switch (opt) {
            case 'a':
                a = atof(optarg);
                break;
            case 'e':
                e = atof(optarg);
                break;
            case 'i':
                i = atof(optarg) * PI / 180.0;
                break;
            case 'w':
                omega = atof(optarg) * PI / 180.0;
                break;
            case 'O':
                Omega = atof(optarg) * PI / 180.0;
                break;
            case 'M':
                M0 = atof(optarg) * PI / 180.0;
                break;
            case 't':
                total_simulation_time = atof(optarg);
                break;
            case 'n':
                n_steps = atoi(optarg);
                break;
            case 'W':
                wait = 1;
                break;
            case 'x':
                example = optarg;
                break;
            default:
                fprintf(stderr, "Usage: %s [-a|--semi-major-axis semi_major_axis] [-e|--eccentricity eccentricity] [-i|--inclination inclination] [-w|--argument-of-periapsis argument_of_periapsis] [-O|--longitude-of-ascending-node longitude_of_ascending_node] [-M|--mean-anomaly mean_anomaly] [-t|--total-time total_time] [-n|--steps steps] [-W|--wait] [-x|--example example_orbit]\n", argv[0]);
                exit(EXIT_FAILURE);
        }
    }

    // Set example orbit if specified
    if (example) {
        set_example_orbit(example, &a, &e, &i, &omega, &Omega, &M0);
    }

    simulate_orbit_kepler(a, e, i, omega, Omega, M0, total_simulation_time, n_steps, wait);

    return 0;
}