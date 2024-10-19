#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>

// Constants
const double G = 6.67430e-11;  // Gravitational constant (m^3 kg^-1 s^-2)
const double M = 5.97219e24;   // Mass of Earth (kg)
const double R_EARTH = 6371000;  // Radius of Earth (meters)
const double ALTITUDE = 400000;  // Altitude of satellite above Earth's surface (meters)
const double SIMULATION_TIME = 60.0;  // Simulation time for one complete orbit (seconds)
const int N_STEPS = 360;        // Number of steps for the simulation

// Function to simulate orbit and print position (X, Y, Z) and velocity to stdout
void simulate_orbit(int wait) {
    double r = R_EARTH + ALTITUDE;  // Total distance from Earth's center to satellite (meters)
    double velocity = sqrt(G * M / r);  // Orbital velocity (m/s)
    double real_orbital_period = 5520.0;  // Real orbital period (seconds), 1h32m = 5520 seconds
    double time_compression_factor = real_orbital_period / SIMULATION_TIME;  // Compression factor
    double real_time_step = real_orbital_period / N_STEPS;  // Real time between each step (seconds)
    double simulated_time_step = real_time_step / time_compression_factor;  // Simulated time step (seconds)

    printf("Simulating satellite in equatorial orbit:\n");
    printf("Orbital radius: %.2f meters\n", r);
    printf("Orbital velocity: %.2f meters/second\n", velocity);
    printf("Real orbital period: %.2f seconds\n", real_orbital_period);
    printf("Time compression factor: %.2f\n", time_compression_factor);
    printf("Simulated time step: %.3f seconds\n\n", simulated_time_step);

    // Main simulation loop
    double simulated_time = 0.0;  // Initialize the simulated time
    for (int i = 0; i < N_STEPS; i++) {
        double time = i * real_time_step;  // Real orbital time
        double angle = 2 * M_PI * (time / real_orbital_period);  // Angle in radians
        double x = r * cos(angle);  // X-coordinate (meters)
        double y = r * sin(angle);  // Y-coordinate (meters)
        double z = 0.0;  // Z-coordinate (meters), always 0 for equatorial orbit

        // Print the current simulated time, position (X, Y, Z), and velocity
        printf("Simulated Time: %.2f s | Position: (X: %.2f m, Y: %.2f m, Z: %.2f m) | Velocity: %.2f m/s\n",
               simulated_time, x, y, z, velocity);

        // Increment the simulated time
        simulated_time += simulated_time_step * time_compression_factor;

        // Sleep to simulate real time if wait is true
        if (wait) {
            struct timespec ts;
            ts.tv_sec = 0;
            ts.tv_nsec = (long)(simulated_time_step * 1e9);  // Convert seconds to nanoseconds
            nanosleep(&ts, NULL);
        }
    }
}

int main(int argc, char *argv[]) {
    int wait = 0;  // Default to not waiting

    // Check for the --wait argument
    if (argc > 1 && strcmp(argv[1], "--wait") == 0) {
        wait = 1;
    }

    simulate_orbit(wait);
    return 0;
}