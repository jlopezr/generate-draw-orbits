import subprocess
import matplotlib.pyplot as plt
import re
import matplotlib

# Use TkAgg backend for interactive plotting
matplotlib.use('TkAgg')

# Regular expression to extract the X, Y, and Z coordinates from the output
regex = re.compile(r"Position: \(X: ([\d\.-]+) m, Y: ([\d\.-]+) m, Z: ([\d\.-]+) m\)")

# Initialize lists to store the X, Y coordinates for plotting
x_vals = []
y_vals = []

# Constants
R_EARTH = 6371000  # Radius of Earth in meters

# Set up the plot
plt.ion()  # Turn on interactive mode for dynamic updates
fig, ax = plt.subplots()
orbit_plot, = ax.plot([], [], 'bo-', label='Satellite Orbit')  # Line for the orbit

# Draw the Earth's surface as a circle
earth_circle = plt.Circle((0, 0), R_EARTH, color='orange', fill=False, label='Earth Surface')
ax.add_artist(earth_circle)

# Set plot limits (adjust as necessary)
ax.set_xlim(-7e6, 7e6)
ax.set_ylim(-7e6, 7e6)
ax.set_xlabel('X (meters)')
ax.set_ylabel('Y (meters)')
ax.set_title('Satellite Equatorial Orbit')
ax.grid(True)
ax.legend()

# Run the C program and capture its output in real-time
with subprocess.Popen(['./satellite'], stdout=subprocess.PIPE, text=True) as proc:
    for line in proc.stdout:
        # Search for the line containing the satellite's position
        match = regex.search(line)
        if match:
            x = float(match.group(1))
            y = float(match.group(2))
            z = float(match.group(3))

            print(f"X: {x}, Y: {y}, Z: {z}")

            # # Append the new position to the lists
            x_vals.append(x)
            y_vals.append(y)

            # Update the plot
            orbit_plot.set_data(x_vals, y_vals)
            ax.relim()  # Recalculate limits
            ax.autoscale_view(True, True, True)  # Autoscale the plot to fit new data
            plt.draw()
            plt.pause(0.01)  # Pause briefly to allow the plot to update
        else:
            print("Error parsing line: " + line)

# Show the final plot when the simulation ends
plt.ioff()
plt.show()