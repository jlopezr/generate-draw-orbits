import sys
import matplotlib.pyplot as plt
import re
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Use TkAgg backend for interactive plotting
matplotlib.use('TkAgg')

# Regular expression to extract the X, Y, and Z coordinates from the input
regex = re.compile(r"Position: \(X: ([\d\.-]+) m, Y: ([\d\.-]+) m, Z: ([\d\.-]+) m\)")

# Initialize lists to store the X, Y, Z coordinates for plotting
x_vals = []
y_vals = []
z_vals = []

# Constants
R_EARTH = 6371000  # Radius of Earth in meters

# Set up the plot
plt.ion()  # Turn on interactive mode for dynamic updates
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
orbit_plot, = ax.plot([], [], [], 'bo-', label='Satellite Orbit', markersize=2)  # Line for the orbit with smaller markers

# Draw the Earth's surface as a sphere
u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
x_sphere = R_EARTH * np.cos(u) * np.sin(v)
y_sphere = R_EARTH * np.sin(u) * np.sin(v)
z_sphere = R_EARTH * np.cos(v)
ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color='orange', label='Earth Surface')

# Set plot limits (adjust as necessary)
ax.set_xlim(-7e6, 7e6)
ax.set_ylim(-7e6, 7e6)
ax.set_zlim(-7e6, 7e6)
ax.set_xlabel('X (meters)')
ax.set_ylabel('Y (meters)')
ax.set_zlabel('Z (meters)')
ax.set_title('Satellite Equatorial Orbit')
ax.legend()

# Flag to indicate if the window is closed
window_closed = False

# Function to handle window close event
def on_close(event):
    global window_closed
    print("Window closed")
    plt.close(fig)
    window_closed = True
    sys.exit(0)

# Connect the close event to the handler
fig.canvas.mpl_connect('close_event', on_close)

# Read from standard input in real-time
for line in sys.stdin:
    if window_closed:
        break

    # Search for the line containing the satellite's position
    match = regex.search(line)
    if match:
        x = float(match.group(1))
        y = float(match.group(2))
        z = float(match.group(3))

        print(f"X: {x}, Y: {y}, Z: {z}")

        # Append the new position to the lists
        x_vals.append(x)
        y_vals.append(y)
        z_vals.append(z)

        # Update the plot
        orbit_plot.set_data(x_vals, y_vals)
        orbit_plot.set_3d_properties(z_vals)
        ax.relim()  # Recalculate limits
        ax.autoscale_view(True, True, True)  # Autoscale the plot to fit new data
        plt.draw()
        fig.canvas.flush_events()  # Force a redraw of the plot
    # else:
    #     print("Error parsing line: " + line)

# Show the final plot when the input ends
plt.ioff()
plt.show()