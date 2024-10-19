import sys
import matplotlib.pyplot as plt
import re
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import argparse

# Use TkAgg backend for interactive plotting
matplotlib.use('TkAgg')

# Argument parser setup
parser = argparse.ArgumentParser(description='Plot satellite orbit in 3D.')
parser.add_argument('--follow', action='store_true', help='Follow the satellite with the view')
args = parser.parse_args()

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
last_point_plot = ax.scatter([], [], [], color='red', s=50, label='Last Point')  # Scatter plot for the last point

# Draw the Earth's surface as a sphere
u, v = np.mgrid[0:2*np.pi:100j, 0:np.pi:50j]
x_sphere = R_EARTH * np.cos(u) * np.sin(v)
y_sphere = R_EARTH * np.sin(u) * np.sin(v)
z_sphere = R_EARTH * np.cos(v)
ax.plot_wireframe(x_sphere, y_sphere, z_sphere, color='orange', alpha=0.3, label='Earth Surface')

# Draw the equator
theta = np.linspace(0, 2 * np.pi, 100)
x_equator = R_EARTH * np.cos(theta)
y_equator = R_EARTH * np.sin(theta)
z_equator = np.zeros_like(theta)
ax.plot(x_equator, y_equator, z_equator, color='green', label='Equator')

# Draw the Greenwich meridian
phi = np.linspace(-np.pi / 2, np.pi / 2, 100)
x_meridian = R_EARTH * np.cos(phi)
y_meridian = np.zeros_like(phi)
z_meridian = R_EARTH * np.sin(phi)
ax.plot(x_meridian, y_meridian, z_meridian, color='yellow', linewidth=2, label='Greenwich Meridian')

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
        last_point_plot._offsets3d = (x_vals[-1:], y_vals[-1:], z_vals[-1:])  # Update the last point

        if args.follow:
            # Calculate the azimuth and elevation angles
            azim = np.degrees(np.arctan2(y, x))
            elev = np.degrees(np.arctan2(z, np.sqrt(x**2 + y**2)))
            ax.view_init(elev=elev, azim=azim)

        ax.relim()  # Recalculate limits
        ax.autoscale_view(True, True, True)  # Autoscale the plot to fit new data
        plt.draw()
        fig.canvas.flush_events()  # Force a redraw of the plot

# Show the final plot when the input ends
plt.ioff()
plt.show()