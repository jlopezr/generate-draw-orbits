import sys
import matplotlib.pyplot as plt
import re
import matplotlib

# Use TkAgg backend for interactive plotting
matplotlib.use('TkAgg')

# Regular expression to extract the X, Y, and Z coordinates from the input
regex = re.compile(r"Position: \(X: ([\d\.-]+) m, Y: ([\d\.-]+) m, Z: ([\d\.-]+) m\)")

# Initialize lists to store the X, Y coordinates for plotting
x_vals = []
y_vals = []

# Constants
R_EARTH = 6371000  # Radius of Earth in meters

# Set up the plot
plt.ion()  # Turn on interactive mode for dynamic updates
fig, ax = plt.subplots()
orbit_plot, = ax.plot([], [], 'bo-', label='Satellite Orbit', markersize=2)  # Line for the orbit with smaller markers
last_point_plot = ax.scatter([], [], color='red', s=50, label='Last Point')  # Scatter plot for the last point

# Draw the Earth's surface as a circle
earth_circle = plt.Circle((0, 0), R_EARTH, color='orange', fill=False, label='Earth Surface')
ax.add_artist(earth_circle)

# Set plot limits (adjust as necessary)
ax.set_xlim(-7e6, 7e6)
ax.set_ylim(-7e6, 7e6)
ax.set_aspect('equal', 'box')
ax.set_xlabel('X (meters)')
ax.set_ylabel('Y (meters)')
ax.set_title('Satellite Equatorial Orbit (View from North Pole)')
ax.grid(True)
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

        # Update the plot
        orbit_plot.set_data(x_vals, y_vals)
        last_point_plot.set_offsets([[x_vals[-1], y_vals[-1]]])  # Update the last point
        ax.relim()  # Recalculate limits
        ax.autoscale_view(True, True, True)  # Autoscale the plot to fit new data
        plt.draw()
        fig.canvas.flush_events()  # Force a redraw of the plot
    # else:
    #     print("Error parsing line: " + line)

# Show the final plot when the input ends
plt.ioff()
plt.show()