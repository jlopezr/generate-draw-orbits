import sys
import matplotlib.pyplot as plt
import re
import matplotlib
from mpl_toolkits.basemap import Basemap
import numpy as np
import signal

# Use TkAgg backend for interactive plotting
matplotlib.use('TkAgg')

# Regular expression to extract the X, Y, and Z coordinates from the input
regex = re.compile(r"Position: \(X: ([\d\.-]+) m, Y: ([\d\.-]+) m, Z: ([\d\.-]+) m\)")

# Initialize lists to store the latitude and longitude for plotting
lats = []
lons = []

# Constants
R_EARTH = 6371000  # Radius of Earth in meters

# Set up the plot
plt.ion()  # Turn on interactive mode for dynamic updates
fig, ax = plt.subplots(figsize=(15, 10))
m = Basemap(projection='mill', ax=ax)

# Draw the Blue Marble image
#m.bluemarble()

# Draw the map boundary and fill the background with blue (for seas)
m.drawmapboundary(fill_color='aqua')

# Fill the continents with orange and set lake color to blue
m.fillcontinents(color='orange', lake_color='aqua')

# Draw coastlines and countries
m.drawcoastlines()
m.drawcountries()

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

# Register the signal handler
signal.signal(signal.SIGINT, on_close)

# Function to convert X, Y, Z to latitude and longitude
def xyz_to_latlon(x, y, z):
    lat = np.degrees(np.arcsin(z / R_EARTH))
    lon = np.degrees(np.arctan2(y, x))
    return lat, lon

# Function to check if the path crosses the 180th meridian
def crosses_180(lon1, lon2):
    return abs(lon1 - lon2) > 180

# Variable to store the last plotted red point
last_red_point = None

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

        # Convert to latitude and longitude
        lat, lon = xyz_to_latlon(x, y, z)

        # Append the new position to the lists
        if lons and crosses_180(lons[-1], lon):
            # If the path crosses the 180th meridian, start a new segment
            lats.append(np.nan)
            lons.append(np.nan)

        lats.append(lat)
        lons.append(lon)

        # Update the plot
        x_map, y_map = m(lons, lats)
        m.plot(x_map, y_map, 'bo-', markersize=1, label='Satellite Orbit')

        # Remove the last red point if it exists
        if last_red_point:
            last_red_point.remove()

        # Plot the last point with a larger marker size
        if len(lats) > 1:
            last_red_point, = m.plot(x_map[-1], y_map[-1], 'ro', markersize=5, label='Current Position')

        plt.draw()
        fig.canvas.flush_events()  # Force a redraw of the plot

# Show the final plot when the input ends
plt.ioff()
plt.show()