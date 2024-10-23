import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# 1. Create a figure and axis
fig, ax = plt.subplots()
x = np.linspace(0, 2 * np.pi, 100)
line, = ax.plot(x, np.sin(x))

# 2. Initialize the plot (initial state of the animation)
def init():
    line.set_ydata(np.sin(x))
    return line,

# 3. Define the update function for each frame
def update(frame):
    line.set_ydata(np.sin(x + frame / 10.0))  # Shift the sine wave
    return line,

# 4. Create the animation object
ani = FuncAnimation(fig, update, frames=np.arange(0, 100), init_func=init, blit=True)

# 5. Show the animation
plt.show()
