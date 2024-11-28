from class_definitions import *

# Setting up the plume
plume1 = plume(1,0)
plume_list = [plume1]

# Setting up the time
time_step = 0.01
end_time = 25.01
time_array = np.arange(0,end_time,time_step)
plotting_times = [0.62,1.25,2.5,6.25,25]

# Setting up the box
# The reason the box needs the time array as input is that it needs to know the length of the vectors it is going to make
box1 = box(plume_list,time_array,a_star_vent=0.2)


# Initally solving the plume
box1.solve_plume()

# Iteration through time
for time in time_array:

    # Plotting functionality
    for plot_time in plotting_times:

        if np.isclose(time,plot_time):
            print('Currently plotting at time:', plot_time)
            box1.plot_stratification(plot_time,analytic=False)
            box1.plot_volume_flux(0,time)
            box1.plot_bouyancy_flux(0,time)
            box1.plot_momentum_flux(0,time)
    
    # main solution steps
    box1.deposit_layer(time_step)
    box1.solve_plume()

plt.legend()
plt.show()