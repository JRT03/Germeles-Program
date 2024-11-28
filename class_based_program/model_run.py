from class_definitions import *

# Setting up the plume
plume1 = plume(1,0)
plume_list = [plume1]

# Setting up the time
time_step = 0.01
end_time = 25.1
time_array = np.arange(0,end_time,time_step)
plotting_times = [0.5,5,25]
# [0.5,5,25]

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

#box1.figd = plt.figure()
#box1.axd = box1.figd.add_subplot()
#box1.axd.grid()
#box1.axd.plot(2*np.linspace(2,3,3),np.linspace(2,3,3),color='black',linestyle='--',label='Analytic')
box1.axd.set_xlim(0,8)
box1.axd.set_ylim(0,1)
box1.axd.set_title("Stratification Within Room",fontsize=14)
box1.axd.set_xlabel("$\\delta$",fontsize=14)
box1.axd.set_ylabel("$\\zeta$",fontsize=14)
#box1.axd.legend(loc='lower right')
box1.figd.savefig('figures//tau3_ventilation',dpi=300)
plt.show()