from class_definitions import *

# Setting up the plume
plume1 = plume(1,0)
plume_list = [plume1]

# Setting up the time
time_step = 0.01
end_time = 20.05
time_array = np.arange(0,end_time,time_step)
plotting_times = [0.1,1,3,10,20]

# Setting up the box
# The reason the box needs the time array as input is that it needs to know the length of the vectors it is going to make
box1 = box(inital_delta=np.array([]), inital_zeta_steps=np.array([]),time_array=time_array,plume_list=plume_list)

# Iterating through time

box1.solve_plume()
for time in time_array:
    
    for plot_time in plotting_times:
        if np.isclose(time,plot_time):
            box1.plot_stratification(plot_time)
            box1.plot_volume_flux(0,time)
            box1.plot_bouyancy_flux(0,time)
    

    box1.deposit_layer(time_step)
    box1.solve_plume()

plt.legend()
plt.show()


# TO DO 
# edit the class definitions to be able to plot stuff (include analyitic bool to also plot analytic solution)
# edit the plume deposition so they deposit at the right level 
# validate the program foruse with starting stratification
# implement the ventilation flow


# FOR SLIDES
# think about how to implement multiple plumes
# make draft slides for the project