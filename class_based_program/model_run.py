from class_definitions import *

# Setting up the plume
plume1 = plume(1,0)
plume_list = [plume1]

# Setting up the time
time_step = 0.01
end_time = 0.1
time_array = np.arange(0,end_time,time_step)

# Setting up the box
box1 = box(inital_delta=np.array([]), inital_zeta_steps=np.array([]),time_array=time_array,plume_list=plume_list)

# Iterating through time
for ti in range(len(time_array)):
    box1.solve_plume()
    box1.deposit_layer(time_step)


box1.strat_printer()


# TO DO 
# edit the class definitions to be able to plot stuff (include analyitic bool to also plot analytic solution)
# edit the plume deposition so they deposit at the right level 
# validate the program foruse with starting stratification
# implement the ventilation flow


# FOR SLIDES
# think about how to implement multiple plumes
# make draft slides for the project