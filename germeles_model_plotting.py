import numpy as np
import matplotlib.pyplot as plt
plt.style.use('dark_background')

# function required for perfoming runge kutta integration
def RK4(func,f_val,y_init,t_step):
    k1 = func(f_val,y_init)
    k2 = func(f_val,y_init + k1*(t_step/2))
    k3 = func(f_val,y_init + k2*(t_step/2))
    k4 = func(f_val, y_init + k3*t_step)

    y_next = y_init + t_step * (k1 + 2*k2 + 2*k3 + k4) / 6

    return y_next

def qm_function(f_val,qm_vector):
    qm_dot = np.zeros(2)
    qm_dot[0] = beta * np.sqrt(qm_vector[1])
    qm_dot[1] = gamma * (qm_vector[0] * f_val)/(qm_vector[1])

    return qm_dot


def find_time_index(values,array):
    indices = []
    for i in range(len(values)):
        j = 0
        found = False
        while j < len(array) and not found:
            if float(array[j]) >= values[i]:
                indices.append(j)
                found = True
            j += 1
        if not found:
            indices.append('n/a')

    return indices

# begin figure
fig, axes = plt.subplots(2, 2, figsize=(10, 8)) 
axes = axes.flatten()
print(type(axes[0]))

beta = 3.131
gamma = 1.440
h = 0.248

# number of spatial steps within the domain
spatial_step_num = 100
# set the number of time steps that will be simulated
time_step_num = 100

# setting up the spatial steps vector
eta = np.zeros(spatial_step_num + 1 + time_step_num) # there is one more value of eta than there are spatial steps
eta[:spatial_step_num+1] = np.linspace(0,1-h,spatial_step_num+1)

#setting up the time vector 
time_vector = np.linspace(0,1,time_step_num+1) # there is one more time value than there are time steps (moves between values)
time_step = time_vector[1]

# finds the indexes of the times that we want to plot
plot_times = [0.2,0.4,0.8,1]
plot_indexes = find_time_index(plot_times,time_vector)

# specifying inital stratification
delta = np.zeros(spatial_step_num + time_step_num) 
delta[:spatial_step_num] = 0


# setting up other variables including applications of initial conditions
f = np.zeros(spatial_step_num + time_step_num) # there is one less value of f than there is etas 
f[0] = 1 - delta[0]
m = np.zeros(spatial_step_num + 1 + time_step_num) # there are the same number of ms as etas 
m[0] = 1
q = np.zeros(spatial_step_num + 1 + time_step_num) # there are the same number of qs as etas
q[0] = 1


# This is the first spatial step size (could replace this with a stenciling operation)
spatial_step_size = eta[1] - eta[0]

# first integration based on beginning f
qm_vector = np.array([q[0],m[0]])
qm_vector_next = RK4(qm_function,f[0],qm_vector,spatial_step_size)
q[1] = qm_vector_next[0]
m[1] = qm_vector_next[1]

# loop to solve the rest of the proflie at zero time
for i in range(1,spatial_step_num):

    #using newly calculated q to find f through next point in space
    f[i] = f[i-1] - q[i]*(delta[i] - delta[i-1])

    spatial_step_size = eta[i] - eta[i-1]

    # calculating q and m at the next point in space using values at previous point (RK4) 
    qm_vector = np.array([q[i],m[i]])
    qm_vector_next = RK4(qm_function,f[i],qm_vector,time_step)
    q[i+1] = qm_vector_next[0]
    m[i+1] = qm_vector_next[1]

# time step loop
for j in range(time_step_num):
    
    # squeezing down of exisitng spatial steps
    eta[:spatial_step_num+1] = eta[:spatial_step_num+1] - h*(q[:spatial_step_num+1] - 1)*time_step


    # specifying delta value of new spatial step added
    delta[spatial_step_num] = delta[spatial_step_num-1] + f[spatial_step_num-1]/q[spatial_step_num]

    # specifying lcoation of new spatial step
    eta[spatial_step_num + 1] = 1 - h*(1-time_vector[j]-time_step)

    # updating the number of spatial steps in domain
    spatial_step_num += 1

    # first integration based on beginning f
    spatial_step_size = eta[1] - eta[0]
    qm_vector = np.array([q[0],m[0]])
    qm_vector_next = RK4(qm_function,f[0],qm_vector,spatial_step_size)
    q[1] = qm_vector_next[0]
    m[1] = qm_vector_next[1]
    
    
    for i in range(1,spatial_step_num):
        
        #using newly calculated q to find f through next point in space
        f[i] = f[i-1] - q[i]*(delta[i] - delta[i-1])
        
        # calculating q and m at the next point in space using values at previous point (RK4) 
        spatial_step_size = eta[i] - eta[i-1]
        qm_vector = np.array([q[i],m[i]])
        qm_vector_next = RK4(qm_function,f[i],qm_vector,spatial_step_size)
        q[i+1] = qm_vector_next[0]
        m[i+1] = qm_vector_next[1]

    # if the time step is within the time steps listed to be plotted then plot the graph
    if j + 1 in plot_indexes:
       plot_position = plot_indexes.index(j+1)
       axes[plot_position].stairs(delta[:spatial_step_num],eta[:spatial_step_num+1],orientation='horizontal',baseline=None)
       axes[plot_position].set_ylim(0,1)
       axes[plot_position].set_xlim(0,0.45)
       axes[plot_position].set_xlabel('$\delta$')
       axes[plot_position].set_ylabel('$\eta$')
       axes[plot_position].set_title('$\\tau$ = ' + str(plot_times[plot_position]))

plt.tight_layout()
plt.show()