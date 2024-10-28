import numpy as np 
import matplotlib.pyplot as plt

def RK4(func,f_val,qm4_init,step):
    k1 = func(f_val,qm4_init)
    k2 = func(f_val,qm4_init + k1*(step/2))
    k3 = func(f_val,qm4_init + k2*(step/2))
    k4 = func(f_val, qm4_init + k3*step)

    qm4_next = qm4_init + step * (k1 + 2*k2 + 2*k3 + k4) / 6

    return qm4_next

def qm4_funciton(f_val,qm4_vector):
    # creating output vector (derivative of input vector)
    qm4_dot = np.zeros(2)
    # updating the output vecotr using the equations provided
    qm4_dot[0] = (qm4_vector[1])**(1/4)
    qm4_dot[1] = 2*f_val*qm4_vector[0]

    return qm4_dot
 

initial_spatial_steps = 100
spatial_steps = initial_spatial_steps

max_time = 20
time_step = 0.1

time = np.arange(0,max_time,time_step)


eta = np.zeros(spatial_steps + len(time)-1) # there is one less time step than there are entries in time array
eta[:spatial_steps+1] = np.linspace(0,1,spatial_steps+1)

# delta: the density profile within the box
# f: bouyancy flux of the plume
# q: volume flux of the plume
# m^2: momentum flux of the plume
delta = np.zeros(len(eta) - 1) # delta is a stepped profile so there is one less value than coordinate positions
f = np.zeros(len(eta) - 1) # f is also a stepped profile
q = np.zeros(len(eta))
m = np.zeros(len(eta))


# Initially the similarity solution for free space solution is applied
f[:spatial_steps] = 1
q[:spatial_steps+1] = (3/10)*((18/5)**(1/3))*eta[:spatial_steps+1]**(5/3)
m[:spatial_steps+1] = (1/2)*((18/5)**(1/3))*eta[:spatial_steps+1]**(2/3)


# First depositon of plume onto top layer
eta[:spatial_steps+1] -= q[:spatial_steps+1]*time_step # The layers are squashed down by q * timestep
front_index = spatial_steps # the index of the front is static and does not change from this point
spatial_steps += 1
eta[spatial_steps] = 1 # Final eta postiion should still be 1 at end of the box
delta[spatial_steps-1] = f[spatial_steps-2]/q[spatial_steps-1] 


for t in range(1,len(time)):
    # Plume fluxes at entrance 
    f[front_index-1] = 1 # This is the step value just before the front
    q[front_index] = (3/10)*((18/5)**(1/3))*eta[front_index]**(5/3)
    m[front_index] = (1/2)*((18/5)**(1/3))*eta[front_index]**(2/3)

    for i in range(front_index,spatial_steps-1):
        # Application of update rule for f based on analytical integration
        f[i] = f[i-1] - q[i]*(delta[i] - delta[i-1])
        
        # Perform RK4 integration to find next value of q and m
        qm4_vector = np.array([q[i],m[i]**4])
        space_step = eta[i+1] - eta[i]
        qm4_next = RK4(qm4_funciton,f[i],qm4_vector,space_step)
        q[i+1] = qm4_next[0]
        m[i+1] = qm4_next[1]**(1/4)

    eta[:spatial_steps+1] -= q[:spatial_steps+1]*time_step
    spatial_steps += 1
    eta[spatial_steps] = 1
    delta[spatial_steps-1] = delta[spatial_steps-2] + f[spatial_steps-2]/q[spatial_steps-1] 
    


            

    


        






