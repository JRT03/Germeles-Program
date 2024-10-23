import numpy as np
import matplotlib.pyplot as plt

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
    
'''
# parameters of setup
Ho = 2.03 # inital height of fluid in tank
Hf = 2.7 # final height in tank
rho_o = 1.0582 # density of tank content before addition of cargo
rho_c = 0.9979 # density of tank cargo
Q_vol = 0.0718 # volumetric rate of addtion of cargo
tankD = 4 # diamter of cylindrical tank
d = 0.4 # diamter of nozzle
alpha = 0.082 # plume entrainment coefficient 
lamb = 1.16 # plume density distribution coefficient


# calculating the non dimensional constants
beta = 4*np.sqrt(2)*alpha*(Hf/d)
gamma = ((np.pi**2)*23173.22835*(1+lamb**2)*(d**4)*Hf*(rho_o - rho_c))/(32*rho_o*(Q_vol**2))
h = 1- (Ho / Hf)

'''

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


spatial_step_size = eta[1] - eta[0]

# first integration based on beginning f
qm_vector = np.array([q[0],m[0]])
qm_vector_next = RK4(qm_function,f[0],qm_vector,spatial_step_size)
q[1] = qm_vector_next[0]
m[1] = qm_vector_next[1]


for i in range(1,spatial_step_num):

    #using newly calculated q to find f through next point in space
    f[i] = f[i-1] - q[i]*(delta[i] - delta[i-1])

    spatial_step_size = eta[i] - eta[i-1]

    # calculating q and m at the next point in space using values at previous point (RK4) 
    qm_vector = np.array([q[i],m[i]])
    qm_vector_next = RK4(qm_function,f[i],qm_vector,time_step)
    q[i+1] = qm_vector_next[0]
    m[i+1] = qm_vector_next[1]


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

plt.stairs(delta,eta)
plt.xlim(0,1)
plt.show()