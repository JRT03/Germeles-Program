import numpy as np
import matplotlib.pyplot as plt

# Right now the difference I can see in the code is that I have used the germeles approach with simple staircase profile for f
# Other code uses halfsteps in zeta in the middle of the RK4 to find the next f

def RK4(func,f_val,qm4_init,step):
    '''
    INPUT
    func: function which returns first derivative of q and m^4 given q and m^4s
    f_val: argument required for the function defined in func
    qm4_init: initial values of q and m^4 in a vector
    step: value over which the quantities were being integrated
    OUTPUT
    qm4_next: the next value of q and m^4 returned in a vectore
    '''
    
    k1 = func(f_val,qm4_init)
    k2 = func(f_val,qm4_init + k1*(step/2))
    k3 = func(f_val,qm4_init + k2*(step/2))
    k4 = func(f_val, qm4_init + k3*step)

    qm4_next = qm4_init + step * (k1 + 2*k2 + 2*k3 + k4) / 6

    return qm4_next

def qm4_funciton(f_val,qm4_vector):
    '''
    INPUT
    f_val: value of f at this point in space required for caclculation of the derivative
    qm4_vector: vector containing the value of q and m^4
    OUTPUT
    qm4_dot: derivatives of q and m^4 
    '''
    
    # creating output vector (derivative of input vector)
    qm4_dot = np.zeros(2)
    # updating the output vecotr using the equations provided
    qm4_dot[0] = (qm4_vector[1])**(1/4)
    qm4_dot[1] = 2*f_val*qm4_vector[0]

    return qm4_dot

def free_q(zeta):
    return (3/10)*((18/5)**(1/3))*(zeta)**(5/3)

def free_m(zeta):
    return (1/2)*((18/5)**(1/3))*(zeta)**(2/3)


# This code is initially programmed as if there is no exisitng stratification
delta = [0] # delta will hold the step values of non-dimensional bouyancy 
f = [1]
q = [free_q(1)]
m = [free_m(1)]
zeta_steps = [0,1] # zeta_steps will hold the respective locations of the steps
d_tau = 0.01
tau = np.arange(0,12,d_tau)
print(tau)

plotting_tau = [0.25,1,2,4]

# PLOTTING INITIAL CONDITINOS
zeta_init = np.linspace(0,1,100)
q_init = free_q(zeta_init)
m_init = free_m(zeta_init)

fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(m_init, zeta_init, label="tau = " + str(tau[0]))
ax1.set_title("Plot of m(zeta)")
ax1.set_xlabel("m")
ax1.set_ylabel("zeta")
ax1.legend()

fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
ax2.plot(q_init, zeta_init, label="tau = " + str(tau[0]), color="orange")
ax2.set_title("Plot of q(zeta)")
ax2.set_xlabel("q")
ax2.set_ylabel("zeta")
# END OF PLOTTING INITAL CONDITIONS


# FIRST LAYER OF PLUME IS DEPOSITED
delta.append(0 + 1/free_q(1)) # Append the bouyancy at this stage
zeta_steps[-1] = 1 - (free_q(1)) * d_tau # This is the location of the first step
zeta_steps.append(1) 
q_i = free_q(zeta_steps[1])
m4_i = free_m(zeta_steps[1])**4
f_i = 1 - q_i*(delta[1] - delta[0])
qm4_vector = np.array([q[0],m4_i])
step = zeta_steps[-1] - zeta_steps[1]
qm4_vector_next = RK4(qm4_funciton,1,qm4_vector,step)

f = [1,f_i]
q = [q_i,qm4_vector_next[0]]
m = [m4_i**(1/4),qm4_vector_next[1]**(1/4)]


fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.stairs(delta,zeta_steps, orientation='horizontal',baseline=None,label='Tau =' +str(tau[1]))
ax3.set_title('Stratification after first plume deposition')
ax3.set_xlabel("delta")
ax3.set_ylabel("zeta")


# NOW FIRST LAYER IS DEPOSITED WE CAN ITERATE THROUGH UNTIL END
for i in range(1,len(tau)):
    # NEW LAYER DEPOSITED
    delta.append(delta[-1] + f[-1]/q[-1]) # add the new delta value
    for l in range(1,len(zeta_steps)):
        zeta_steps[l] = zeta_steps[l] - q[l-1]*d_tau # update the step positionos
    zeta_steps.append(1)

    # NOW PLUME CAN BE SOLVED AGAIN
    q = [free_q(zeta_steps[1])]
    m = [free_m(zeta_steps[1])]
    f = [1]

    for j in range(1,len(zeta_steps)-1):
        
        
        f.append(f[-1] - q[-1]*(delta[j] - delta[j-1]))
        qm4_vector = np.array([q[-1],m[-1]**(1/4)])
        step = zeta_steps[j+1] - zeta_steps[j]
        qm4_vector_next = RK4(qm4_funciton,f[-1],qm4_vector,step)
        q.append(qm4_vector_next[0])
        m.append((qm4_vector_next[1])**(1/4))

    for k in range(len(plotting_tau)):
        if np.isclose(plotting_tau[k],tau[i]+d_tau,1e-6):
            print('im_runnin')
            print(tau[i]+d_tau)
            zeta_coord = np.zeros(100 + len(zeta_steps)-2)
            zeta_coord[0:100] = np.linspace(0,zeta_steps[1],100)
            m_array = np.zeros(len(zeta_coord))
            q_array = np.zeros(len(zeta_coord))
            m_array[0:100] = free_m(zeta_coord[0:100])
            q_array[0:100] = free_q(zeta_coord[0:100])

            for o in range(2,len(zeta_steps)):
                zeta_coord[99 + o -1] = zeta_steps[o]
                m_array[99 + o -1] = m[o-1]
                q_array[99 + o -1] = q[o-1]
            
            # analytic bit
            
            front_pos = ((1 + 0.2*((18/5)**(1/3)))*(tau[i]+d_tau))**(-3/2)
            front_delta = (10/3)*((5/18)**(1/3))
            f_hat = (1 - front_pos**(5/3))/(1 - front_pos)
            delta_inf = np.zeros(len(zeta_coord))
            delta_inf[100:] = 5*((5/18)**(1/3))*(zeta_coord[100:]**(-2/3))*(1 - zeta_coord[100:]*10/39 - (zeta_coord[100:]**2)*155/8112)
            analytical_const = 5*((5/18)**(1/3))*(((front_pos**(-2/4) - 1)/(1 - front_pos)) + 3*f_hat*((1-front_pos**(1/3))/(1- front_pos) - (5/8)*(1 - front_pos**(4/3))/(1 - front_pos) - (155/56784)*(1-front_pos**(7/3))/(1-front_pos)))
            delta_analytic = np.zeros(len(zeta_coord))
            delta_analytic[100:] = (f_hat**(2/3))*delta_inf[100:] - analytical_const


            #print(m)

            ax1.plot(m_array, zeta_coord, label="tau = " + str(plotting_tau[k]) )
            ax2.plot(q_array, zeta_coord, label="tau = " + str(plotting_tau[k]))
            ax3.stairs(delta,zeta_steps, orientation='horizontal',baseline=None,label='Tau =' +str(plotting_tau[k]))
            ax3.plot(delta_analytic,zeta_coord)

            

ax3.stairs(delta,zeta_steps, orientation='horizontal',baseline=None,label='Tau =' +str(tau[-1]+d_tau))

ax1.legend()
ax2.legend()
ax3.legend()
plt.show()