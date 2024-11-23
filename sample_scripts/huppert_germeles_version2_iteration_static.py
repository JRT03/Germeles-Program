import numpy as np
import matplotlib.pyplot as plt

# Right now the difference I can see in the code is that I have used the germeles approach with simple staircase profile for f
# Other code uses halfsteps in zeta in the middle of the RK4 to find the next f


def zeta_integration(step,f_inital,q_inital,m_inital,delta,prev_delta):
    zeta_halfstep = step / 2
    
    q_halfstep = np.sqrt(q_inital**2 + zeta_halfstep*2*m_inital*q_inital)
    m_halfstep = np.sqrt(m_inital**2 + zeta_halfstep*q_inital*f_inital/(m_inital**2))
    f_halfstep = f_inital - zeta_halfstep*q_inital*(delta - prev_delta)/step

    q_final = np.sqrt(q_inital**2 + step*2*m_halfstep*q_halfstep)
    m_final = np.sqrt(m_inital**2 + step*q_halfstep*f_halfstep/(m_halfstep**2))
    f_final = f_inital - q_halfstep*(delta - prev_delta)

    return q_final, m_final, f_final 


def zeta_integrationRK4(step,f_inital,q_inital,m_inital,delta,prev_delta):
    half_step = step/2

    f_final = f_inital - q_inital * (delta - prev_delta)

    k1q =  m_inital
    k1m =  (q_inital * f_inital) / (2 * m_inital**3)
    

    #f2 = f_inital - (q_inital + k1q*half_step) * (delta - prev_delta) * half_step/step 

    k2q =  m_inital + k1m*half_step
    k2m = (q_inital + k1q * half_step) * (f_final) / (2 * (m_inital + k1m * half_step)**3)

    #f3 = f_inital - (q_inital + k2q * half_step) * (delta - prev_delta) * half_step/step

    k3q =  m_inital + k2m * half_step
    k3m = (q_inital + k2q * half_step) * (f_final) / (2 * (m_inital + k2m * half_step)**3)

    #f4 = f_inital - (q_inital + k3q * step) * (delta - prev_delta)
    #f_final = f_inital - (q_inital + k1q * step) * (delta - prev_delta)
    k4q =  m_inital + k3m*step
    k4m = (q_inital + k3q * step) * (f_final) / (2 * (m_inital + k3m * half_step)**3)


    q_final = q_inital + step*(k1q + 2*k2q + 2*k3q + k4q)/6
    m_final = m_inital + step*(k1m + 2*k2m + 2*k3m + k4m)/6
    return q_final, m_final,f_final


def plot_momentum(zeta_steps_array, m_array,zeta_size,tau):
    global fig1
    global ax1
    if fig1 is None:
        fig1 = plt.figure()
        ax1 = fig1.add_subplot()
        ax1.set_title("Plume Momentum Flux Through Room")
        ax1.set_xlabel("m")
        ax1.set_ylabel("$\\zeta$")
        ax1.grid()

    plotting_zeta = np.zeros(100 + zeta_size-3)
    plotting_zeta[0:100] = np.linspace(0,zeta_steps_array[1],100)
    plotting_zeta[99:] = zeta_steps_array[2:zeta_size]


    plotting_m = np.zeros(len(plotting_zeta))
    plotting_m[0:100] = free_m(plotting_zeta[0:100])
    plotting_m[99:] = m_array[1:zeta_size-1]
    ax1.plot(plotting_m, plotting_zeta , label="tau = " + str(tau))

def plot_bouyancy_flux(zeta_steps_array, q_array,zeta_size,tau):
    global fig2
    global ax2
    if fig2 is None:
        fig2 = plt.figure()
        ax2 = fig2.add_subplot()
        ax2.set_title("Plume Bouyancy Flux Through Room")
        ax2.set_xlabel("q")
        ax2.set_ylabel("$\\zeta$")
        ax2.grid()

    ax2.set_ylabel("zeta")

    plotting_zeta = np.zeros(100 + zeta_size-3)
    plotting_zeta[0:100] = np.linspace(0,zeta_steps_array[1],100)
    plotting_zeta[99:] = zeta_steps_array[2:zeta_size]


    plotting_q = np.zeros(len(plotting_zeta))
    plotting_q[0:100] = free_q(plotting_zeta[0:100])
    plotting_q[99:] = q_array[1:zeta_size-1]
    ax2.plot(plotting_q, plotting_zeta , label="tau = " + str(tau))


def plot_stratification(zeta_steps_array, delta_array, zeta_size, tau):
    global fig3
    global ax3
    if fig3 is None:
        fig3 = plt.figure()
        ax3 = fig3.add_subplot()
        ax3.set_title("Stratification Within Room")
        ax3.set_xlabel("$\\delta$")
        ax3.set_ylabel("$\\zeta$")
        ax3.grid()

    ax3.stairs(delta_array[:zeta_size-1], zeta_steps_array[:zeta_size],orientation='horizontal',baseline=None,label='$\\tau$ = ' + str(tau))




def free_q(zeta):
    return (3/10)*((18/5)**(1/3))*(zeta)**(5/3)

def free_m(zeta):
    return (1/2)*((18/5)**(1/3))*(zeta)**(2/3)


fig1 = None
ax1 = None

fig2 = None
ax2 = None

fig3 = None
ax3 = None

# This code is initially programmed as if there is no exisitng stratification
delta = [0] # delta will hold the step values of non-dimensional bouyancy 
f = [1]
q = [free_q(1)]
m = [free_m(1)]
zeta_steps = [0,1] # zeta_steps will hold the respective locations of the steps
d_tau = 0.01
tau = np.arange(0,12,d_tau)

plotting_tau = [0.25,1,2,4]

delta_array = np.zeros(len(tau)+1)
f_array = np.zeros(len(tau)+1)
q_array = np.zeros(len(tau)+1)
m_array = np.zeros(len(tau)+1)
f_array = np.zeros(len(tau)+1)
zeta_steps_array = np.zeros(len(tau)+2)
zeta_steps_array[-1] = 1
zeta_steps_array[-2] = 0
f_array[0] = 1
q_array[0] = free_q(1)
m_array[0] = free_m(1)


# PLOTTING INITIAL CONDITINOS
zeta_init = np.linspace(0,1,100)
q_init = free_q(zeta_init)
m_init = free_m(zeta_init)
'''
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
'''
# END OF PLOTTING INITAL CONDITIONS


# FIRST LAYER OF PLUME IS DEPOSITED
delta.append(0 + 1/free_q(1)) # Append the bouyancy at this stage
zeta_steps[-1] = 1 - (free_q(1)) * d_tau # This is the location of the first step
zeta_steps.append(1) 
q_i = free_q(zeta_steps[1])
m_i = free_m(zeta_steps[1])


step = zeta_steps[-1] - zeta_steps[1]
q_next, m_next, f_next = zeta_integrationRK4(step,1,q_i,m_i,delta[1],delta[0])


f = [1,f_next]
q = [q_i,q_next]
m = [m_i,m_next]


delta_array[1] = 0 + 1/free_q(1)
zeta_steps_array[1:2] = 1 - free_q(1)*d_tau
zeta_steps_array[2] = 1

f_array[0] = 1
q_array[0] = free_q(zeta_steps_array[1])
m_array[0] = free_m(zeta_steps_array[1])

step_ar = zeta_steps_array[2] - zeta_steps_array[1]
q_array[1],m_array[1],f_array[1] = zeta_integrationRK4(step,1,q_array[0],m_array[0],delta[1],delta[0])
 

fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.stairs(delta,zeta_steps, orientation='horizontal',baseline=None,label='Tau =' +str(tau[1]))
ax3.set_title('Stratification of Space')
ax3.set_xlabel("delta")
ax3.set_ylabel("zeta")


zeta_size = 3

# NOW FIRST LAYER IS DEPOSITED WE CAN ITERATE THROUGH UNTIL END
for i in range(1,len(tau)):
    # NEW LAYER DEPOSITED
    delta.append(delta[-1] + f[-1]/q[-1]) # add the new delta value
    for l in range(1,len(zeta_steps)):
        zeta_steps[l] = zeta_steps[l] - q[l-1]*d_tau # update the step positionos
    zeta_steps.append(1)

    delta_array[i+1] = delta_array[i] + f_array[i]/q_array[i]
    zeta_steps_array[1:i+2] = zeta_steps_array[1:i+2] - q_array[0:i+1]*d_tau
    zeta_steps_array[i+2] = 1
    zeta_size += 1


    # NOW PLUME CAN BE SOLVED AGAIN
    q = [free_q(zeta_steps[1])]
    m = [free_m(zeta_steps[1])]
    f = [1]

    q_array[0] = free_q(zeta_steps_array[1])
    m_array[0] = free_m(zeta_steps_array[1])
    f_array[0] = 1


    for j in range(1,len(zeta_steps)-1):
        
        
        step = zeta_steps[j+1] - zeta_steps[j]
        
        q_next, m_next, f_next = zeta_integrationRK4(step,f[j-1],q[-1],m[j-1],delta[j],delta[j-1])
        
        q.append(q_next)
        m.append(m_next)
        f.append(f_next)

    for j in range(1,i+2): 
        step_ar = zeta_steps[j+1] - zeta_steps[j]
        q_array[j],m_array[j],f_array[j] = zeta_integrationRK4(step_ar,f_array[j-1],q_array[j-1],m_array[j-1],delta_array[j],delta_array[j-1])


    for k in range(len(plotting_tau)):
        if np.isclose(plotting_tau[k],tau[i]+d_tau,1e-6):
            plot_momentum(zeta_steps_array,m_array,zeta_size, tau[i]+d_tau)
            plot_bouyancy_flux(zeta_steps_array,q_array,zeta_size,tau[i] + d_tau)
            plot_stratification(zeta_steps_array,delta_array,zeta_size,tau[i]+d_tau)


    '''
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
        
            front_pos = (1 + 0.2*((18/5)**(1/3))*(tau[i]+d_tau))**(-3/2) # this is zeta_0
            front_delta = (10/3)*((5/18)**(1/3))
            f_hat = (1 - front_pos**(5/3))/(1 - front_pos)
            delta_inf = np.zeros(len(zeta_coord))
            delta_inf[100:] = 5*((5/18)**(1/3))*(zeta_coord[100:]**(-2/3))*(1 - zeta_coord[100:]*10/39 - (zeta_coord[100:]**2)*155/8112)
            analytical_const = ( 
            5*((5/18)**(1/3))*((front_pos**(-2/3) - 1)/(1 - front_pos)
            + 3*(f_hat**(2/3))*((1 - front_pos**(1/3))/(1- front_pos)
            - (5/78)*(1 - front_pos**(4/3))/(1 - front_pos) 
            - (155/56784)*(1 - front_pos**(7/3))/(1- front_pos)))
            )
            
            delta_analytic = np.zeros(len(zeta_coord))
            delta_analytic[100:] = -((f_hat**(2/3))*delta_inf[100:] - analytical_const)
            delta_analytic[99] = (10/3)*(5/18)**(1/3)


            ax1.plot(m_array, zeta_coord, label="tau = " + str(plotting_tau[k]) )
            ax2.plot(q_array, zeta_coord, label="tau = " + str(plotting_tau[k]))
            ax3.stairs(delta,zeta_steps, orientation='horizontal',baseline=None,label='Tau =' +str(plotting_tau[k]))
            ax3.plot(delta_analytic[99:],zeta_coord[99:],linestyle='--',color='black')
        '''


plot_momentum(zeta_steps_array,m_array,zeta_size,tau[-1]+d_tau)

            

ax3.stairs(delta,zeta_steps, orientation='horizontal',baseline=None,label='Tau =' +str(tau[-1]+d_tau))

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

front_pos = (1 + 0.2*((18/5)**(1/3))*(tau[i]+d_tau))**(-3/2) # this is zeta_0
front_delta = (10/3)*((5/18)**(1/3))
f_hat = (1 - front_pos**(5/3))/(1 - front_pos)
delta_inf = np.zeros(len(zeta_coord))
delta_inf[100:] = 5*((5/18)**(1/3))*(zeta_coord[100:]**(-2/3))*(1 - zeta_coord[100:]*10/39 - (zeta_coord[100:]**2)*155/8112)
analytical_const = ( 
5*((5/18)**(1/3))*((front_pos**(-2/3) - 1)/(1 - front_pos)
+ 3*(f_hat**(2/3))*((1 - front_pos**(1/3))/(1- front_pos)
- (5/78)*(1 - front_pos**(4/3))/(1 - front_pos) 
- (155/56784)*(1 - front_pos**(7/3))/(1- front_pos)))
)

delta_analytic = np.zeros(len(zeta_coord))
delta_analytic[100:] = -((f_hat**(2/3))*delta_inf[100:] - analytical_const)
delta_analytic[99] = (10/3)*(5/18)**(1/3)
ax3.plot(delta_analytic[99:],zeta_coord[99:],linestyle='--',color='black')

if ax1 is not None:
    ax1.legend()

if ax2 is not None:
    ax2.legend()

if ax3 is not None:
    ax3.legend()

plt.show()