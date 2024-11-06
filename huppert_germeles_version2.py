import numpy as np
import matplotlib.pyplot as plt

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
tau = 0
d_tau = 0.01

# PLOTTING INITIAL CONDITINOS

zeta_init = np.linspace(0,1,100)
q_init = free_q(zeta_init)
m_init = free_m(zeta_init)

fig1 = plt.figure()
ax1 = fig1.add_subplot(1, 1, 1)
ax1.plot(m_init, zeta_init, label="m(zeta)")
ax1.set_title("Plot of m(zeta)")
ax1.set_xlabel("m")
ax1.set_ylabel("zeta")
ax1.legend()

fig2 = plt.figure()
ax2 = fig2.add_subplot(1, 1, 1)
ax2.plot(q_init, zeta_init, label="q(zeta)", color="orange")
ax2.set_title("Plot of q(zeta)")
ax2.set_xlabel("q")
ax2.set_ylabel("zeta")
ax2.legend()

# END OF PLOTTING INITAL CONDITIONS



# FIRST LAYER OF PLUME IS DEPOSITED
delta.append(0 + 1/free_q(1)) # Append the bouyancy at this stage
zeta_steps[-1] = 1 - (free_q(1)) * d_tau # This is the location of the first step
zeta_steps.append(1) 
tau += d_tau 
q_i = free_q(zeta_steps[1])
m4_i = free_m(zeta_steps[1])**4
f_i = 1 - q_i*(delta[1] - delta[0])
qm4_vector = np.array([q[0],m4_i])
step = zeta_steps[-1] - zeta_steps[1]
qm4_vector_next = RK4(qm4_funciton,1,qm4_vector,step)

print("f at top of domain: ",f_i)
print("q at top of domain: ",qm4_vector_next[0])
print("m at top of domain: ",qm4_vector_next[1]**(1/4))

f = [1,f_i]
q = [q_i,qm4_vector_next[0]]
m = [m4_i**(1/4),qm4_vector_next[1]**(1/4)]
print(f)
print(q)
print(m)


fig3 = plt.figure()
ax3 = fig3.add_subplot(1,1,1)
ax3.stairs(delta,zeta_steps, orientation='horizontal',baseline=None)
ax3.set_title('Stratification after first plume deposition')
ax3.set_xlabel("delta")
ax3.set_ylabel("zeta")


#zeta_second = np.zeros(101)
#zeta_second[0:100] = np.linspace(0,zeta_steps[1],100)
#zeta_second[100] = 1

#q_second = np.zeros(101)
#q_second[0:100] = free_q(zeta_second[0:100])
#q_second[100] = qm4_vector_next[0]

#ax2.plot(q_second,zeta_second)


# NOW FIRST LAYER IS DEPOSITED WE CAN ITERATE THROUGH UNTIL END
while tau < 2*d_tau:
    # CALCULATE THE PLUME CHARACTERISTICS
    f = [1]
    for i in range(len(delta)):
      pass

    tau += d_tau
    



plt.show()