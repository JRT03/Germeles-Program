import numpy as np
import matplotlib.pyplot as plt

class plume:
    def __init__(self,bouyant_flux,height):
        # Handling the bouancy flux input
        if not isinstance(bouyant_flux, (int,float)):
            raise Exception('Bouyancy flux specified is not a number')
        self.bouyant_flux = bouyant_flux
        
        # Handling the plume bouyancy flux input
        if not isinstance(height, (int,float)):
            raise Exception('Plume height specified is not a number')
        if height != 0:
            raise Exception('Plume height is not currently accounted for in program (set to zero)')
        self.height = height
        
        # Plume arrays are initalised here to be created when box environement defined
        self.f_array = None
        self.q_array = None
        self.m_array = None

        # Figures are initalised here as none so they can be created when plotting is called for
        self.figm = None
        self.figf = None
        self.figq = None


class box:
    def __init__(self,inital_delta,inital_zeta_steps,time_array,plume_list):
    
        
        # Handling the intial stratification delta values
        if not isinstance(inital_delta, np.ndarray):
            raise Exception('Inital delta should be an array or list')
        self.delta_array = np.zeros(len(inital_delta)+len(time_array)+1)
        self.delta_array[1:len(inital_delta)+1] = inital_delta[:]

        # Handling the intial stratification zeta steps
        if not isinstance(inital_zeta_steps, np.ndarray):
            raise Exception('Inital delta should be an array or list')
        
        self.zeta_steps_array = np.zeros(len(inital_zeta_steps) + len(time_array) + 2)
        self.zeta_steps_array[1:len(inital_zeta_steps)+1] = inital_zeta_steps[:]
        self.zeta_steps_array[len(inital_zeta_steps)+1] = 1

        self.step_num = len(inital_zeta_steps)

        self.plume_list = plume_list
        for plume_obj in self.plume_list:
            if not isinstance(plume_obj, plume):
                raise Exception('You have a non-plume object in your plume list')

            plume_obj.f_array = np.zeros(len(inital_delta) + len(time_array) + 1)
            plume_obj.m_array = np.zeros(len(inital_delta) + len(time_array) + 1)
            plume_obj.q_array = np.zeros(len(inital_delta) + len(time_array) + 1)
        
        # Detla figure is initalised here so it can be created when plotting is called
        self.figd = None
   
    
    @staticmethod
    def free_q(zeta):
        return (3/10)*((18/5)**(1/3))*(zeta)**(5/3)
    
    @staticmethod
    def free_m(zeta):
        return (1/2)*((18/5)**(1/3))*(zeta)**(2/3)

    @staticmethod       
    def zeta_integrationRK4(step,f_inital,q_inital,m_inital,delta,prev_delta):
        half_step = step/2

        f_final = f_inital - q_inital * (delta - prev_delta)

        k1q =  m_inital
        k1m =  (q_inital * f_inital) / (2 * m_inital**3)

        k2q =  m_inital + k1m*half_step
        k2m = (q_inital + k1q * half_step) * (f_final) / (2 * (m_inital + k1m * half_step)**3)

        k3q =  m_inital + k2m * half_step
        k3m = (q_inital + k2q * half_step) * (f_final) / (2 * (m_inital + k2m * half_step)**3)

        k4q =  m_inital + k3m*step
        k4m = (q_inital + k3q * step) * (f_final) / (2 * (m_inital + k3m * half_step)**3)

        q_final = q_inital + step*(k1q + 2*k2q + 2*k3q + k4q)/6
        m_final = m_inital + step*(k1m + 2*k2m + 2*k3m + k4m)/6
        return q_final, m_final,f_final
        
    def solve_plume(self):
        for plume_obj in self.plume_list:
            plume_obj.q_array[0] = self.free_q(self.zeta_steps_array[1])
            plume_obj.m_array[0] = self.free_m(self.zeta_steps_array[1])
            plume_obj.f_array[0] = 1


            for i in range(self.step_num):
                step = self.zeta_steps_array[i+2] - self.zeta_steps_array[i+1]
                plume_obj.q_array[i+1],plume_obj.m_array[i+1],plume_obj.f_array[i+1] = self.zeta_integrationRK4(step,plume_obj.f_array[i],plume_obj.q_array[i],plume_obj.m_array[i],self.delta_array[i+1],self.delta_array[i])


    def deposit_layer(self,time_step):
        for plume_obj in self.plume_list:
            f_negative_index = None
            for i in range(self.step_num+1):
                if plume_obj.f_array[i] < 0:
                    f_negative_index = i
                    break
            
            if f_negative_index is None:
                self.delta_array[self.step_num+1] = self.delta_array[self.step_num] + plume_obj.f_array[self.step_num]/plume_obj.q_array[self.step_num]
                self.zeta_steps_array[1:self.step_num+2] = self.zeta_steps_array[1:self.step_num+2] - plume_obj.q_array[0:self.step_num+1]*time_step
                self.zeta_steps_array[self.step_num+2] = 1
                self.step_num +=1
            
            else:
                #print('begin solve')
                #print('f_negative index: ',f_negative_index)
                f_negative_zeta = self.zeta_steps_array[f_negative_index] + (self.zeta_steps_array[f_negative_index+1]- self.zeta_steps_array[f_negative_index])*(plume_obj.f_array[f_negative_index-1])/(plume_obj.f_array[f_negative_index-1] - plume_obj.f_array[f_negative_index])
                #print('f_negative zeta:',f_negative_zeta)
                #print('previous delta:',self.delta_array[:self.step_num+1])
                self.delta_array[f_negative_index+1:self.step_num+2] = self.delta_array[f_negative_index:self.step_num+1]
                self.delta_array[f_negative_index+1] = self.delta_array[f_negative_index]
                #print('new delta:',self.delta_array[:self.step_num+1])
                
                #print('begin alteration of arrays')
                #print('first zeta steps', print(self.zeta_steps_array[:self.step_num+2]))
                self.zeta_steps_array[1:f_negative_index+1] = self.zeta_steps_array[1:f_negative_index+1] - plume_obj.q_array[0:f_negative_index]*time_step
                #print('after squeezing below',print(self.zeta_steps_array[:self.step_num+2]))
                
                
                self.zeta_steps_array[f_negative_index+1:self.step_num+3] = self.zeta_steps_array[f_negative_index:self.step_num+2]
                #print('after zeta shifting along',self.zeta_steps_array[:self.step_num+3])
                self.zeta_steps_array[f_negative_index+1] = f_negative_zeta
                #print('after new zeta added',self.zeta_steps_array[:self.step_num+3])
                self.step_num += 1
                


    def plot_stratification(self,tau,analytic):
        if self.figd is None:
            self.figd = plt.figure()
            self.axd = self.figd.add_subplot()
            self.axd.set_title("Stratification Within Room")
            self.axd.set_xlabel("$\\delta$")
            self.axd.set_ylabel("$\\zeta$")
            self.axd.grid()
        
        if analytic:
            if len(self.plume_list) != 1:
                print("WARNING Analytical solution not intended for use with multiple plumes")

            zeta_coord = np.linspace(self.zeta_steps_array[1],1,100)
            front_pos = (1 + 0.2*((18/5)**(1/3))*(tau))**(-3/2) # this is zeta_0
            f_hat = (1 - front_pos**(5/3))/(1 - front_pos)
            delta_inf = np.zeros(len(zeta_coord))
            delta_inf = 5*((5/18)**(1/3))*(zeta_coord**(-2/3))*(1 - zeta_coord*10/39 - (zeta_coord**2)*155/8112)
            analytical_const = ( 
            5*((5/18)**(1/3))*((front_pos**(-2/3) - 1)/(1 - front_pos)
            + 3*(f_hat**(2/3))*((1 - front_pos**(1/3))/(1- front_pos)
            - (5/78)*(1 - front_pos**(4/3))/(1 - front_pos) 
            - (155/56784)*(1 - front_pos**(7/3))/(1- front_pos)))
            )
            delta_analytic = np.zeros(len(zeta_coord))
            delta_analytic = -((f_hat**(2/3))*delta_inf - analytical_const)
            self.axd.plot(delta_analytic,zeta_coord,linestyle='--',color='black')

        self.axd.stairs(self.delta_array[:self.step_num+1],self.zeta_steps_array[:self.step_num+2],orientation='horizontal',baseline=None,label='$\\tau$ = ' + str(tau))
        self.axd.legend(loc='lower right')

    def plot_volume_flux(self,plume_index,tau):
        plume_obj = self.plume_list[plume_index]
        if plume_obj.figq is None:
            plume_obj.figq = plt.figure()
            plume_obj.axq = plume_obj.figq.add_subplot()
            plume_obj.axq.set_title("Plume Index: " + str(plume_index))
            plume_obj.axq.set_xlabel("$q$")
            plume_obj.axq.set_ylabel("$\\zeta$")
            plume_obj.axq.grid()
    
        plotting_zeta = np.zeros(100 + self.step_num)
        plotting_zeta[:100] = np.linspace(0,self.zeta_steps_array[1],100)
        plotting_zeta[100:] = self.zeta_steps_array[2:self.step_num+2]

        plotting_q = np.zeros(len(plotting_zeta))
        plotting_q[:100] = self.free_q(plotting_zeta[:100])
        plotting_q[100:] = plume_obj.q_array[1:self.step_num+1]
        
        plume_obj.axq.plot(plotting_q,plotting_zeta,label='$\\tau = $' + str(tau))
        plume_obj.axq.legend()

    def plot_momentum_flux(self,plume_index,tau):
        plume_obj = self.plume_list[plume_index]
        if plume_obj.figm is None:
            plume_obj.figm = plt.figure()
            plume_obj.axm = plume_obj.figm.add_subplot()
            plume_obj.axm.set_title("Plume Index: " + str(plume_index))
            plume_obj.axm.set_xlabel("$m$")
            plume_obj.axm.set_ylabel("$\\zeta$")
            plume_obj.axm.grid()
    
        plotting_zeta = np.zeros(100 + self.step_num)
        plotting_zeta[:100] = np.linspace(0,self.zeta_steps_array[1],100)
        plotting_zeta[100:] = self.zeta_steps_array[2:self.step_num+2]

        plotting_m = np.zeros(len(plotting_zeta))
        plotting_m[:100] = self.free_m(plotting_zeta[:100])
        plotting_m[100:] = plume_obj.m_array[1:self.step_num+1]
        
        plume_obj.axm.plot(plotting_m,plotting_zeta,label='$\\tau = $' + str(tau))
        plume_obj.axm.legend()


    
    def plot_bouyancy_flux(self,plume_index,tau):
        plume_obj = self.plume_list[plume_index]
        if plume_obj.figf is None:
            plume_obj.figf = plt.figure()
            plume_obj.axf = plume_obj.figf.add_subplot()
            plume_obj.axf.set_title("Plume Index: " + str(plume_index))
            plume_obj.axf.set_xlabel("$f$")
            plume_obj.axf.set_ylabel("$\\zeta$")
            plume_obj.axf.grid()
    
        plotting_zeta = np.zeros(100 + self.step_num)
        plotting_zeta[:100] = np.linspace(0,self.zeta_steps_array[1],100)
        plotting_zeta[100:] = self.zeta_steps_array[2:self.step_num+2]

        plotting_f = np.zeros(len(plotting_zeta))
        plotting_f[:100] = 1
        plotting_f[100:] = plume_obj.f_array[1:self.step_num+1]
        
        plume_obj.axf.plot(plotting_f,plotting_zeta,label='$\\tau = $' + str(tau))
        plume_obj.axf.legend()


    def plume_printer(self):
        for plume_obj in self.plume_list:
            plume_obj.print_full_state()

    def strat_printer(self):
        print('delta array:',self.delta_array)
        print('zeta steps:',self.zeta_steps_array)
        print('number of steps:',self.step_num)
        