import numpy as np
import matplotlib.pyplot as plt

class plume:
    def __init__(self,bouyant_flux,height):
        # Hndling the bouancy flux input
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
    
    
    def print_full_state(self):
        print('f array: ',self.f_array)
        print('m array: ',self.m_array)
        print('q array: ',self.q_array)


class box:
    def __init__(self,inital_delta,inital_zeta_steps,time_array,plume_list):
    
        
        # Handling the intial stratification delta values
        if not isinstance(inital_delta, np.ndarray):
            raise Exception('Inital delta should be an array or list')
        self.delta_array = np.zeros(len(inital_delta)+len(time_array)+1)
        self.delta_array[0:len(inital_delta)] = inital_delta[0:len(inital_delta)]

        # Handling the intial stratification zeta steps
        if not isinstance(inital_zeta_steps, np.ndarray):
            raise Exception('Inital delta should be an array or list')
        
        self.zeta_steps_array = np.zeros(len(inital_zeta_steps) + len(time_array) + 2)
        self.zeta_steps_array[0:len(inital_zeta_steps)] = inital_zeta_steps[0:len(inital_zeta_steps)]
        self.zeta_steps_array[len(inital_zeta_steps)+1] = 1

        self.step_num = len(inital_zeta_steps)

        self.plume_list = plume_list
        for plume_obj in self.plume_list:
            if not isinstance(plume_obj, plume):
                raise Exception('You have a non-plume object in your plume list')

            plume_obj.f_array = np.zeros(len(inital_delta) + len(time_array) + 1)
            plume_obj.m_array = np.zeros(len(inital_delta) + len(time_array) + 1)
            plume_obj.q_array = np.zeros(len(inital_delta) + len(time_array) + 1)
   
    
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
            self.delta_array[self.step_num+1] = self.delta_array[self.step_num] + plume_obj.f_array[self.step_num]/plume_obj.q_array[self.step_num]
            self.zeta_steps_array[1:self.step_num+2] = self.zeta_steps_array[1:self.step_num+2] - plume_obj.q_array[0:self.step_num+1]*time_step
            self.zeta_steps_array[self.step_num+2] = 1
            self.step_num +=1


    def plume_printer(self):
        for plume_obj in self.plume_list:
            plume_obj.print_full_state()

    def strat_printer(self):
        print('delta array:',self.delta_array)
        print('zeta steps:',self.zeta_steps_array)
        print('number of steps:',self.step_num)
        