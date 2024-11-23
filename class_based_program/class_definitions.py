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
        

    def print_params(self):
        print("bouyancy flux:", self.bouyant_flux)
        print("plume height:",self.height)


class box:
    def __init__(self,time_step,time_length,inital_delta,inital_zeta_steps,plume_list):
        
        # Handling the time array
        try:
            self.time_array = np.arange(0,time_length,time_step)
        except:
            print('You have invalid combination of time length and time step')
        
        # Handling the intial stratification delta values
        if not isinstance(inital_delta, (np.ndarray,list)):
            raise Exception('Inital delta should be an array or list')
        self.delta_array = np.zeros(len(inital_delta)+len(self.time_array)+1)
        self.delta_array[0:len(inital_delta)] = inital_delta[0:len(inital_delta)]

        # Handling the intial stratification zeta steps
        if not isinstance(inital_zeta_steps, (np.ndarray,list)):
            raise Exception('Inital delta should be an array or list')
        
        if len(inital_delta) != (len(inital_zeta_steps) - 1):
            raise Exception('The length of intial_zeta_steps should be one more than length of inital_delta')
        self.zeta_steps_array = np.zeros(len(inital_zeta_steps))
        self.zeta_steps_array[0:len(inital_zeta_steps)] = inital_zeta_steps[0:len(inital_zeta_steps)]


        self.plume_list = plume_list
        for plume_obj in self.plume_list:
            if not isinstance(plume_obj, plume):
                raise Exception('You have a non-plume object in your plume list')

        

    def printer(self):
        for i in range(len(self.plume_list)):    
            self.plume_list[i].print_params()
        