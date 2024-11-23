from class_definitions import *


plume1 = plume(1,0)

plume_list = [plume1]


box1 = box(time_step=0.01,time_length=1,inital_delta=np.array([1,2,3,4]), inital_zeta_steps=np.array([1,2,3,4,5]),plume_list=plume_list)
box1.printer()

