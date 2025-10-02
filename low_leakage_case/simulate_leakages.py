import numpy as np


num_ant=30
low_amp_leak=0.01
high_amp_leak=0.07

mean_phase_leak=50
sigma_phase_leak=10

leak_amp_XY=np.random.uniform(low_amp_leak,high_amp_leak,num_ant)

leak_phase_XY=np.random.normal(mean_phase_leak,sigma_phase_leak,num_ant)*np.pi/180

leak_amp_YX=np.random.uniform(low_amp_leak,high_amp_leak,num_ant)

leak_phase_YX=np.random.normal(mean_phase_leak,sigma_phase_leak,num_ant)*np.pi/180

np.savetxt("antenna_leakages.txt",np.array([leak_amp_XY,leak_phase_XY,leak_amp_YX,leak_phase_YX]).T)
