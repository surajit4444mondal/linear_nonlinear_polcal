import numpy as np
from casatools import table


def determine_antenna_leak_matrix(num_ant):#,low_amp_leak,high_amp_leak,mean_phase_leak, sigma_phase_leak):
    
    #leak_amp=np.random.uniform(low_amp_leak,high_amp_leak,num_ant)
    #np.abs(np.random.normal(mean_amp_leak,sigma_amp_leak,num_ant))
    #leak_phase=np.random.normal(mean_phase_leak,sigma_phase_leak,num_ant)*np.pi/180
    
    data=np.loadtxt("antenna_leakages.txt")
    
    leak_amp=data[:,0]
    leak_phase=data[:,1]
    
    ant_leaks=leak_amp*np.exp(1j*leak_phase)
    
    #np.savetxt("antenna_leakages.txt",np.array([leak_amp,leak_phase]).T)
    
    ant_leak_matrix=np.zeros((num_ant,2,2),dtype=complex)
    
    for i in range(num_ant):
        ant_leak_matrix[i,:,:]=np.array([[1,ant_leaks[i]],[np.conj(ant_leaks[i]),1]],dtype=complex)
    
    return ant_leak_matrix
    

def get_crosshand_phase_matrix(cross_phase=50):
    cross_phase_rad=cross_phase*np.pi/180
    
    return np.array([[np.exp(1j*cross_phase_rad),0],[0,1]],dtype=complex)
    

def get_ant_jones_matrix(poldistort,polrot):
    return np.matmul(poldistort,polrot)
    
def get_visibilities(num_ant,source_stokes):#,low_amp_leak,\
                  #  high_amp_leak,mean_phase_leak, \
                  #  sigma_phase_leak):
                    
    num_baseline=int(num_ant*(num_ant-1)/2)
    vis=np.zeros((num_baseline,4),dtype=complex)
    
    source_RR=(source_stokes[0]+source_stokes[3])
    source_LL=(source_stokes[0]-source_stokes[3])
    source_RL=(source_stokes[1]+1j*source_stokes[2])
    source_LR=(source_stokes[1]-1j*source_stokes[2])
    
    source_model_matrix=np.array([[source_RR,source_RL],\
                        [source_LR,source_LL]],dtype=complex)

    leak_matrix=determine_antenna_leak_matrix(num_ant)#,low_amp_leak,\
                                              #  high_amp_leak,mean_phase_leak, \
                                              #  sigma_phase_leak)
                                                
    crossphase_matrix=get_crosshand_phase_matrix()
    
    k=0
    for i in range(num_ant):
        for j in range(i+1,num_ant):
            jones_matrix_ant1=get_ant_jones_matrix(leak_matrix[i],crossphase_matrix)
            jones_matrix_ant2=get_ant_jones_matrix(leak_matrix[j],crossphase_matrix)
            
            jones_ant2_adj=jones_matrix_ant2.conj().T
            
            vis_mat=np.matmul(np.matmul(jones_matrix_ant1,source_model_matrix),jones_ant2_adj)
            
            vis[k,:]=vis_mat.flatten()
            
            del vis_mat, jones_matrix_ant1, jones_matrix_ant2, jones_ant2_adj
            
            k+=1
    return vis
    



msname='CASA/polarized_source_I50_Q1_U5_V0.ms'

source_stokes=[50,1,5,0] ### IQUV
vis=get_visibilities(num_ant=30,source_stokes=source_stokes)#,low_amp_leak=0.01,\
                   # high_amp_leak=0.07,mean_phase_leak=50, \
                   # sigma_phase_leak=10)

vis=np.expand_dims(np.swapaxes(vis,0,1),axis=1)

tb=table()
tb.open(msname,nomodify=False)
try:
    data=tb.getcol('DATA')
    tb.putcol('DATA',vis)
    tb.flush()
finally:
    tb.close()    
    
    
    
    
    
    
    

