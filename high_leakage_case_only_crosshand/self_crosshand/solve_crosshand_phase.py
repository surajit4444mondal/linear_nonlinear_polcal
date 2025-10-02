import numpy as np
from casatools import table
import matplotlib.pyplot as plt
from scipy.optimize import minimize

def get_data(msname):
    tb=table()
    tb.open(msname)
    try:
        data=tb.getcol('DATA')
        flag=tb.getcol('FLAG')
    finally:
        tb.close()
    
    pos=np.where(flag==True)
    data[pos]=np.nan
    
    return data

def get_leak_corrected_data(msname,leakage_file):
    
    leakage=np.loadtxt(leakage_file)
    
    leakage_XX=leakage[:,0]+leakage[:,1]*1j
    leakage_XY=leakage[:,2]+1j*leakage[:,3]
    leakage_YX=leakage[:,4]+1j*leakage[:,5]
    leakage_YY=leakage[:,6]+1j*leakage[:,7]
    
    flagged_ant=np.array([3,19,20,29])
    leakage_XX[flagged_ant]=np.nan
    leakage_XY[flagged_ant]=np.nan
    leakage_YX[flagged_ant]=np.nan
    leakage_YY[flagged_ant]=np.nan
    
    
    data=get_data(msname)
    corrected_data=np.zeros_like(data)
    
    shape=data.shape
    
    num_ant=30

    k=0
    for i in range(num_ant):
        
        leak_matrix1=np.array([[leakage_XX[i],leakage_XY[i]],[leakage_YX[i],leakage_YY[i]]],dtype=complex)
        inv_leak1=np.linalg.inv(leak_matrix1)
        for j in range(i+1,num_ant):
            leak_matrix2=np.array([[leakage_XX[j],leakage_XY[j]],[leakage_YX[j],leakage_YY[j]]],dtype=complex)
            leak_matrix2_adj=np.conj(leak_matrix2).T
            inv_leak2=np.linalg.inv(leak_matrix2_adj)
            
            vis_mat=np.array([[data[0,0,k],data[1,0,k]],[data[2,0,k],data[3,0,k]]],dtype=complex)
            
            corrected_data[:,0,k]=np.matmul(np.matmul(inv_leak1,vis_mat),inv_leak2).flatten()
            k+=1
    return corrected_data

def compute_IQUV(vis_data):
    avg_data=np.nanmean(vis_data,axis=(1,2))
    Idata=0.5*(avg_data[0]+avg_data[3]).real
    Vdata=0.5*(avg_data[0]-avg_data[3]).real
    Qdata=0.5*(avg_data[1]+avg_data[2]).real
    Udata=0.5*(avg_data[1]-avg_data[2]).imag
    
    return [Idata,Qdata,Udata,Vdata]
    
def rotate_QUV(params,observed_stokes,models,return_corrected=False):
   
    num_sources=len(observed_stokes)
    sources_RL_basis=[None]*num_sources

    
    theta=params[0]
    epsilon=0#params[1]
    phi=0#params[2]
    
    theta_matrix=np.array([[np.exp(1j*theta),0],\
                            [0,1]],dtype=complex)

    epsilon_matrix=np.array([[np.cos(epsilon),1j*np.sin(epsilon)],\
                            [1j*np.sin(epsilon), np.cos(epsilon)]],dtype=complex)
    

    phi_matrix=np.array([[np.cos(phi),np.sin(phi)],\
                        [-np.sin(phi),np.cos(phi)]],dtype=complex)
                        
    jones_matrix=np.matmul(theta_matrix,\
                            np.matmul(epsilon_matrix,phi_matrix))
    
    jones_matrix_inv=np.linalg.inv(jones_matrix)
    
    jones_matrix_adj=np.conj(jones_matrix).T
    
    jones_matrix_adj_inv=np.linalg.inv(jones_matrix_adj)

                    

    for j,source_stokes in enumerate(observed_stokes):

        source_RL=np.array([[source_stokes[0]+source_stokes[3],\
                            source_stokes[1]+1j*source_stokes[2]],\
                            [source_stokes[1]-1j*source_stokes[2],\
                            source_stokes[0]-source_stokes[3]]],dtype=complex)
        
        
        sources_RL_basis[j]=np.matmul(np.matmul(jones_matrix_inv,source_RL),jones_matrix_adj_inv)
        


    
    IQUV_corrected=np.zeros((num_sources,4))
    sum1=0
    
    #print ("==================================================")

    for j,(model,corrected) in enumerate(zip(models,sources_RL_basis)):
        #print (corrected)
        IQUV_corrected[j,:]=compute_IQUV(np.expand_dims(np.array(corrected).flatten(),(1,2)))     
          
        sum1+=np.sum(IQUV_corrected[j,:]-np.array(model)[:])**2
        
        #print (IQUV_corrected[j,:])
    #print ("==================================================")
    #print (sum1)

    if return_corrected:
        return IQUV_corrected
    
    return (sum1)

def solve_unitary_matrices(models,observed_stokes):
    res1=minimize(rotate_QUV,[0],args=(observed_stokes,models,\
                    ),method='Nelder-Mead',\
                    bounds=[[-3.14159,3.14159]])
    if res1.success:
        thetas=res1.x
    else:
        thetas=[np.nan,np.nan,np.nan]
    return thetas
            
msnames=['polarized_source_I50_Q1_U5_V0.ms','polarized_source_I50_Q5_U1_V0.ms']

models=[[50,1,5,0],[50,5,1,0]]

IQUV_corrected=[[],[]]

for j,msname in enumerate(msnames):
    corrected_vis=get_leak_corrected_data(msname,'quartical_ant_leaks.txt')
    IQUV_corrected[j]=compute_IQUV(corrected_vis)

#print (IQUV_corrected)    


thetas=solve_unitary_matrices(models,IQUV_corrected)

print (np.array(thetas)*180/np.pi)
#thetas=[53*np.pi/180,0.1,0.1]

corrected=rotate_QUV(thetas,IQUV_corrected,models,return_corrected=True)

print (corrected)

            
            
            
