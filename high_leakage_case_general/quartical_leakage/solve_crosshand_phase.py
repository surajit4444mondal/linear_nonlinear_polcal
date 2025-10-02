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

def get_leak_corrected_data(msname,leakage):
    
    data=get_data(msname)
    corrected_data=np.zeros_like(data)
    
    shape=data.shape
    
    num_ant=30

    k=0
    for i in range(num_ant):
        leak_matrix1=np.array([[1,leakage[i]],[np.conj(leakage[i]),1]],dtype=complex)
        inv_leak1=np.linalg.inv(leak_matrix1)
        for j in range(i+1,num_ant):
            leak_matrix2=np.array([[1,leakage[j]],[np.conj(leakage[j]),1]],dtype=complex)
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
    
def rotate_QU(params,Q,U,Qmodel=None,Umodel=None,return_corrected=False):
    '''
    This function rotates the Stokes vector in the UV plane by an angle theta. Theta, whe positive,
    implies rotation is in counterclockwise direction.
    :param params: This can be a float/lmfit parameter
    :param U: Observed U
    :param V: observed V
    :param I: observed I. Observed I is assumed to be the true I as well. 
                            Required if, return_corrected is False
    :param Umodel: Umodel is the predicted U due to the primary beam leakage
                   of an unpolarised source. Required if, return_corrected is False
    :param Vmodel: Vmodel is the predicted U due to the primary beam leakage
                   of an unpolarised source. Required if, return_corrected is False
    :param return_corrected: If True, the corrected U and V will be returned. By,
                                corrected, we just mean the rotated U and V. Default
                                is False
    if params is a float/list/ndarray, this function will return the sum of squared residuals
    if not, it returns array of residuals. Lmfit uses the second option, whereas scipy.minimize
    uses the first option.
    '''
    theta=params
    if not isinstance(theta,float):
        try:
            theta=params['theta'].value
        except IndexError:
            pass
    Qcor=Q*np.cos(theta)+U*np.sin(theta)
    Ucor=U*np.cos(theta)-Q*np.sin(theta)
    if return_corrected:
        return Qcor,Ucor
    
    if isinstance(params,np.ndarray) or isinstance(params,list) or isinstance(params,float):
        
        sum1=np.nansum((Qcor-Qmodel)**2+(Ucor-Umodel)**2)  ### This difference is coming from the way
                                                      ### lmfit uses the minimising function and 
   #                                                   ### how the scipy.minimize uses it.
        
    else:
        
        sum1=np.abs(Ucor-Umodel)+np.abs(Qcor-Qmodel)
        
    return (sum1)

def solve_crosshand_phase(Qmodel,Umodel,obsQ,obsU):
    res1=minimize(rotate_QU,0,args=(obsQ,obsU,\
                    Qmodel,Umodel),method='Nelder-Mead',\
                    bounds=[[-3.14159,3.14159]])
    if res1.success:
        crosshand_theta=res1.x
    else:
        crosshand_theta=np.nan
    return crosshand_theta
            
msname='polarized_source_I50_Q1_U5_V0.ms'
Umodel=5
Qmodel=1

leakage_data=np.loadtxt("quartical_ant_leaks.txt")

complex_leak=leakage_data[:,0]*np.exp(1j*leakage_data[:,1])

corrected_vis=get_leak_corrected_data(msname,complex_leak)

IQUV_corrected=compute_IQUV(corrected_vis)

crosshand_theta=solve_crosshand_phase(Qmodel,Umodel,IQUV_corrected[1],IQUV_corrected[2])

print (crosshand_theta*180/np.pi)

            
            
            
