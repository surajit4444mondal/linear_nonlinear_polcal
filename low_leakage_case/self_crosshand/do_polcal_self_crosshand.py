import numpy as np
from casatools import table

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
    
def get_corrected_data(msname,leakage,crosshand_phase):
    
    data=get_data(msname)
    corrected_data=np.zeros_like(data)
    
    shape=data.shape
    
    num_ant=30

    k=0
    
    crosshand_matrix_inv=np.array([[np.exp(-1j*crosshand_phase),0],\
                                   [0,1]],dtype=complex)
                                   
    
    for i in range(num_ant):
        leak_matrix1=np.array([[1,leakage[i]],[np.conj(leakage[i]),1]],dtype=complex)
        inv_leak1=np.linalg.inv(leak_matrix1)
        inv_jones1=np.matmul(crosshand_matrix_inv,inv_leak1)
        for j in range(i+1,num_ant):
            leak_matrix2=np.array([[1,leakage[j]],[np.conj(leakage[j]),1]],dtype=complex)
            leak_matrix2_adj=np.conj(leak_matrix2).T
            inv_leak2=np.linalg.inv(leak_matrix2_adj)
            
            inv_jones2=np.matmul(inv_leak2,np.conj(crosshand_matrix_inv).T)
            
            vis_mat=np.array([[data[0,0,k],data[1,0,k]],[data[2,0,k],data[3,0,k]]],dtype=complex)
            
            corrected_data[:,0,k]=np.matmul(np.matmul(inv_jones1,vis_mat),inv_jones2).flatten()
            k+=1
    return corrected_data

def write_corrected_data(msname,corrected_data):
    tb=table()
    
    pos=np.where(np.isnan(corrected_data)==True)
    
    flag=np.zeros(corrected_data.shape,dtype=bool)
    
    flag[pos]=True
    
    corrected_data[pos]=0
    
    
    tb.open(msname,nomodify=False)
    try:
        tb.putcol('CORRECTED_DATA',corrected_data)
        tb.putcol('FLAG',flag)
        print ("corrected data and flags successfully updated")
    finally:
        tb.close()
        

msname='polarized_source_I50_Q1_U5_V0.ms'

leakage_data=np.loadtxt("quartical_ant_leaks.txt")
crosshand_phase=50.00131058*np.pi/180

complex_leak=leakage_data[:,0]*np.exp(1j*leakage_data[:,1])

corrected_vis=get_corrected_data(msname,complex_leak,crosshand_phase)

write_corrected_data(msname,corrected_vis)
