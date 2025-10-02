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

def get_corrected_data(msname,leakage_file,thetas):
    
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
    
    theta=thetas[0]
    epsilon=0#thetas[1]
    phi=0#thetas[2]
    
    theta_matrix=np.array([[np.exp(1j*theta),0],\
                            [0,1]],dtype=complex)

    epsilon_matrix=np.array([[np.cos(epsilon),1j*np.sin(epsilon)],\
                            [1j*np.sin(epsilon), np.cos(epsilon)]],dtype=complex)
    

    phi_matrix=np.array([[np.cos(phi),np.sin(phi)],\
                        [-np.sin(phi),np.cos(phi)]],dtype=complex)
                        
    polrot_matrix=np.matmul(theta_matrix,\
                            np.matmul(epsilon_matrix,phi_matrix))
                            
    
    
    
    
    
    data=get_data(msname)
    corrected_data=np.zeros_like(data)
    
    shape=data.shape
    
    num_ant=30

    k=0
    for i in range(num_ant):
        
        polconv_matrix1=np.array([[leakage_XX[i],leakage_XY[i]],[leakage_YX[i],leakage_YY[i]]],dtype=complex)
        jones_matrix1=np.matmul(polconv_matrix1,polrot_matrix)
        inv_jones_matrix1=np.linalg.inv(jones_matrix1)
        for j in range(i+1,num_ant):
            polconv_matrix2=np.array([[leakage_XX[j],leakage_XY[j]],[leakage_YX[j],leakage_YY[j]]],dtype=complex)
            jones_matrix2=np.matmul(polconv_matrix2,polrot_matrix)
            jones_matrix2_adj=np.conj(jones_matrix2).T
            inv_jones_matrix2=np.linalg.inv(jones_matrix2_adj)
            
            vis_mat=np.array([[data[0,0,k],data[1,0,k]],[data[2,0,k],data[3,0,k]]],dtype=complex)
            
            corrected_data[:,0,k]=np.matmul(np.matmul(inv_jones_matrix1,vis_mat),inv_jones_matrix2).flatten()
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


msnames=['unpolarized_source_I10.ms',\
        'polarized_source_I50_Q1_U5_V0.ms','polarized_source_I50_Q5_U1_V0.ms',\
        'polarized_source_I50_Q4_U15_V0.ms','polarized_source_I50_Q15_neg_U5_V0.ms']

leakage_file='quartical_ant_leaks.txt'
thetas=np.array([50.00131058])*np.pi/180


for msname in msnames:
    corrected_vis=get_corrected_data(msname,leakage_file,thetas)
    write_corrected_data(msname,corrected_vis)
