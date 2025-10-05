import numpy as np
import matplotlib.pyplot as plt
from casatools import table
import os

def get_data(msname):
    tb=table()
    tb.open(msname)
    try:
        data=tb.getcol('DATA')
        corrected_data=tb.getcol('CORRECTED_DATA')
        flag=tb.getcol('FLAG')
    finally:
        tb.close()
    
    pos=np.where(flag==True)
    data[pos]=np.nan
    corrected_data[pos]=np.nan
    
    return data,corrected_data
    
def compute_IQUV(vis_data):
    avg_data=np.nanmean(vis_data,axis=(1,2))
    Idata=0.5*(avg_data[0]+avg_data[3]).real
    Vdata=0.5*(avg_data[0]-avg_data[3]).real
    Qdata=0.5*(avg_data[1]+avg_data[2]).real
    Udata=0.5*(avg_data[1]-avg_data[2]).imag
    

    return [Idata,Qdata,Udata,Vdata]


def get_polfrac_polangle(IQUV_data):
    I=IQUV_data[0]
    Q=IQUV_data[1]
    U=IQUV_data[2]
    
    pol_frac=np.sqrt(Q**2+U**2)/I
    
    pol_ang=np.arctan2(U,Q)
    
    if pol_ang<0:
        pol_ang+=np.pi
    
    #if pol_ang>np.pi/2:
    #    pol_ang-=np.pi
    
    return pol_frac,pol_ang


os.chdir('self_crosshand')
    
msname2=['polarized_source_I50_Q1_U5_V0.ms',\
        'polarized_source_I50_Q5_U1_V0.ms',\
        'polarized_source_I50_Q4_U15_V0.ms',\
        'polarized_source_I50_Q15_neg_U5_V0.ms']
        
model_flux=[[50,1,5,0],\
            [50,5,1,0],\
            [50,4,15,0],\
            [50,15,-5,0]]  #IQUV


fig,ax=plt.subplots(nrows=1,ncols=4,figsize=[12,4],subplot_kw={'projection':'polar'})

for j,(msname1,model1) in enumerate(zip(msname2,model_flux)):

    data,corrected_data=get_data(msname1)

    IQUV_data=compute_IQUV(data)
    IQUV_corrected=compute_IQUV(corrected_data)
    print (IQUV_corrected)
    polfrac_data,polang_data=get_polfrac_polangle(IQUV_data)
    polfrac_corrected,polang_corrected=get_polfrac_polangle(IQUV_corrected)
    polfrac_model,polang_model=get_polfrac_polangle(model1)
    #print (polang_model*180/np.pi,polang_corrected*180/np.pi)
    ax[j].plot([polang_data,0],[polfrac_data,0],'r',label='data')
    ax[j].plot([polang_corrected,0],[polfrac_corrected,0],'b',label='corrected')
    ax[j].plot([polang_model,0],[polfrac_model,0],'k',label='model')
    
    ax[j].legend()
    ax[j].grid('on')

#fig.savefig("compare_calibrated_values.pdf")
plt.show()
    






