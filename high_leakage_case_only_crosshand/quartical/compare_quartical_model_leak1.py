import xarray,dask.array as da,dask
import numpy as np,os
from daskms.experimental.zarr import xds_from_zarr,xds_to_zarr
import matplotlib.pyplot as plt
from scipy.linalg import polar


def get_polconversion(leak_matrix):
    leak_XY=np.zeros(num_ant,dtype=complex)
    leak_YX=np.zeros(num_ant,dtype=complex)

    leak_ref_ant=leak_matrix[0,0:4].reshape((2,2))

    inv_leak_ref_ant=np.linalg.inv(leak_ref_ant)

    unitary,hermitian=polar(leak_ref_ant,side='left')

    
    #inv_hermitian=np.linalg.inv(hermitian)



    for i in range(num_ant):
        ant_leak=leak_matrix[i,0:4].reshape((2,2))
        ant_leak_normalized=np.matmul(ant_leak,inv_leak_ref_ant)
        hermitian_leak=np.matmul(ant_leak_normalized,hermitian)
        leak_XY[i]=hermitian_leak[0,1]
        leak_YX[i]=hermitian_leak[1,0]
    
    return leak_XY,leak_YX
        

caltable='gains_leakage.qc'
soltype='L'
num_ant=30

gains = xds_from_zarr(caltable+"::"+soltype)
axis_names=gains[0].attrs['GAIN_AXES']
gain_data=gains[0].gains.to_numpy()
gain_flag=gains[0].gain_flags.to_numpy()
gain_flag=gain_flag.astype('bool')

for i in range(gain_data.shape[-1]):
    gain_data[...,i][gain_flag]=np.nan
    
quartical_leak_XY,quartical_leak_YX=get_polconversion(gain_data[0,0,:num_ant,0,:])

ant_leak=np.loadtxt('../antenna_leakages.txt')
model_leak_matrix=np.zeros((num_ant,4),dtype=complex)
model_leak_matrix[:,1]=ant_leak[:,0]*np.exp(1j*ant_leak[:,1])
model_leak_matrix[:,2]=ant_leak[:,2]*np.exp(1j*ant_leak[:,3])

model_leak_matrix[:,0]=1
model_leak_matrix[:,3]=1

model_leak_XY,model_leak_YX=get_polconversion(model_leak_matrix)

model_leak_XY_real=model_leak_XY.real
model_leak_XY_imag=model_leak_XY.imag

model_leak_XY_real_equal=np.linspace(np.min(model_leak_XY_real),np.max(model_leak_XY_real),5)
model_leak_XY_imag_equal=np.linspace(np.min(model_leak_XY_imag),np.max(model_leak_XY_imag),5)



model_leak_YX_real=model_leak_YX.real
model_leak_YX_imag=model_leak_YX.imag

model_leak_YX_real_equal=np.linspace(np.min(model_leak_YX_real),np.max(model_leak_YX_real),5)
model_leak_YX_imag_equal=np.linspace(np.min(model_leak_YX_imag),np.max(model_leak_YX_imag),5)





fig,ax=plt.subplots(nrows=2,ncols=2,figsize=[12,6])
ax=ax.flatten()

ax[0].plot(model_leak_XY_real,quartical_leak_XY.real,'ro')
ax[0].plot(model_leak_XY_real_equal,model_leak_XY_real_equal,'b')
ax[0].set_xlabel("Model XY leak real")
ax[0].set_ylabel("Solved XY leak real") 

ax[1].plot(model_leak_XY_imag,(-1j*quartical_leak_XY).real,'ro')
ax[1].plot(model_leak_XY_imag_equal,model_leak_XY_imag_equal,'b')
ax[1].set_xlabel("Model XY leak imag")
ax[1].set_ylabel("Solved XY leak imag") 

ax[2].plot(model_leak_YX_real,quartical_leak_YX.real,'ro')
ax[2].plot(model_leak_YX_real_equal,model_leak_YX_real_equal,'b')
ax[2].set_xlabel("Model YX leak real")
ax[2].set_ylabel("Solved YX leak real") 

ax[3].plot(model_leak_YX_imag,(-1j*quartical_leak_YX).real,'ro')
ax[3].plot(model_leak_YX_imag_equal,model_leak_YX_imag_equal,'b')
ax[3].set_xlabel("Model YX leak imag")
ax[3].set_ylabel("Solved YX leak imag") 


fig.savefig("quartical_solved_leakage_comparison.png")
plt.show()
