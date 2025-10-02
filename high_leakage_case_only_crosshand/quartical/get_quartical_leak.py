import xarray,dask.array as da,dask
import numpy as np,os
from daskms.experimental.zarr import xds_from_zarr,xds_to_zarr
import matplotlib.pyplot as plt
from scipy.linalg import polar

def get_polconversion(leak_matrix):
    leak_XY=np.zeros(num_ant,dtype=complex)
    leak_YX=np.zeros(num_ant,dtype=complex)
    leak_XX=np.zeros(num_ant,dtype=complex)
    leak_YY=np.zeros(num_ant,dtype=complex)

    leak_ref_ant=leak_matrix[0,0:4].reshape((2,2))

    inv_leak_ref_ant=np.linalg.inv(leak_ref_ant)

    unitary,hermitian=polar(leak_ref_ant,side='left')

    
    inv_hermitian=np.linalg.inv(hermitian)



    for i in range(num_ant):
        ant_leak=leak_matrix[i,0:4].reshape((2,2))
        ant_leak_normalized=np.matmul(ant_leak,inv_leak_ref_ant)
        hermitian_leak=np.matmul(ant_leak_normalized,hermitian)
        leak_XY[i]=hermitian_leak[0,1]
        leak_YX[i]=hermitian_leak[1,0]
        leak_XX[i]=hermitian_leak[0,0]
        leak_YY[i]=hermitian_leak[1,1]
    
    return leak_XX,leak_XY,leak_YX,leak_YY

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
    


ant_leaks_XX,ant_leaks_XY,ant_leaks_YX,ant_leaks_YY=get_polconversion(gain_data[0,0,:num_ant,0,:])



real_XY=np.real(ant_leaks_XY)
imag_XY=np.imag(ant_leaks_XY)

real_YX=np.real(ant_leaks_YX)
imag_YX=np.imag(ant_leaks_YX)

real_XX=np.real(ant_leaks_XX)
imag_XX=np.imag(ant_leaks_XX)

real_YY=np.real(ant_leaks_YY)
imag_YY=np.imag(ant_leaks_YY)

#plt.plot((1j*ant_leaks).real,'o-')
#plt.show()

np.savetxt("quartical_ant_leaks.txt",np.array([real_XX,imag_XX,\
                                                real_XY,imag_XY,\
                                                real_YX,imag_YX,\
                                                real_YY,imag_YY]).T)
