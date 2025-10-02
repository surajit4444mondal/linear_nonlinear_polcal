import xarray,dask.array as da,dask
import numpy as np,os
from daskms.experimental.zarr import xds_from_zarr,xds_to_zarr
import matplotlib.pyplot as plt

caltable='gains_leakage.qc'
soltype='L'

gains = xds_from_zarr(caltable+"::"+soltype)
axis_names=gains[0].attrs['GAIN_AXES']
gain_data=gains[0].gains.to_numpy()
gain_flag=gains[0].gain_flags.to_numpy()
gain_flag=gain_flag.astype('bool')

for i in range(gain_data.shape[-1]):
    gain_data[...,i][gain_flag]=np.nan
    
ant_leaks_XY=gain_data[0,0,:30,0,1]
ant_leaks_YX=gain_data[0,0,:30,1,0]

amp_XY=np.abs(ant_leaks_XY)
ang_XY=np.angle(ant_leaks_XY)
amp_YX=np.abs(ant_leaks_YX)
ang_YX=np.angle(ant_leaks_YX)

#plt.plot((1j*ant_leaks).real,'o-')
#plt.show()

np.savetxt("quartical_ant_leaks.txt",np.array([amp_XY,ang_XY,amp_YX,ang_YX]).T)
