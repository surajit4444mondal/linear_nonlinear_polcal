from casatools import table
import numpy as np
import matplotlib.pyplot as plt


def get_data(caltable):
    tb=table()
    tb.open(caltable)
    try:
        data=tb.getcol('CPARAM')
        flag=tb.getcol('FLAG')
    finally:
        tb.close()
    
    #pos=np.where(flag==True)
    #data[pos]=np.nan
    return data[0,0,:30],data[1,0,:30]


caltable='unpolarized_source_I10.leak'

ant_leak=np.loadtxt('../antenna_leakages.txt')
model_leak_XY=ant_leak[:,0]*np.exp(1j*ant_leak[:,1])

model_leak_ref=model_leak_XY[0]

model_leak_YX=ant_leak[:,2]*np.exp(1j*ant_leak[:,3])

model_leak_XY-=model_leak_ref

model_leak_real=model_leak_XY.real
model_leak_imag=model_leak_XY.imag

model_leak_real_equal=np.linspace(np.min(model_leak_real),np.max(model_leak_real),5)
model_leak_imag_equal=np.linspace(np.min(model_leak_imag),np.max(model_leak_imag),5)

model_leak_YX+=np.conj(model_leak_ref)

model_leak_adj_real=model_leak_YX.real
model_leak_adj_imag=model_leak_YX.imag



model_leak_adj_real_equal=np.linspace(np.min(model_leak_adj_real),np.max(model_leak_adj_real),5)
model_leak_adj_imag_equal=np.linspace(np.min(model_leak_adj_imag),np.max(model_leak_adj_imag),5)

solved_leak_XY,solved_leak_YX=get_data(caltable)

fig,ax=plt.subplots(nrows=2,ncols=2,figsize=[12,6])
ax=ax.flatten()

ax[0].plot(model_leak_real,solved_leak_XY.real,'ro')
ax[0].plot(model_leak_real_equal,model_leak_real_equal,'b')
ax[0].set_xlabel("Model XY leak real")
ax[0].set_ylabel("Solved XY leak real") 

ax[1].plot(model_leak_imag,(-1j*solved_leak_XY).real,'ro')
ax[1].plot(model_leak_imag_equal,model_leak_imag_equal,'b')
ax[1].set_xlabel("Model XY leak imag")
ax[1].set_ylabel("Solved XY leak imag") 

ax[2].plot(model_leak_adj_real,solved_leak_YX.real,'ro')
ax[2].plot(model_leak_adj_real_equal,model_leak_adj_real_equal,'b')
ax[2].set_xlabel("Model YX leak real")
ax[2].set_ylabel("Solved YX leak real") 

ax[3].plot(model_leak_adj_imag,(-1j*solved_leak_YX).real,'ro')
ax[3].plot(model_leak_adj_imag_equal,model_leak_adj_imag_equal,'b')
ax[3].set_xlabel("Model YX leak imag")
ax[3].set_ylabel("Solved YX leak imag") 


fig.savefig("casa_solved_leakage_comparison.png")
plt.show()
