import numpy as np
import matplotlib.pyplot as plt


ant_leak=np.loadtxt('../antenna_leakages.txt')
model_leak=ant_leak[:,0]*np.exp(1j*ant_leak[:,1])

model_leak_XY=np.zeros_like(model_leak)
model_leak_YX=np.zeros_like(model_leak)

model_leak_XY[:]=model_leak[:]


model_leak_XY_real=model_leak_XY.real
model_leak_XY_imag=model_leak_XY.imag

model_leak_XY_real_equal=np.linspace(np.min(model_leak_XY_real),np.max(model_leak_XY_real),5)
model_leak_XY_imag_equal=np.linspace(np.min(model_leak_XY_imag),np.max(model_leak_XY_imag),5)

model_leak_YX=np.conj(model_leak)

model_leak_YX_real=model_leak_YX.real
model_leak_YX_imag=model_leak_YX.imag

model_leak_YX_real_equal=np.linspace(np.min(model_leak_YX_real),np.max(model_leak_YX_real),5)
model_leak_YX_imag_equal=np.linspace(np.min(model_leak_YX_imag),np.max(model_leak_YX_imag),5)


quartical_leak=np.loadtxt("quartical_ant_leaks.txt")

quartical_leak_XY=quartical_leak[:,0]*np.exp(1j*quartical_leak[:,1])

quartical_leak_YX=quartical_leak[:,2]*np.exp(1j*quartical_leak[:,3])



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


#fig.savefig("quartical_solved_leakage_comparison.png")
plt.show()
