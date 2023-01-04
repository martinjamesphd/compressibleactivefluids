def rhs(omega_ft):#evaluates the right hand side in Fourier space (nonlinearity evaluated in real space)
 return -1j*(k_x_half*np.fft.rfft2(np.fft.irfft2(-1j*k_y_half*omega_ft/k_squared)*np.fft.irfft2(omega_ft))+k_y_half*np.fft.rfft2(np.fft.irfft2(1j*k_x_half*omega_ft/k_squared)*np.fft.irfft2(omega_ft)))
 
def export_and_plot(omega):
 im=ax.imshow(omega.T,interpolation='nearest',origin='lower',extent=[0,L,0,L],cmap='seismic')
 cb=fig.colorbar(im)
 plt.pause(0.0001)
 cb.remove()

import numpy as np
import matplotlib.pyplot as plt

global scale 
global L #dimesions of the system LxL=(2pi/scale)x(2pi/scale)
global N #number of grid points
global alpha #parameter
global end_time 
global dt
global k_x_half
global k_y_half
global k_squared

fig,ax=plt.subplots()
plt.ion()

scale=1
L=2*np.pi/scale
N=256
alpha=0
end_time=1000
dt=0.1
nu=0.0001

k_x_half=np.zeros((N,N//2+1))
k_y_half=np.zeros((N,N//2+1))

for i in range(N//2+1):
 for j in range(N//2+1):
  k_x_half[i,j]=i*scale
  k_y_half[i,j]=j*scale
  
for i in range(N//2+1,N):
 for j in range(N//2+1):
  k_x_half[i,j]=(i-N)*scale
  k_y_half[i,j]=j*scale
  
time=0
index=0

omega=(np.random.rand(N,N)-0.5)*3
omega_ft=np.fft.rfft2(omega)

k_matrix=-nu*(k_x_half**2+k_y_half**2)-alpha
k_matrix=np.exp(k_matrix*dt/2)
k_squared=k_x_half**2+k_y_half**2
k_squared[0,0]=1

while time<end_time:
 omega_ft_temp=k_matrix*(omega_ft+0.5*dt*rhs(omega_ft))
 omega_ft=k_matrix*(k_matrix*omega_ft+dt*rhs(omega_ft_temp))
 
 if(index%10==0):
  export_and_plot(np.fft.irfft2(omega_ft))
  print (time)

 index=index+1
 time=time+dt
