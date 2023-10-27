import matplotlib
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from uedge import bbb, com
gs=gridspec.GridSpec(2, 2)

plt.figure(10)
plt.subplot(gs[0,0])
CS = plt.contour(com.zm[:,:,0], com.rm[:,:,0], bbb.te/bbb.ev)
plt.clabel(CS, inline=1, fontsize=10)
params = {'mathtext.default': 'regular' }          
plt.rcParams.update(params)
plt.title('T$\mathregular{_e}$ [ev]')
plt.ylabel('R [m]')
plt.grid(True)


plt.subplot(gs[0,1])
CS = plt.contour(com.zm[:,:,0], com.rm[:,:,0], bbb.ti/bbb.ev)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('T$\mathregular{_i}$ [ev]')
plt.grid(True)


plt.subplot(gs[1,0])
CS = plt.contour(com.zm[:,:,0], com.rm[:,:,0], bbb.ni[:,:,0]/1e20)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('N$\mathregular{_i}$/1e20 [m-3]')
plt.xlabel('Z [m]')
plt.ylabel('R [m]')
plt.grid(True)


plt.subplot(gs[1,1])
CS = plt.contour(com.zm[:,:,0], com.rm[:,:,0], bbb.up[:,:,0]/1e3)
plt.clabel(CS, inline=1, fontsize=10)
plt.title('U$\mathregular{_p}$/1e3 [m/s]')
plt.xlabel('Z [m]')
plt.grid(True)


plt.show()
