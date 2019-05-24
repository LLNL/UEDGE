def plotcontour(iso=False, title="UEDGE data"):

#import matplotlib
#import matplotlib.gridspec as gridspec
    gs=gridspec.GridSpec(2, 2)

    plt.figure(10)
    plt.subplot(gs[0,0])
    CS = plt.contour(com.rm[:,:,0], com.zm[:,:,0], bbb.te/bbb.ev)
    plt.clabel(CS, inline=1, fontsize=10)
    params = {'mathtext.default': 'regular' }          
    plt.rcParams.update(params)
    plt.title('T$\mathregular{_e}$ [ev]')
    plt.ylabel('Z [m]')
    plt.grid(True)


    plt.subplot(gs[0,1])
    CS = plt.contour(com.rm[:,:,0], com.zm[:,:,0], bbb.ti/bbb.ev)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('T$\mathregular{_i}$ [ev]')
    plt.grid(True)


    plt.subplot(gs[1,0])
    CS = plt.contour(com.rm[:,:,0], com.zm[:,:,0], bbb.ni[:,:,0]/1e20)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('N$\mathregular{_i}$/1e20 [m-3]')
    plt.xlabel('R [m]')
    plt.ylabel('Z [m]')
    plt.grid(True)


    plt.subplot(gs[1,1])
    CS = plt.contour(com.rm[:,:,0], com.zm[:,:,0], bbb.up[:,:,0]/1e3)
    plt.clabel(CS, inline=1, fontsize=10)
    plt.title('U$\mathregular{_p}$/1e3 [m/s]')
    plt.xlabel('R [m]')
    plt.grid(True)


    #####-why does not work???
    #
    #if (iso):
    #    plt.axes().set_aspect('equal', 'datalim')
    #else:
    #    plt.axes().set_aspect('auto', 'datalim')


    plt.show()

#====================================================================#
