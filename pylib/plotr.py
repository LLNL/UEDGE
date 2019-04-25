# Plotting radial profiles of UEDGE data
#
# Usage example:
# plotr(bbb.te/ev, title="Te [eV]")
#
# First coding: MVU, 17-jul-17
#=========================================================#

def plotr(v, title="UEDGE data"):

    fig,ax = plt.subplots(1)

    plt.plot(com.rm[1,:,0],v[1,:])
    plt.plot(com.rm[1,:,0],v[com.nx,:])

    plt.plot(com.rm[1,:,0],v[1,:])
    plt.plot(com.rm[1,:,0],v[com.nx,:])

    plt.xlabel('R [m]')
    fig.suptitle(title)
    plt.grid(True)

    plt.show()

#=========================================================#
