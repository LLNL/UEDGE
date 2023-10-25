from uedge import bbb
import numpy as np

def ModTanh(xarg, a0, a1, a2, a3, a4):
    zeta = (a2-xarg)/a3
    res = a0 + a1*((1.+a4*zeta)*np.exp(zeta) - np.exp(-zeta))/(np.exp(zeta)+np.exp(-zeta))
    return(res)


def mass_balance(ix=1, jy=1, ispec=0):

    #-flux entering from east side                                                                                    
    flux_east = -bbb.fnix[ix,jy,ispec]

    #-flux entering from west side                                                                                    
    flux_west =  bbb.fnix[ix-1,jy,ispec]

    #-flux entering from south side                                                                                   
    flux_south =  bbb.fniy[ix,jy-1,ispec]

    #-flux entering from north side                                                                                   
    flux_north = -bbb.fniy[ix,jy,ispec]


    print("Fluxes entering the cell")
    print("flux_east:", flux_east, "1/s")
    print("flux_west:", flux_west, "1/s")
    print("flux_north:", flux_north, "1/s")
    print("flux_south:", flux_south, "1/s")

    print("ion density source in the cell:", bbb.psor[ix,jy,ispec], "1/s")

    print("Sum of entering density fluxes:", flux_east+flux_west+flux_north+flux_south, "1/s")
