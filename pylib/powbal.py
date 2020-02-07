def powbal(iix,iiy):
    #
    # examining power balance in a cell (iix,iiy)
    # developed from Gary Porter's script calc_netcelleflux
    # Note: make sure to do first UEDGE> engbal(pcoree+pcorei)
    #
    # First version: MVU 6-feb-2020
    #============================================================================#

    from uedge import bbb

    #-Hydrogen radiation

    #hrad=bbb.psor[iix,iiy,0]*(bbb.eeli[iix,iiy]-bbb.ebind*bbb.ev) + \
    #      bbb.erlrc[iix,iiy] 

    #-this is a shorter way to do the same thing
    hrad=bbb.erliz[iix,iiy] + bbb.erlrc[iix,iiy] 
    print("Hydrogen radiation=", hrad, " Wt")


    if (bbb.isimpon==0):
        irad=0.0
    else:
        irad=prad[iix,iiy]*vol[iix,iiy]   # impurity radiation
    print("Impurity radiation=", irad, " Wt")


    #- net poloidal
    pflux = (bbb.feex[iix,iiy]+bbb.feix[iix,iiy]-bbb.feex[iix-1,iiy]-bbb.feix[iix-1,iiy])   
    print("Divergence of poloidal thermal current=", pflux, " Wt")

    #- net radial
    rflux = (bbb.feey[iix,iiy]+bbb.feiy[iix,iiy]-bbb.feey[iix,iiy-1]-bbb.feiy[iix,iiy-1])   
    print("Divergence of radial thermal current=", rflux, " Wt")

    mismatch=abs(pflux+rflux-hrad-irad)/(abs(hrad)+abs(irad)+abs(pflux)+abs(rflux))
    print("Relative mismatch of energy fluxes in the cell=", mismatch)
