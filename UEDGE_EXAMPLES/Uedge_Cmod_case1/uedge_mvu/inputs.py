from uedge import *

def export_inputs(fname='default_inputs.py'):
    f = open(fname, 'w')

###-Geometry parameters-###
    f.write("\n###-GEOMETRY-###\n")
    f.write("bbb.mhdgeo="+str(bbb.mhdgeo)+"\n")
    f.write("com.geometry='"+(str((com.geometry).astype('U13')[0])).strip()+"'\n")
    f.write("com.isudsym="+str(com.isudsym)+"\n")


###-Grid parameters-###
    f.write("\n###-GRID-###\n")
    f.write("com.nx="+str(com.nx)+"\n")
    f.write("com.ny="+str(com.ny)+"\n")

    f.write("flx.psi0min1="+str(flx.psi0min1)+"\n")
    f.write("flx.psi0min2="+str(flx.psi0min2)+"\n")
    f.write("flx.psi0sep="+str(flx.psi0sep)+"\n") 
    f.write("flx.psi0max="+str(flx.psi0max)+"\n") 
    f.write("flx.psi0max_inner="+str(flx.psi0max_inner)+"\n") 
    
    f.write("flx.alfcy="+str(flx.alfcy)+"\n") 
    f.write("grd.slpxt="+str(grd.slpxt)+"\n") 
    f.write("bbb.ngrid="+str(bbb.ngrid)+"\n") 
    f.write("flx.altsearch="+str(flx.altsearch)+"\n") 
    f.write("flx.istchkon="+str(flx.istchkon)+"\n") 
    f.write("flx.xcutoff1="+str(flx.xcutoff1)+"\n") 
    f.write("flx.ycutoff1="+str(flx.ycutoff1)+"\n") 

    f.write("grd.isspnew="+str(grd.isspnew)+"\n")
    f.write("grd.isztest[0]="+str(grd.isztest[0])+"\n")
    f.write("grd.isztest[1]="+str(grd.isztest[1])+"\n")
    f.write("grd.rstrike[0]="+str(grd.rstrike[0])+"\n")
    f.write("grd.rstrike[1]="+str(grd.rstrike[1])+"\n")
    f.write("grd.zstrike[0]="+str(grd.zstrike[0])+"\n")
    f.write("grd.zstrike[1]="+str(grd.zstrike[1])+"\n")


###-Solver parameters-###
    f.write("\n###-SOLVER-###\n")
    f.write("bbb.svrpkg='"+(str((bbb.svrpkg).astype('U13')[0])).strip()+"'\n")
    f.write("bbb.premeth='"+(str((bbb.premeth).astype('U13')[0])).strip()+"'\n")
    f.write("bbb.mfnksol="+str(bbb.mfnksol)+"\n")
    f.write("bbb.lfililut="+str(bbb.lfililut)+"\n")
    f.write("bbb.lenplufac="+str(bbb.lenplufac)+"\n")
    f.write("bbb.lenpfac="+str(bbb.lenpfac)+"\n")
    f.write("bbb.rlx="+str(bbb.rlx)+"\n")

    ###bbb.del = 1.e-7 # invalid syntax for python
    ###setattr(bbb, 'del', 1e-7)
    f.write("setattr(bbb,'del',"+str(getattr(bbb,'del'))+")\n")

    f.write("bbb.epscon1="+str(bbb.epscon1)+"\n")
    f.write("bbb.n0g[0]="+str(bbb.n0g[0])+"\n")

    f.write("bbb.methn="+str(bbb.methn)+"\n")
    f.write("bbb.methu="+str(bbb.methu)+"\n")
    f.write("bbb.methe="+str(bbb.methe)+"\n")
    f.write("bbb.methi="+str(bbb.methi)+"\n")
    f.write("bbb.methg="+str(bbb.methg)+"\n")




    f.close()
    print("Saved in file ", fname)


export_inputs()
