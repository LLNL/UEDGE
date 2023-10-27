import matplotlib.pyplot as plt
import numpy as np
from uedge import *
from uedge.uedgeplots import *



exec(open("rd_d3dHsm_in.py").read())

#-show some results
plotmesh(z_min=0., z_max=3., r_min=0.8, r_max=2.5)
plotmeshval(bbb.te/bbb.ev, z_min=0., z_max=3., r_min=0.8, r_max=2.5, title="Te [ev]")
plotmeshval(bbb.ni, z_min=0., z_max=3., r_min=0.8, r_max=2.5, title="Ni [m-3]")
plotmeshval(bbb.ng, z_min=0., z_max=3., r_min=0.8, r_max=2.5, title="Ng [m-3]")
