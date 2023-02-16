# Set params for inner (1) and outer (2) plate line-segments

from uedge import *   #needed to define variables within this file
grd.nplate1 = 5
grd.gchange("Mmod",0)

grd.rplate1=[\
   1.600E+00,  1.27300E+00,
   1.15310E+00,  1.01600E+00,  1.01600E+00]

grd.zplate1=[\
   2.34100E-01,  2.34100E-01,
   2.34100E-01,  3.71200E-01,  1.0000]

grd.nplate2 = 10
grd.gchange("Mmod,0")

grd.rplate2=[\
   2.13690E+00,  1.78570E+00,  1.76800E+00,  1.76800E+00,
   1.68100E+00,  1.67500E+00,  1.67200E+00,  1.67200E+00,  
   1.55500E+00,  1.21200E+00]
    
grd.zplate2=[\
   6.28600E-01,  4.25600E-01,  3.89300E-01,  3.46000E-01,
   3.46000E-01,  3.43000E-01,  3.37000E-01,  2.34100E-01,  
   2.34100E-01,  2.34100E-01]
