# Setup file to run time-dependently using dtreal
# Change dtreal for starting dt and savefname to change pfb file name
# Once variables are set, read rdrundt to execute a time-dependent run
from uedge import *

i_stor = 0
nfe_tot = 0
savefn = "savedt.hdf5"  # name of hdf5 savefile written every timestep

bbb.rdtphidtr = 1e20	# ratio dtphi/dtreal
bbb.ismfnkauto = 1		# if =1, mfnksol=3 for dtreal<dtmfnk3, otherwise=-3
bbb.dtmfnk3 = 5.e-4		# dtreal for mfnksol sign change if ismfnkauto=1
bbb.mult_dt = 3.4		# factor expanding dtreal after ii2max steps
bbb.ii1max = 500		# number of changes to dtreal
bbb.ii2max = 5		    # number of timesteps at current dtreal
bbb.itermxrdc = 7	    # value of itermx used by rdcontdt
bbb.incpset = 7		    # iterations until Jacobian is recomputed
bbb.ftol_dt = 1.e-5		# fnrm tolerance for the time-dependent steps
bbb.ftol_min = 1e-9		# value of fnrm where time advance will stop
bbb.dt_tot = 0.		    # tot time accumulated for run (output, not input)
bbb.t_stop = 100.		# value of dt_tot (sec) where calculation will stop
bbb.dt_max = 100.		# maximum time step for dtreal
bbb.dt_kill = 1e-14		# min allowed time step; rdcontdt stops if reached
bbb.deldt_min = 0.04	# minimum relative change allowed for model_dt > 0
bbb.initjac = 0		    # if=1, calc initial Jac upon reading rdcontdt
bbb.numrevjmax = 2		# number of dt reductions before Jac recalculated
bbb.numfwdjmax = 1		# number of dt increases before Jac recalculated
###bbb.ismmaxuc = 1		   # =1 for intern calc mmaxu; =0,set mmaxu & dont chng
bbb.irev = -1		    # flag to allow reduced dt advance after cutback
bbb.rlx = 0.9		    # max. change in variable at each linear iteration
bbb.itermx = 7	 	    # max. number of linear iterations allowed
bbb.tstor_s = 1e-5		# beginning time for storing solution
bbb.tstor_e = 1e-3		# ending time for storing solution
bbb.n_stor = 0		    # number of linearly spaced storage points
bbb.ipt = 1		   	    # index of variable; value printed at step
			            # if ipt not reset from unity, ipt=idxte(nx,iysptrx+1)

