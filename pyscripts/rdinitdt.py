# Setup file to run time-dependently using dtreal
# Change dtreal for starting dt and savefname to change pfb file name
# Once variables are set, read rdrundt to execute a time-dependent run

#integer ii1,ii2,irev,ii1max,ii2max,initjac
#real8 dtreal_sav, dt_tot, mult_dt, ftol_dt, t_stop, dt_max, dtmfnk3
#real8 tstor_s, tstor_e, dt_stor, deldt_min, ftol_min, rdtphidtr
#integer n_stor, i_stor, ipt, ismfnkauto, nfe_tot
#character*5 savefname
i_stor = 0
dt_tot = 0.
nfe_tot = 0
savefname = "it333" 	# change these 5 charac. for pfb name pfdt_"savefname"
                        # also creates pftstor_... for time-dependent arrays
##dtreal = 1.e-9		# starting value of dtreal
rdtphidtr = 1e20	# ratio dtphi/dtreal
ismfnkauto = 1		# if =1, mfnksol=3 for dtreal<dtmfnk3, otherwise=-3
dtmfnk3 = 1e-4		# dtreal for mfnksol sign change if ismfnkauto=1
mult_dt = 3.4		# factor expanding dtreal after ii2max steps
ii1max = 100		# number of changes to dtreal
ii2max = 5		# number of timesteps at current dtreal
ftol_dt = 1.e-4		# fnrm tolerance for the time-dependent steps
ftol_min = 1e-9		# value of fnrm where time advance will stop
dt_tot = 0.		# total time accumulated for run (output, not input)
t_stop = 100.		# value of dt_tot (sec) where calculation will stop
dt_max = 100.		# maximum time step for dtreal
deldt_min = 0.04	# minimum relative change allowed for model_dt > 0
initjac = 0		# if=1, calc initial Jac upon reading rdcontdt
###ismmaxuc = 1		# =1 for internal calc mmaxu; =0,set mmaxu & don't chng
irev = -1		# flag to allow reduced dt advance after cutback
bbb.rlx = 0.9		# max. change in variable at each linear iteration
bbb.itermx = 10		# max. number of linear iterations allowed
tstor_s = 1e-3		# beginning time for storing solution
tstor_e = 4e-2		# ending time for storing solution
n_stor = 0		# number of linearly spaced storage points
ipt = 1			# index of one variable whose value printed at each out
			# if ipt not reset from unity, ipt=idxte(nx,iysptrx+1)

