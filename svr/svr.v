svr
 
***** UOA_routines:

xnewuoa(n:integer, npt:integer, x:double, rhobeg:double,
       rhoend:double, maxfun:integer, iprint:integer)      subroutine
######################################################################
#
# NOTE: THIS SUBROUTINE IS NOT TO BE CALLED DIRECTLY, INVOKE IT WITH
#       SCRIPT FUNCTION: "newuoa".
#
######################################################################
#
# Powell's (new) Unconstrained Optimization Algorithum, NEWUOA, seeks
# the least value of a function of many variables, by a trust region
# method that forms quadratic models by interpolation. There can be
# some freedom in the interpolation conditions, which is taken up by
# minimizing the Frobenius norm of the change to the second derivative
# of the quadratic model, beginning with a zero matrix.
#
# n      : Number of variables (n >= 2).
# npt    : Number of interpolation conditions, must be in the
#            interval: [n+2,(n+1)(n+2)/2].
# x      : Initial values of the variables must be set in x(1:n); on
#            return x will contain the values that give the least
#            calculated F.
# rhobeg and rhoend: Initial and final values of a trust region radius,
#            so both must be positive with rhoend<=rhobeg. Typically
#            rhobeg should be about one tenth of the greatest expected
#            change to a variable; rhoend should indicate the accuracy
#            that is required in the final values of the variables.
# maxfun : Upper bound on the number of calls of calfun.
# iprint : 0, 1, 2 or 3; controls the amount of printing. Specifically,
#            there is no output if iprint==0 and there is output only at
#            the return if iprint==1. Otherwise, each new value of rho
#            is printed, with the best vector of variables so far and
#            the corresponding value of the objective function. Each
#            new value of F with its variables are output if iprint==3.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# NEWUOA documentation may be obtained from the Cambridge Department
# of Applied Mathematics and Theoretical Physics/Numerical Analysis
# Department in Report NA2004/08 by M. J. D. Powell, available on the
# web at www.damtp.cam.ac.uk/user/na/na.html.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

***** UOA:
n_uoa		integer /0/
	# Number of independent variables
npt_uoa		integer /0/
	# Number of interpolation points
rhobeg_uoa	real
	# Trust-region radius, RHOBEG
rhoend_uoa	real
	# Trust-region radius, RHOEND
rho_uoa		real
	# Current value of trust-region radius, RHO
m_uoa		integer /0/
	# Length of w_uoa (storage), >= (npt+11)*(npt+n)+n*(3*n+11)/2)
f_uoa		real
	# Value of object function
x_uoa(n_uoa)	_real
	# Values of independent variables
w_uoa(m_uoa)	_real
	# Total workspace required by NEWUOA

