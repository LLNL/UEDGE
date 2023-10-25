#Maxim,There is a function rsa(te, ne, tau, k) that gives the sigma-v
#[m**3/sec] for te [Joules] and ne [/m**3].  The tau is a (real) dummy
#variable and k is the initial charge state [0].  The results depend on
#which rate tables you're using.  Usually we have istabon=10 so we're
#using the tables from Daren Stotler.e.g., at the UEDGE prompt type
#rsa(10*ev, 1e20, 0,0) yields 1.07249D-14Marv


integer nn=3
integer nt=1000

##-range of density
double nmin=1e19
double nmax=1e21

##-range of temperature
double tmin=1.0
double tmax=1000.

##-arrays of density and temperature
double narr=nmin+(nmax-nmin)*(iota(nn)-1)/(nn-1)
double tarr=tmin+(tmax-tmin)*(iota(nt)-1)/(nt-1)

##-array of ionization rate
double rsarr(nn,nt)

integer in, it

do in=1,nn
	do it=1,nt

		#in
		#it
		#narr(in)
		#tarr(it)


		double tnow=ev*tarr(it)
		double nnow=narr(in)
		##rsarr(in,it)=rsa(10*ev, 1e20, 0,0)
		rsarr(in,it)=rsa(tnow, nnow, 0,0)	
	enddo
enddo

nf
plot rsarr(1,) tarr thick=3 scale=loglin
plot rsarr(2,) tarr color=red thick=2
plot rsarr(3,) tarr color=green

#========================================================================#
