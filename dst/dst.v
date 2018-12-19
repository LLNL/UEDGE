dst

#   array dimensions
{
nxm	= integer	100
nym	= integer	50
frdim	= integer	200	# number of dust radius intervals in the size distribution
fvdim	= integer	200	# number of dust speed intervals in the speed distribution
ffrdim	= integer	10000	# number of dust radius intervals in the cdf of radius
ffvdim	= integer	10000	# number of dust speed intervals in the cdf of speed
nimpm	= integer	2	# max number of impurity atom species
nimpzm	= integer	30	# max number of impurity ion species of various charge states
ncellm	= integer	400	# max number of boundary cells (2*nxm+4*nym)
naa	= integer	1	#
}


***** dustcom :

#   input geometry data
#   coordinates are in cm
#   magnetic field is in Tesla
#   Use (dustdim) first

rrm1(nxm,nym)	real	[cm]
rrm2(nxm,nym)	real	[cm]
rrm3(nxm,nym)	real	[cm]
rrm4(nxm,nym)	real	[cm]
zzm1(nxm,nym)	real	[cm]
zzm2(nxm,nym)	real	[cm]
zzm3(nxm,nym)	real	[cm]
zzm4(nxm,nym)	real	[cm]
brm1(nxm,nym)	real	[Tesla]
brm2(nxm,nym)	real	[Tesla]
brm3(nxm,nym)	real	[Tesla]
brm4(nxm,nym)	real	[Tesla]
bzm1(nxm,nym)	real	[Tesla]
bzm2(nxm,nym)	real	[Tesla]
bzm3(nxm,nym)	real	[Tesla]
bzm4(nxm,nym)	real	[Tesla]
btm1(nxm,nym)	real	[Tesla]
btm2(nxm,nym)	real	[Tesla]
btm3(nxm,nym)	real	[Tesla]
btm4(nxm,nym)	real	[Tesla]
vol(nxm,nym)	real	[cm^3]
msize(nxm,nym)	real
nx	integer
ny	integer
ixpt1(1:2)	integer
ixpt2(1:2)	integer
iysptrx(1:2)	integer
ixtop	integer
ixrb1	integer
ixlb2	integer
nw	integer
kmesh	character*5

#   arrays of cell neighbors and surfaces
ixcrs(nxm,nym,4)	integer
iycrs(nxm,nym,4)	integer
kcrs(nxm,nym,4)		integer

#   dust source mesh
nispcl(nxm,nym)		integer
nnispcl(nxm,nym)	integer
wes0(nxm,nym)	real
ndst(nxm,nym)	real
sbnd(nxm,nym)	real
fluxc(nxm,nym)	real
twall(nxm,nym)	real
score(nxm)	real
ssptr(nxm)	real

#   size and speed distributions
rdist(frdim)		real
dndist(frdim)		real
dgdist(frdim)		real
frdist(ffrdim+1)	real
fvdist(ffvdim+1)	real

#   global parameters of dust source
dstcflx		real	[atom/s]	# carbon total flux from walls as dust
dstflux		real	[part/s]	# total flux of injected dust particles
srcwi		real
srcwo		real
srcdib		real
srcdob		real
srcpib		real
srcpob		real
srcdit		real
srcdot		real
srcpit		real
srcpot		real	# partial sources

#   input plasma parameters
nem(nxm,nym)	real	[cm^-3]	# electron density
tem(nxm,nym)	real	[eV]	# electron temperature
exm(nxm,nym)	real	[V/m]
eym(nxm,nym)	real	[V/m]	# electric field components
npm(nxm,nym)	real	[cm^-3]	# proton density
tim(nxm,nym)	real	[eV]	# ion temperature
vpm(nxm,nym)	real	[cm/s]	# plasma parallel velocity
vym(nxm,nym)	real	[cm/s]	# plasma perpendicular velocity
nam(nxm,nym)	real	[cm^-3]	# atom density
tam(nxm,nym)	real	[eV]	# atom temperature
vpam(nxm,nym)	real	[cm/s]	# atom parallel velocity
vyam(nxm,nym)	real	[cm/s]	# atom perpendicular velocity
nmm(nxm,nym)	real	[cm^-3]	# molecule density
tmm(nxm,nym)	real	[eV]	# molecule temperature
vpmm(nxm,nym)	real	[cm/s]	# molecule parallel velocity
vymm(nxm,nym)	real	[cm/s]	# molecule perpendicular velocity
erm(nxm,nym)	real	[V/m]
ezm(nxm,nym)	real	[V/m]	# electric field components
pxm(nxm,nym)	real	[Pa/m]
pym(nxm,nym)	real	[Pa/m]
prm(nxm,nym)	real	[Pa/m]
pzm(nxm,nym)	real	[Pa/m]	# pressure components
vvrm(nxm,nym)	real
vvzm(nxm,nym)	real
vcvm(nxm,nym)	real

#    blob parameters
blfm(nxm,nym)	real	[1/s]	# collision frequency with blobs in the cell
tmfp	real	[s]	# characteristic time for collisions with blobs
nblb	real	[cm^-3]	# plasma density in blob
tblb	real	[erg]	# plasma temperature in blob

#   input impurity parameters:
nizm(nxm,nym,nimpzm)	real	[cm^-3]	# impurity Z ion density
tizm(nxm,nym,nimpzm)	real	[eV]	# impurity Z ion te,perature
vpizm(nxm,nym,nimpzm)	real	[cm/s]	# impurity Z ion parallel velocity
vyizm(nxm,nym,nimpzm)	real	[cm/s]	# impurity Z ion perpendicular velocity
niam(nxm,nym,nimpm)	real	[cm^-3]	# impurity atom density
tiam(nxm,nym,nimpm)	real	[eV]	# impurity atom temperature
vpiam(nxm,nym,nimpm)	real	[cm/s]	# impurity atom parallel velocity
vyiam(nxm,nym,nimpm)	real	[cm/s]	# impurity atom perpendicular velocity
mz(nimpzm)	real
mznrm(nimpzm)	real
zpot(nimpzm)	real
mia(nimpm)	real
mianrm(nimpm)	real
iz(nimpzm)	integer
zimp(nimpzm)	integer
znuc(nimpzm)	integer
ia(nimpm)	integer
zatm(nimpm)	integer
nimp	integer
nimpz	integer

#   surface fluxes parameters:
fliz(nxm,nym,nimpzm)	real	[1/s]	# impurity ion flux
fniz(nxm,nym,nimpzm)	real	[cm^-3]	# impurity ion density
flia(nxm,nym,nimpm)	real	[1/s]	# impurity atom flux
fnia(nxm,nym,nimpm)	real	[cm^-3]	# impurity atom density
flp(nxm,nym)	real	[1/s]	# plasma ion flux
fnp(nxm,nym)	real	[cm^-3]	# plasma ion density
flm(nxm,nym)	real	[1/s]	# molecule flux
fnm(nxm,nym)	real	[cm^-3]	# molecule density
fla(nxm,nym)	real	[1/s]	# atom flux
fna(nxm,nym)	real	[cm^-3]	# atom density

tp1cl(nxm,nym)	integer
tp2cl(nxm,nym)	integer
tp3cl(nxm,nym)	integer
typecl(nxm,nym)	integer
ixcl(ncellm)	integer
iycl(ncellm)	integer

#   3D vectors:
crizm(nxm,nym,nimpzm)	real	[cm/s]
czizm(nxm,nym,nimpzm)	real	[cm/s]
ctizm(nxm,nym,nimpzm)	real	[cm/s]	# impurity ion velocity components
criam(nxm,nym,nimpm)	real	[cm/s]
cziam(nxm,nym,nimpm)	real	[cm/s]
ctiam(nxm,nym,nimpm)	real	[cm/s]	# impurity atom velocity components
crm(nxm,nym)	real	[cm/s]
czm(nxm,nym)	real	[cm/s]
ctm(nxm,nym)	real	[cm/s]	# plasma velosity components
frm(nxm,nym)	real	[V/m]
fzm(nxm,nym)	real	[V/m]
ftm(nxm,nym)	real	[V/m]	# electric field components
cram(nxm,nym)	real	[cm/s]
czam(nxm,nym)	real	[cm/s]
ctam(nxm,nym)	real	[cm/s]	# atom velocity components
crmm(nxm,nym)	real	[cm/s]
czmm(nxm,nym)	real	[cm/s]
ctmm(nxm,nym)	real	[cm/s]	# molecule velocity components
pprm(nxm,nym)	real	[dyne/cm^2]
ppzm(nxm,nym)	real	[dyne/cm^2]
pptm(nxm,nym)	real	[dyne/cm^2]	# plasma pressure gradient components


***** dustcgm :
#   geometry data in a current cell (ix,iy)
#   include 'dustdim' first

icgx	integer
icgy	integer		# current cell indices
rc1	real
zc1	real
rc2	real
zc2	real
rc3	real
zc3	real
rc4	real
zc4	real
r12	real
z12	real
d12	real
r13	real
z13	real
d13	real
r34	real
z34	real
d34	real
r42	real
z42	real
d42	real
ppt	real		# distance to the boundary
rrt	real
zzt	real		# intersection point at boundary
rrtin	real
zztin	real		# the point shifted inside cell
rrtout	real
zztout	real		# the point shifted outside cell
vnr	real
vnz	real
vnt	real		# normal vector to the intersecting boundary

#   plasma data in a current cell (ix,iy)
ne	real	[cm^-3]	# electron density
te	real	[eV]	# electron temperature
np	real	[cm^-3]	# proton density
ti	real	[eV]	# proton temperature
teerg	real	[erg]
tierg	real	[erg]	# the temperatures in ergs
uzp	real
urp	real
uthp	real	[cm/s]	# proton velocity
er	real	[V/m]
ez	real	[V/m]
eth	real	[V/m]	# electric field components
na	real	[cm^-3] # atom density
ta	real	[eV]	# atom temperature
taerg	real	[erg]
uza	real	[cm/s]
ura	real	[cm/s]
utha	real	[cm/s]	# atom velocity
nm	real	[cm^-3]	# molecule density
tm	real	[eV]	# molecule temperature
uzm	real	[cm/s]
urm	real	[cm/s]
uthm	real	[cm/s]	# molecule velocity
prp	real	[dyne/cm^2]
pzp	real	[dyne/cm^2]
pthp	real	[dyne/cm^2]	# ion pressure gradient
vcv	real	[cm/s]	# convective speed, cm/s
blf	real	[1/s]	# blob collision frequency, 1/s

#   impurity data in a current cell (ix,iy)
nia(nimpm)	real	[CM^-3]	# impurity density
tia(nimpm)	real	[eV]	# impurity temperature
niz(nimpzm)	real	[cm^-3]	# impurity ion density
tiz(nimpzm)	real	[eV]	# impurity ion temperature
tiaerg(nimpm)	real	[erg]
tizerg(nimpzm)	real	[erg]
uzia(nimpm)	real	[cm/s]
uria(nimpm)	real	[cm/s]
uthia(nimpm)	real	[cm/s]	# impurity atom velocity
uziz(nimpzm)	real	[cm/s]
uriz(nimpzm)	real	[cm/s]
uthiz(nimpzm)	real	[cm/s]	# impurity ion velocity
iia(nimpm)	integer
iiz(nimpzm)	integer	# impurity atomic numbers


***** dustcntrl :

#   controll flags
jcol		integer
jstt		integer
jthe		integer
jsee		integer
jimp		integer
jsrc		integer
jdsr		integer
jrds		integer
jrdf		integer
jtheps		integer


***** dustinp :

#   input switches for dust trajectory
nisp		integer
winisp		integer
wonisp		integer
dibnisp		integer
dobnisp		integer
pibnisp		integer
pobnisp		integer
ditnisp		integer
dotnisp		integer
pitnisp		integer
potnisp		integer
jtrj		integer
jrartrj		integer
jtrjfm		integer
inipnt		integer
inivel		integer
irndsk		integer

#   dust source in a boundary cell (ixs,iys)
cells		integer
ixs		integer
iys		integer

#   file names
fnamgmt		character*256
fnamplz		character*256
fnamimp		character*256
namfo		character*256
outdir		character*256

#   input data to make dust trajectory
rd0		real
vp0		real
td0		real
xliqud0		real
radst		real
wrd		real
rpow		real
maxrd		real
minrd		real
minlogrd	real
maxlogrd	real
dlogrd		real
fdlogrd		real
dvd		real
fdvd		real
vdt		real
maxvd		real
unorm1		real
uphi1		real
dirvc		real
unorm0		real
uphi0		real
rini		real
zini		real
vrini		real
vzini		real
vtini		real
rblb		real
vblb		real
machb		real
albws		real
mlbws		real
rnvs		real
mrvs		real
tcs		real
albwl		real
mlbwl		real
rnvl		real
mrvl		real
tcl		real
lambda		real
ksip		real
ksiq		real
ksifp		real
ksifa		real
ksirn		real
ksicp		real
ksib		real
ipot		real
gsplit		real
mp		real
mpnrm		real
invmpnrm	real
tw0		real
tw04		real
tw05		real
dstsc		real
dstrc		real
flxic		real
stp0		real
ip		integer
zip		integer

#   material data
cwfunc		real
tensf		real
stens		real
md		real
Tmm0		real
tsub		real
xmeltene	real
zmat		integer
imater		integer
canmelt		integer


***** dustout :
#   output data
#   Use (dustcom) first 
#
# winj = initial statweight
# wcor = statweight delivered into the core plasma
# wspr = statweight crossing separatrix (LCFS) toward the core
# wabs = statweight absorbed by plasma in the simulated domain
# werr = statweight rejected due to internal errors
#
# 'w' - wall
# 'd' - divertor plate
# 'p' - private flux region
# 'i' - inner
# 'o' - outer
# 'b' - bottom
# 't' - top
# 'dep**' - deposited flux
# 'flx**' - incoming flux
# 'rfl**' - reflected flux
#
# xcore ,wcore   = distributions over the core interface
# xsptr ,wsptr   = distributions over the separatrix (LCFS)
# xwall ,wwall   = distributions over the wall
# xdivib,wdivib  = distributions over the inner bottom divertor plate
# xdivob,wdivob  = distributions over the outer bottom divertor plate
# xdivit,wdivit  = distributions over the inner top divertor plate
# xdivot,wdivot  = distributions over the outer top divertor plate
# xpfwl ,wpfwl   = distributions over the dome wall
#
# pasiob - flux passed through the PF cut from inner to outer bottom divertor
# pasoib - flux passed through the PF cut from outer to inner bottom divertor
# pasiot - flux passed through the PF & SOL cut from inner to outer top divertor
# pasoit - flux passed through the PF & SOL cut from outer to inner top divertor
#
# xsts,tsts,wsts - profiles of  dust parameters in cells
# nsts,rsts,msts
# tdsts,qsts,ests
# vrsts,vzsts,vtsts,vsts
# xlsts
# ltsts
#

winj		real
wcor		real
wspr		real
werr		real
wabs		real
pasiob		real
pasoib		real
pasiot		real
pasoit		real
depwi		real
depwo		real
depdib		real
depdob		real
depdit		real
depdot		real
deppib		real
deppob		real
deppit		real
deppot		real
flxwi		real
flxwo		real
flxdib		real
flxdob		real
flxdit		real
flxdot		real
flxpib		real
flxpob		real
flxpit		real
flxpot		real
rflwi		real
rflwo		real
rfldib		real
rfldob		real
rfldit		real
rfldot		real
rflpib		real
rflpob		real
rflpit		real
rflpot		real

xsts(nxm,nym)		real
tsts(nxm,nym)		real
wsts(nxm,nym)		real
nsts(nxm,nym)		real
rsts(nxm,nym)		real
msts(nxm,nym)		real
tdsts(nxm,nym)		real
qsts(nxm,nym)		real
ests(nxm,nym)		real
vrsts(nxm,nym)		real
vzsts(nxm,nym)		real
vtsts(nxm,nym)		real
vsts(nxm,nym)		real
xlsts(nxm,nym)		real
ltsts(nxm,nym)		real

xwall(nxm)		real
wwall(nxm)		real
gwall(nxm)		real
rwall(nxm)		real
mwall(nxm)		real
tdwall(nxm)		real
qwall(nxm)		real
ewall(nxm)		real
vrwall(nxm)		real
vzwall(nxm)		real
vtwall(nxm)		real
vwall(nxm)		real
xlwall(nxm)		real
ltwall(nxm)		real

xcore(nxm)		real
wcore(nxm)		real
gcore(nxm)		real
rcore(nxm)		real
mcore(nxm)		real
tdcore(nxm)		real
qcore(nxm)		real
ecore(nxm)		real
vrcore(nxm)		real
vzcore(nxm)		real
vtcore(nxm)		real
vcore(nxm)		real
xlcore(nxm)		real
ltcore(nxm)		real
xsptr(nxm)		real
wsptr(nxm)		real
gsptr(nxm)		real
rsptr(nxm)		real
msptr(nxm)		real
tdsptr(nxm)		real
qsptr(nxm)		real
esptr(nxm)		real
vrsptr(nxm)		real
vzsptr(nxm)		real
vtsptr(nxm)		real
vsptr(nxm)		real
xlsptr(nxm)		real
ltsptr(nxm)		real

xdivib(nym)		real
wdivib(nym)		real
gdivib(nym)		real
rdivib(nym)		real
mdivib(nym)		real
tddivib(nym)		real
qdivib(nym)		real
edivib(nym)		real
vrdivib(nym)		real
vzdivib(nym)		real
vtdivib(nym)		real
vdivib(nym)		real
xldivib(nym)		real
ltdivib(nym)		real

xdivit(nym)		real
wdivit(nym)		real
gdivit(nym)		real
rdivit(nym)		real
mdivit(nym)		real
tddivit(nym)		real
qdivit(nym)		real
edivit(nym)		real
vrdivit(nym)		real
vzdivit(nym)		real
vtdivit(nym)		real
vdivit(nym)		real
xldivit(nym)		real
ltdivit(nym)		real

xdivob(nym)		real
wdivob(nym)		real
gdivob(nym)		real
rdivob(nym)		real
mdivob(nym)		real
tddivob(nym)		real
qdivob(nym)		real
edivob(nym)		real
vrdivob(nym)		real
vzdivob(nym)		real
vtdivob(nym)		real
vdivob(nym)		real
xldivob(nym)		real
ltdivob(nym)		real

xdivot(nym)		real
wdivot(nym)		real
gdivot(nym)		real
rdivot(nym)		real
mdivot(nym)		real
tddivot(nym)		real
qdivot(nym)		real
edivot(nym)		real
vrdivot(nym)		real
vzdivot(nym)		real
vtdivot(nym)		real
vdivot(nym)		real
xldivot(nym)		real
ltdivot(nym)		real

xpfwl(nxm)		real
wpfwl(nxm)		real
gpfwl(nxm)		real
rpfwl(nxm)		real
mpfwl(nxm)		real
tdpfwl(nxm)		real
qpfwl(nxm)		real
epfwl(nxm)		real
vrpfwl(nxm)		real
vzpfwl(nxm)		real
vtpfwl(nxm)		real
vpfwl(nxm)		real
xlpfwl(nxm)		real
ltpfwl(nxm)		real

rdfnc(nxm,nym,0:frdim)		real
vdfnc(nxm,nym,0:fvdim)		real
vdfnctr(nxm,nym,0:fvdim)	real
vdfncpl(nxm,nym,0:fvdim)	real
rdfnct(0:frdim)			real
vdfnct(0:fvdim)			real
vdfnctrt(0:fvdim)		real
vdfncplt(0:fvdim)		real
rdnod(0:frdim)			real
vdnod(0:fvdim)			real
vtrdnod(0:fvdim)		real


***** dustcur :
#   Use (dustdim) first

ffz		real
ffr		real
fft		real
ffzz(nimpzm)	real
ffrz(nimpzm)	real
fftz(nimpzm)	real
gamz(nimpzm)	real
gamia(nimpm)	real
game		real
gamp		real
gamthe		real
gamsee		real
gama		real
aui		real
auz(nimpzm)	real


***** dusttrj :

rb		real
zb		real
thb		real
vrb		real
vzb		real
vtb		real
rdb		real
tdb		real
msdb		real
tb		real
spb		real
stb		real
fpb		real
fmb		real
qpb		real
qmb		real
fidb		real
wtb		real
r		real
rr		real
z		real
zz		real
th		real
thth		real
t		real
tt		real
sp		real
st		real
spp		real
stt		real
vr		real
vvr		real
vz		real
vvz		real
vt		real
vvt		real
rd		real
rdd		real
td		real
tdd		real
msd		real
msdd		real
fpd		real
fmd		real
qpd		real
qmd		real
fid		real
wt		real
fpdd		real
fmdd		real
qpdd		real
qmdd		real
fidd		real
wtt		real
xliqudb		real
xliqud		real
xliqudd		real
lghtb		real
lght		real
lghtd		real

kaa		integer
ra(naa)		real
za(naa)		real
tha(naa)	real
vra(naa)	real
vza(naa)	real
vta(naa)	real
rda(naa)	real
tda(naa)	real
msda(naa)	real
wta(naa)	real
tta(naa)	real
spa(naa)	real
sta(naa)	real
fpa(naa)	real
fma(naa)	real
qpa(naa)	real
qma(naa)	real
fida(naa)	real
xliquda(naa)	real
lghta(naa)	real

kjm		integer
kjms		integer


***** dst_sub :

dustt		subroutine	# main DUSTT routine

