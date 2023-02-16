c!include "../mppl.h"
c-----------------------------------------------------------------------
      subroutine turb_diffus (btotyf, lte, lpi, teyf, tiyf, neyf,
     .   ted, tid, densd, mhyd, zavg, linelen,
     .   chi, lmodechin, hmodechin)
c ... Calculate anomalous radial diffusivity caused by turbulence.

      implicit none

c ... Input arguments:
      real btotyf     # total B at y-face
      real lte, lpi   # radial scale lengths of Te and ion pressure
      real teyf, tiyf, neyf   # Te, Ti, and ne at y-face
      real ted, tid, densd    # Te, Ti, and ne at divertor plate
      real mhyd       # mass of hydrogenic ions
      real zavg       # mean ion charge
      real linelen    # field-line length typical of SOL

c ... Output arguments:
      real chi                    # diffusivity (SI units)
      real lmodechin, hmodechin   # normalized diffusivities

c ... Common blocks:
      Use(Phyvar)           # ev
      Use(Turbulence)       # gammasi,kappabar,lambdap,
                            # suppress_lmode,maxmag_lmode,
                            # nky,kybeg,kyend,kya,kyb,
                            # iprint_lmode,tol_lmode,islmodebeta,
                            # gradvconst

c ... Local variables:
      real rhos
      real cubrtnu, epsilon, csed
      real betad, kt
      real gammamax, kymax

c ... Calculate quantities dependent on variables at plates.
      rhos = sqrt(mhyd * ted) / (ev * btotyf)
      cubrtnu = (2. * densd * (1. + gammasi) *
     .   sqrt(1. + (tid / ted)) * zavg * lte /
     .   (neyf * lambdap * linelen)) ** 0.333333333333333
c ... The 1.e-7 factor below accounts for use of MKS rather than cgs
      betad = 8.e-7 * Pi * densd * (ted + tid) / btotyf**2
      kt = (1. + gammasi) * zavg * 
     .                      sqrt(0.5 * betad * lambdan) / cubrtnu**2

c ... Calculate normalized L-mode diffusivity lmodechin.
      if (suppress_lmode .eq. YES) then
         lmodechin = 0.
         gammamax = 0.
         kymax = 0.
      else
         call lmode_chi_norm (kappabar, lte, rhos,
     .      cubrtnu, tiyf/ev, ted/ev, zavg, lpi, lambdap, maxmag_lmode,
     .      nky, kybeg, kyend, kya, kyb, tol_lmode, iprint_lmode,
     .      islmodebeta, kt, lmodechin, gammamax, kymax)
      endif

c ... Calculate normalized H-mode diffusivity hmodechin.
      epsilon = rhos / lte
      call hmode_chi_norm (gradvconst, cubrtnu, epsilon,
     .   lambdap, hmodechin)

c ... Convert from normalized to dimensional diffusivity chi.
      csed = sqrt(ted / mhyd)
      call turb_chi (lmodechin, hmodechin, rhos, csed, lte,
     .   lambdap, cubrtnu, chi)

      return
      end
c****** end of subroutine turb_diffus ************
c-----------------------------------------------------------------------
      subroutine lmode_chi_norm (kappabar, lte, rhos,
     .   cubrtnuarg, ti0, ted, zavg, lpi, lambdap,
     .   maxmag, nky, kybeg, kyend, kya, kyb, tol, iprint,
     .   islmodebeta, kt, lmodechin, gammamax, kymax)

c ... Calculate normalized L-mode transport coefficient chi by
c     maximizing growth rate as a function of ky.

      implicit none

c ... Input arguments:
      real kappabar     # field-line-averaged curvature [1/m]
      real lte          # L_Te = Ted / (dTed / dr0) [m]
      real rhos         # ion gyroradius at Ted [m]
      real cubrtnuarg   # cube root of collisionality nu
      real ti0          # Ti at "mid-plane" [eV]
      real ted          # Te at divertor plate [eV]
      real zavg         # average Z
      real lpi          # Pi / (dPi / dr) at "mid-plane" [m]
      real lambdap      # e (dPhi0 / dr0) / (dTed / dr0)
      real maxmag       # max magnitude of parab. step in bracketing ky
      integer nky       # number of ky's used in maximizing growth rate
      real kybeg        # lower limit of acceptable kymax [none]
      real kyend        # upper limit of acceptable kymax [none]
      real kya          # one initial point in search for kymax [none]
      real kyb          # other initial point in search for kymax [none]
      real tol          # abs & rel tolerance in search for kymax
      integer iprint
      integer islmodebeta   # =1 to turn on finite-beta correction
      real kt           # parameter in finite-beta correction [none]

c ... Output arguments:
      real lmodechin
      real gammamax     # maximum growth rate
      real kymax        # ky at maximum growth rate

c ... Common block:
      Use(Turbulence_comm)   # ssqthsqavg,kxconst,epsilon,turbdelta,
                             # cubrtnu,bcoef0,ccoef1,ccoef2,ccoef3

c ... Local variables:
      real third, fivo6
      complex ci
      real kappa, kxsq
      real a, b, c, fa, fb, fc
      integer gcona, gconb, gn

c ... External function:
      external brent, lmode_funct
      real brent, lmode_funct

c ... Set constants.
      third = 1./3.
      fivo6 = third + 0.5
      ci = cmplx(0., 1.)

c ... Copy an argument to the common block for communication to the
c     minimization routine.
      cubrtnu = cubrtnuarg

c ... Calculate variables that are independent of ky.
      kappa = 2. * (2. * ti0 / ted) * kappabar * lte
      epsilon = rhos / lte
      turbdelta = (ti0 / ted) * lte / (zavg * lpi * lambdap)
      
c ... Calculate terms and factors that enter coefficients in normalized
c     dispersion relation for L-mode turbulence.
      bcoef0 = third * ci * (lambdap * cubrtnu)**2
      ccoef1 = 2. * bcoef0 * turbdelta
      ccoef2 = fivo6 * lambdap * cubrtnu**2
      ccoef3 = kappa / (cubrtnu * lambdap)**2

c ... Find maximum growth rate by varying ky using Numerical Recipes
c     routines.
      a = kya
      b = kyb
      call mnbrak (iprint, maxmag, a, b, c, fa, fb, fc, lmode_funct)
      gammamax = -brent (iprint, nky, a, b, c, lmode_funct, tol, kymax)
      if (kymax .lt. kybeg) then
         write(*,90) '*** Max. growth rate for L-mode turbulence',
     .      ' found at ky < kybeg = ', kybeg
         call xerrab("")
      endif
      if (kymax .gt. kyend) then
         write(*,90) '*** Max. growth rate for L-mode turbulence',
     .      ' found at ky > kyend = ', kyend
         call xerrab("")
      endif
 90   format(a,a,f6.3)

c ... Compute normalized L-mode chi.
      kxsq = ssqthsqavg * kymax**2 + kxconst * (epsilon / cubrtnu)**2
      lmodechin = max(0., gammamax) / kxsq

c ... Apply finite-k_parallel (finite-beta) correction, if selected.
      if (islmodebeta .eq. 1) then
         gcona = 0.346203
         gconb = 0.0008
         if (kt .lt. 3.) then
            gn = exp(-gcona * kt**2 / (1. + gconb * kt**4))
         else
            gn = 0.
         endif
         lmodechin = lmodechin * gn
      endif

      return
      end
c-----------------------------------------------------------------------
      real function lmode_funct (ky)

c ... Return -1 times the growth rate, given ky.

      implicit none

c ... Input variable:
      real ky                # normalized poloidal wave number

c ... Common block:
      Use(Turbulence_comm)   # ssqthsqavg,kxconst,epsilon,turbdelta,
                             # cubrtnu,bcoef0,ccoef1,ccoef2,ccoef3

c ... Local variables:
      real kysq, kxsq, ksq
      complex ci
      complex bcoef, ccoef
      complex omega(2)

c ... Set constant.
      ci = cmplx(0., 1.)

c ... Calculate local variables that are dependent on ky.
      kysq = ky**2
      kxsq = ssqthsqavg * kysq + kxconst * (epsilon / cubrtnu)**2
      ksq = kxsq + kysq

c ... Calculate coefficients in normalized dispersion relation for
c     L-mode turbulence.
      bcoef = bcoef0 + 0.5 * (turbdelta * ky + ci / ksq)
      ccoef = ccoef1 * ky -
     .        (ccoef2 + ci * ky - ccoef3 * kysq) / ksq

c ... Calculate normalized roots of the L-mode dispersion relation
      call lmode_roots (bcoef, ccoef, omega)

c ... Return -1 times growth rate.
      lmode_funct = -aimag(omega(1))

      return
      end
c-----------------------------------------------------------------------
      subroutine lmode_roots (bcoef, ccoef, omega)

c ... Calculate normalized roots of the L-mode dispersion relation,
c     given the coefficients.  Root with larger growth rate is returned
c     in omega(1).

      implicit none

c ... Input arguments:
      complex bcoef, ccoef

c ... Output argument:
      complex omega(2)

c ... Local variable:
      complex discrim

c ... Use quadratic formula to evaluate roots for quadratic of the form
c        omega**2 + 2*b*omega + c = 0
      discrim = sqrt(bcoef**2 - ccoef)
      omega(1) = -bcoef + discrim
      omega(2) = -bcoef - discrim

c ... Interchange roots if necessary to get faster growing root in
c     omega(1).
      if (aimag(omega(1)) .lt. aimag(omega(2))) then
         discrim = omega(1)
         omega(1) = omega(2)
         omega(2) = discrim
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine hmode_chi_norm (gradvconst, cubrtnu, epsilon, lambdap,
     .   hmodechin)

c ... Calculate normalized H-mode transport coefficient chi.

      implicit none

c ... Input arguments:
      real gradvconst   # factor involving rad. grad. of v(parallel)
      real cubrtnu      # cube root of collisionality nu
      real epsilon      # rhos / L_Te
      real lambdap      # e (dPhi0 / dr0) / (dTed / dr0)

c ... Output argument:
      real hmodechin    # normalized chi for H-mode turbulence

      hmodechin = gradvconst * cubrtnu / (epsilon * lambdap)

      return
      end
c-----------------------------------------------------------------------
      subroutine turb_chi (lmodechin, hmodechin, rhos, csed, lte,
     .   lambdap, cubrtnu, chi)

c ... Calculate turbulent transport coefficient chi by multiplying
c     normalization constants by sum of normalized chi's from turbulence
c     presumed to dominate in L and H modes.

      implicit none

c ... Input arguments:
      real lmodechin    # normalized chi for L-mode turbulence
      real hmodechin    # normalized chi for H-mode turbulence
      real rhos         # ion gyroradius at Ted [m]
      real csed         # sound speed cs at divertor plate [m/s]
      real lte          # Ted / (dTed / dr0) [m]
      real lambdap      # e (dPhi0 / dr0) / (dTed / dr0)
      real cubrtnu      # cube root of collisionality nu

c ... Output argument:
      real chi          # turbulent transport coefficient chi (SI units)

      chi = rhos**2 * (csed / lte) * (lambdap / cubrtnu) *
     .         (lmodechin + hmodechin)

      return
      end
