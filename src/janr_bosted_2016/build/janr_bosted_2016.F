      program janr_bosted_2016

c     Originally from program /home/bosted/janr_UIM/janr_fixed.f
c     based on Inna Aznauryan JANR release contained in JANE.tar
c     dated Sept 25 2014 (also in above directory).  Peter Bosted
c     removed MINUIT fitting code and condensed separate subroutine
c     files into this package.  My only change was to restore #include
c     directives to make the code more readable.  This code was used
c     to generate predictions of beam-target and target asymmetry
c     observables for CLAS eg1-dvcs, eg1b, eg1c, eg4, e16 experiments.
c     This explains the hard-coded kinematics below.  Note in this
c     particular snapshot, the code was set up to be run separately for
c     each Q2 point, which is input using the environment variable QSQ.
      
      IMPLICIT NONE
      SAVE
      
#include "inpar.inc"
#include "gpar.inc"
#include "res_func.inc"

c hel ampl.
      REAL h1(3,2),h2(3,2),h3(3,2),h4(3,2),h5(3,2),h6(3,2)
       COMMON/HELICITY/h1,h2,h3,h4,h5,h6

      INTEGER i,j,iw,iq,it,ipi
      REAL sinth,theta,im_phi_save(6,3),phi,wfact
      REAL kgamma,qpion,k22,q22,E2pion,Eeta
      REAL W, costh, phid
      REAL Q2, S
c bins for e16
      real*8 q2e16(5)/1.72, 2.05, 2.44, 2.91, 3.48/
c bins for eg1b
      real*8 q2eg4(42)/
     >    0.0038, 0.0045, 0.0054, 0.0064, 0.0077, 
     >    0.0092, 0.0110, 0.0131, 0.0156, 0.0187, 
     >    0.0223, 0.0266, 0.0317, 0.0379, 0.0452, 
     >    0.0540, 0.0645, 0.0770, 0.0919, 0.1100, 
     >    0.1310, 0.1560, 0.1870, 0.2230, 0.2660, 
     >    0.3170, 0.3790, 0.4520, 0.5400, 0.6450, 
     >    0.7700, 0.9190, 1.1000, 1.3100, 1.5600, 
     >    1.8700, 2.2300, 2.6600, 3.1700, 3.7900, 
     >    4.52,   5.00/
c combined bins for eg4
      real q2minexcl(14)/0.00, 0.00646,     
     >    .0110, .0187, .0317, .0540, .0919, 0.156, .266, .452, 
     >    0.770, 1.31, 2.23, 3.79/
       real q2maxexcl(14)/0.00646,
     >     .0110, .0187, .0317, .0540, .0919, 0.156, .266, .452, 
     >    0.770, 1.31, 2.23, 3.79, 6.45/ ! maximum is about 1.0

c bins for eg1b
      real*8 cthmin(12)/-1., -.8, -.6, -.4, -.2, 0,
     >                   .2,  .4,  .6, .7, .8, .9/
      real*8 cthmax(12)/     -.8, -.6, -.4, -.2, 0,
     >                   .2,  .4,  .6, .7, .8, .9, 1./
c bins for eg1c
      real*8 cthminc(6)/-0.20, 0.20, 0.44, 0.63, 0.78, 0.90/
      real*8 cthmaxc(6)/       0.20, 0.44, 0.63, 0.78, 0.90, 1.0/
      character*3 firun
      character*17 fname
      integer iq2

      call getenv('QSQ',firun)
      read(firun,'(i3)') iq2
c      q2=q2eg4(iq2)
      q2 = min(4.99, 0.2 * iq2)

      write(6,'(''using q2='',f7.4,1x,i2.2)') q2,iq2
      write(fname,'(''janrpi0p'',i2.2,''.dat'')') iq2
c      write(fname,'(''janrpi0peg4'',i2.2,''.dat'')') iq2
      open(unit=7,file=fname)
      write(fname,'(''janrpipn'',i2.2,''.dat'')') iq2
c      write(fname,'(''janrpipneg4'',i2.2,''.dat'')') iq2
      open(unit=8,file=fname)

      call janr_init
      call mult_ana(q2)
      call read_multables
      call high_multables

      w = 1.5
ccc      q2 = 1.0
      costh = 0.
      phid = 0.

c      do iq=1,42
c      q2 = q2eg4(iq)
c for e16
c      do iq=1,5
c      q2 = q2e16(iq)
c for eg4
c      do iq=1,11
c      do iq=11,11
c       q2 = (q2minexcl(iq) + q2maxexcl(iq))/2.
c for eg1-dvcs
cc      do iq=1,25
      do iq=1,1
cc       q2 = 0.2 * iq
       CALL JANR_ini_point(Q2)

c for eg1b, eg1c      
      do iw=1,46
      w = 1.1 + 0.02*(iw-1)
c for e16
c      do iw=1,22
c      w = 1.1 + 0.04*(iw-1)
c for eg4
c      do iw=1,31
c       w = 1.1 + 0.03*(iw-1)

c for eg1dvcs
      do it=1,6
       costh = (cthminc(it) + cthmaxc(it))/2.
c for eg1b
c      do it = 1,12
c      costh = (cthmin(it) + cthmax(it))/2.
c for e16
c      do it = 1,8
c      costh = -1.0 + 2.0/8.*(it-0.5)
c for eg4
c      do it = 1,16
c      costh = -1.0 + 2.0/16.*(it-0.5)

      phid = 90.

      theta=acos(costh)
      theta=theta*180./pi
      sinth=sqrt(1.-costh**2)
      phi=phid*pi/180.

      s=w**2
      w2=2.*w
      E=s-mn22
      Egamma=(E-Q2)/w2
      Epion=(E+mp22)/w2
      q22=Epion**2-mp22
c      write(6,'(10f7.3)') theta,sinth,phi,s,w2,e,egamma,epion,q22

      k22=(w-Egamma)**2-mn22
      qpion=sqrt(q22)
      kgamma=sqrt(k22)
      qk=qpion*kgamma*costh-Epion*Egamma
      t=mp22-Q2+2.*qk
      t_min=mp22-Q2+2.*(qpion*kgamma-Epion*Egamma)
      u=2.*mn22+mp22-Q2-s-t
      E2pion=(E+m2pion**2)/w2
      q2pion=E2pion**2-m2pion**2
      Eeta=(E+meta**2)/w2
      qeta=Eeta**2-meta**2
c      write(6,'(10f7.3)') k22,qpion,kgamma,qk,t,
c     >  u,e2pion,q2pion,eeta,qeta

      CALL PHI_AMPL (S,costh,Q2,kgamma,qpion)

      CALL hel_ampl(S,Q2,theta,kgamma)
c      write(6,'(6e12.4)') ((h1(i,j),j=1,2),i=1,3)
c      write(6,'(6e12.4)') ((h2(i,j),j=1,2),i=1,3)
c      write(6,'(6e12.4)') ((h3(i,j),j=1,2),i=1,3)
c      write(6,'(6e12.4)') ((h4(i,j),j=1,2),i=1,3)
c      write(6,'(6e12.4)') ((h5(i,j),j=1,2),i=1,3)
c      write(6,'(6e12.4)') ((h6(i,j),j=1,2),i=1,3)
c      CALL cross_sections(S,phi,epsilon,sinth,costh)

      call response_functions (s)

      wfact = 1.
      if(w.gt.1.232) wfact = (w -1.132)/100.
      do ipi=1,2
       write(6+ipi,'(f11.6,18F8.4)')
     >  sigmat(ipi) * wfact * q2,
     >  sigmaL(ipi)/sigmaT(ipi),
     >  sigmaTL(ipi)/sigmaT(ipi),
     >  sigmaTT(ipi)/sigmaT(ipi),
     >  sigmaTLZ(ipi)/sigmaT(ipi),
     >  sigmaTTZ(ipi)/sigmaT(ipi),
     >  sigmaTLPZ(ipi)/sigmaT(ipi),
     >  sigmaTTPZ(ipi)/sigmaT(ipi),
     >  sigmaTLP(ipi)/sigmaT(ipi),
     >  sigmaTLX(ipi)/sigmaT(ipi),
     >  sigmaTLY(ipi)/sigmaT(ipi),
     >  sigmaTTX(ipi)/sigmaT(ipi),
     >  sigmaTTY(ipi)/sigmaT(ipi),
     >  sigmaTLpX(ipi)/sigmaT(ipi),
     >  sigmaTLpY(ipi)/sigmaT(ipi),
     >  sigmaTTpX(ipi)/sigmaT(ipi),
     >  sigmaLY(ipi)/sigmaT(ipi),
     >  sigmaTY(ipi)/sigmaT(ipi)
      enddo

      enddo ! loop on costh
      enddo ! loop on W
      enddo ! loop on Q2
      RETURN
      END
c
c------------------------------------------------------------------------------
      SUBROUTINE JANR_INIT()
      IMPLICIT NONE
      SAVE

#include "janr_data.inc"
#include "delta_dat.inc"
#include "inpar.inc"
#include "fit_par.inc"
#include "janrtab.inc"
#include "multipoles.inc"

      REAL at2,at3,at4,at5,at8,at9     

      CHARACTER*10 dummych
      INTEGER i,LINPUT
      data LINPUT/21/
      integer ipow, itest

      CALL lang_def
      write(6,*)'Reading the input file'
      OPEN (UNIT=LINPUT,FILE='janr_input',STATUS='OLD')
      READ (LINPUT,51)dummych,fitopt 
      Write(6,*)'Action is    ',fitopt
 51   FORMAT(A10,A10)
      READ (LINPUT,*)
      READ (LINPUT,*) Pi_type2fit, obsv2fit
C     
      IF((Pi_type2fit.LE.0).OR.(Pi_type2fit.GE.4)) THEN
         WRITE(*,*)'You choice of type of Pion is not correct'
         WRITE(*,*)'Possible choices are [1/2/3] for [Pi0/Pi+/Pi-]'
         RETURN
      ENDIF

      do ipow = 1, 4 
         itest = amod( float(obsv2fit), 10.**(ipow) ) / 10**(ipow-1) 
         if ( itest .GT. 9 ) then
            write(*,*)'Your choice of observables is incorrect'
            write(*,*)'Possible choices are 1, 10, 11, 100, 101 etc '
            return
         endif
      enddo

      READ (LINPUT,*) Ebeam
      READ (LINPUT,*) WMIN,WMAX,WSTEP
      READ (LINPUT,*) Q2MIN,Q2MAX,Q2STEP
      READ (LINPUT,*) COSMIN,COSMAX,COSSTEP
      READ (LINPUT,*) PHIMIN,PHIMAX,PHISTEP
      READ (LINPUT,*) IntAccur
      READ (LINPUT,*) tabfile
      READ (LINPUT,*) inffile
 52   FORMAT(A10,A60)
      READ (LINPUT,*)
      READ (LINPUT,*) MAXEVENTS
      READ (LINPUT,*) MAXITER
      READ (LINPUT,*) Randomseed
c     
      READ (LINPUT,*)
      READ (LINPUT,54)maxcalls
      READ (LINPUT,54)tolorance 
      READ (LINPUT,54)errorup
 54   format(f8.2)
      READ (LINPUT,*) parms
      If( parms.Le.MAXpar .And. parms.GT.0 ) Then
         READ (LINPUT,*)
         Do i=1,parms
            read(LINPUT,53)plist(i),pname(i), start_value(i)
     +           ,step_size(i),par_stat(i),low_lim(i),up_lim(i)
 53         FORMAT(i2,1x,A10,2f12.2,i10,2f10.4)
         EndDo
      Else
         Write(6,*)'Wrong Number of parameters    ',parms
         Return
      EndIf
      READ (LINPUT,*)
      READ (LINPUT,52)datform,Datafile
      READ (LINPUT,*)
      READ (LINPUT,*)lmini,lmino,lminf
      READ (LINPUT,*)chimax
      READ (LINPUT,*)Nbad
      If( Nbad.gt.0 )Then
         Do i=1,Nbad
            READ (LINPUT,*)ibad(1,i),ibad(2,i)
         EndDo
      EndIf

c     Read SAID tables
      open(UNIT=67,file='multipols/mp33.dat',status='old')
      do i=1,51
      read(67,*,end=110)wsa(i),mp33r(i),mp33i(i),at2,at3,at4
	 wsa(i)=wsa(i)/1000.
	 mp33r(i)=mp33r(i)/52.437
	 mp33i(i)=mp33i(i)/52.437
      enddo
 110  continue
      close(67)
      open(UNIT=67,file='multipols/s11.dat',status='old')
      do i=1,51
      read(67,*,end=111)wsa(i),at2,at3,at4,at5,s11r(i),s11i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 111  continue
      close(67)
      open(UNIT=67,file='multipols/s31.dat',status='old')
      do i=1,51
      read(67,*,end=112)wsa(i),at2,at3,at4,at5,s31r(i),s31i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 112  continue
      close(67)
      open(UNIT=67,file='multipols/p11.dat',status='old')
      do i=1,51
      read(67,*,end=113)wsa(i),at2,at3,at4,at5,p11r(i),p11i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 113  continue
      close(67)
      open(UNIT=67,file='multipols/p13.dat',status='old')
      do i=1,51
      read(67,*,end=114)wsa(i),at2,at3,at4,at5,p13r(i),p13i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 114  continue
      close(67)
      open(UNIT=67,file='multipols/p31.dat',status='old')
      do i=1,51
      read(67,*,end=115)wsa(i),at2,at3,at4,at5,p31r(i),p31i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 115  continue
      close(67)
      open(UNIT=67,file='multipols/p33.dat',status='old')
      do i=1,51
      read(67,*,end=116)wsa(i),at2,at3,at4,at5,p33r(i),p33i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 116  continue
      close(67)
      open(UNIT=67,file='multipols/d13.dat',status='old')
      do i=1,51
      read(67,*,end=117)wsa(i),at2,at3,at4,at5,d13r(i),d13i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 117  continue
      close(67)
      open(UNIT=67,file='multipols/d15.dat',status='old')
      do i=1,51
      read(67,*,end=118)wsa(i),at2,at3,at4,at5,d15r(i),d15i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 118  continue
      close(67)
      open(UNIT=67,file='multipols/d33.dat',status='old')
      do i=1,51
      read(67,*,end=119)wsa(i),at2,at3,at4,at5,d33r(i),d33i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 119  continue
      close(67)
      open(UNIT=67,file='multipols/f15.dat',status='old')
      do i=1,51
      read(67,*,end=120)wsa(i),at2,at3,at4,at5,f15r(i),f15i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 120  continue
      close(67)
      open(UNIT=67,file='multipols/f37.dat',status='old')
      do i=1,51
      read(67,*,end=121)wsa(i),at2,at3,at4,at5,f37r(i),f37i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 121  continue
      close(67)
      open(UNIT=67,file='multipols/pim.dat',status='old')
      do i=1,13
         read(67,*,end=122)wsa(i),pimr(i),pimi(i)
         wsa(i)=wsa(i)/1000.
      enddo
 122  continue
      close(67)
      open(UNIT=67,file='multipols/sim.dat',status='old')
      do i=1,13
         read(67,*,end=123)wsa(i),simr(i),simi(i)
         wsa(i)=wsa(i)/1000.
      enddo
 123  continue
      close(67)

      write(6,'(''Read in SAID multipoles'')')

      RETURN
      END
c
c------------------------------------------------------------------------------

      SUBROUTINE LANG_DEF
      IMPLICIT NONE
      SAVE
#include "resonances.inc"
#include "fit_par.inc"

      INTEGER i,k

      Langul(1)=1
      Langul(2)=2
      Langul(3)=3
      Langul(4)=1  
      Langul(5)=0
      Langul(6)=3
      Langul(7)=2
      Langul(8)=1    
      Langul(9)=0
      Langul(10)=1
      Langul(11)=1
      Langul(12)=2
      Langul(13)=3
C 
      Langul(21)=1
      Langul(22)=0
      Langul(23)=2  
      Langul(24)=2
      Langul(25)=0
      Langul(26)=3    
      Langul(27)=1
      Langul(28)=3
      Langul(29)=2
      Langul(30)=1
      Langul(31)=1
      Langul(32)=2  
      Langul(33)=3
C	
c  Lprime(i) is the degree of photon momenta in Breit-Wigner formula 
c  for E(l+), E(l-):
c  Lprime(i)=l for E(l+), 
c  Lprime(i)=l-2 for E(l-) at l>1,
c  Lprime(i)=l for E(l-) at l=1.    	
      Lprime(1)=1    	
      Lprime(2)=0
      Lprime(3)=3
      Lprime(4)=1
      Lprime(5)=0
      Lprime(6)=1
      Lprime(7)=2
      Lprime(8)=1
      Lprime(9)=0
      Lprime(10)=1
      Lprime(11)=1
      Lprime(12)=0    	
      Lprime(13)=1
C	
      Lprime(21)=1
      Lprime(22)=0
      Lprime(23)=0
      Lprime(24)=0
      Lprime(25)=0
      Lprime(26)=1
      Lprime(27)=1
      Lprime(28)=3
      Lprime(29)=2
      Lprime(30)=1
      Lprime(31)=1
      Lprime(32)=0    	
      Lprime(33)=1
      
c     eta are branching ratios to pion-nucleon channel
      eta(1)=1.
      eta(2)=0.15
      eta(3)=0.42
      eta(4)=0.2
      eta(5)=0.25
      eta(6)=0.1
      eta(7)=0.12
      eta(8)=0.01
      eta(9)=0.13
      eta(10)=0.16
      eta(11)=0.18
      eta(12)=0.01
      eta(13)=0.01
c
      eta(21)=0.6
      eta(22)=0.485
      eta(23)=0.5
      eta(24)=0.1
      eta(25)=0.75
      eta(26)=0.65
      eta(27)=0.15
      eta(28)=0.06
      eta(29)=0.4
      eta(30)=0.15
      eta(31)=0.2
      eta(32)=0.1
      eta(33)=0.01
c
      RETURN
      END
c
c------------------------------------------------------------------------------


      SUBROUTINE JANR_INI_point(Q2)
      IMPLICIT NONE
      SAVE
      
#include "gpar.inc"
#include "resonances.inc"
#include "amplitudes.inc"
#include "fit_par.inc"

      REAL Q2
      INTEGER i,j,k1
      REAL mpp,mp0
      REAL GM,mr,Gamma,etta,c,kr,lr,jr,qr,Er,a
      REAL A12,A32,S12,RSM,REM
       REAL AA,BB,CC,XX,YY,ZZ,SINS,SIND,COSS,COSD

      pi=2.*acos(0.E0)
      a2=sqrt(2.E0)
      mn=0.93827
      mpp=0.13957
      mp0=0.13498
        mp=(mpp+mp0)/2.
      mn22=mn**2
      mp22=mp**2
      mn2=2.*mn
      meta=0.547
      m2pion=2.*mp
c
        COSS=0.866
        SINS=0.5
        COSD=0.99
        SIND=0.11
c
CCcccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C-All resonances with (****),(***),(**) and masses below 2 GeV are included:
c
C-Resonances with isospin 3/2:
c
c-(1)	P33(1232)
c-(2)	D33(1700)
c-(3)         F37(1950)
c-(4)	P33(1600)
c-(5)	S31(1620)
c-(6)         F35(1905)
c-(7)         D35(1930)
c-(8)         P31(1925) missing resonance
c-(9)         S31(1900)
c-(10)        P33(1920)+ missing resonance (1975)
c-(11)        P31(1910)
c-(12)        D33(1940)
c-(13)        F35(2000)
c
C-Resonances with isospin 1/2:
c
c-(21)	P11(1440)
c-(22) 	S11(1535)
c-(23) 	D13(1520)
c-(24) 	D13(1700)
c-(25) 	S11(1650)
c-(26) 	F15(1680)
c-(27) 	P13(1720)
c-(28)         F17(1990)
c-(29) 	D15(1675)
c-(30)	P11(1710)
c-(31)         P13(1900)
c-(32)         D13(2080)
c-(33)         F15(2000)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC	
C 	I=3/2       c=sqrt(2./3.)
C        AM3(6)=1.5*(A12-(lr+2.)*A32)/(2.*(lr+1.))
C        AE3(6)=1.5*(lr*A32+A12)/(2.*(lr+1.))
C        AS3(6)=-1.5*S12*(2.*a*c)/sqrt(2.)/(lr+1.)
C 	I=1/2       c=-sqrt(1./3.)
C        AM1(1,7)=(A12-(lr+2.)*A32)/(2.*(lr+1.))
C        AE1(1,7)=(lr*A32+A12)/(2.*(lr+1.))
C        AS1(1,7)=-S12*(2.*a*c)/sqrt(2.)/(lr+1.)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do  i=2,40
         XRES(i)=0.5
      end do
         XRES(1)=0.2
C
c	For I=3/2 AM3=M(3/2)
c	For I=1/2 AM1(1,i)=M(1/2)/3.+M(0)
c	For I=1/2 AM1(2,i)=-M(1/2)/3.+M(0)
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C	P33(1232)
C       GM corresponds to G_M(Ash) given in Fig.10, PR C80 (2009) 055203:
        GM=exp(-0.21*Q2)/(1.+0.027334*Q2-0.008607*Q2**2)
        GM=3.*GM/(1.+Q2/0.71)**2
        resm=cc1*GM*(sqrt(((2.3933+Q2)/2.46)**2-0.88)*6.638)*1.02
c      
      AM3(1)=resm
c
C       REM and RSM are parametrizations of the results from PR C80 (2009) 055203:
        REM=-0.02*cc2
      AE3(1)=REM*AM3(1)
c
        RSM=-6.066-8.5639*Q2+2.3705755*Q2**2
        RSM=cc3*(RSM+5.807*sqrt(Q2)-0.75445*Q2**2*sqrt(Q2))/100.
      AS3(1)=RSM*AM3(1)
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C       For the resonances P11(1440),S11(1535), and D13(1520),
C       we present amplitudes, which correspond to fits to the data given 
C       in Figs.13,14,15, PR C80 (2009) 055203
C       For the resonances S11(1650),S31(1620), D13(1700) and D33(1700), we give the
C	that correspond to SQTM
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C P11(1440)
	mr=1.44
        Gamma=0.35
        etta=0.6
       c=-sqrt(1./3.)
        A12=-68.7866+21.3966*Q2+79.8415*sqrt(Q2)
        A12=A12/(1.-0.7178*Q2**2+0.5663*Q2**2*sqrt(Q2))
        A12=A12+12.*(1.-2./(1.+cc4**2))/(1.+0.25*Q2)
        S12=31.19227+3.53338*Q2
        S12=S12/(1.-0.278265*Q2**2+0.3677575*Q2**2*sqrt(Q2))
        S12=S12+8.*(1.-2./(1.+cc5**2))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/(2.*pi*qr*mr*Gamma))/50.6
        AM1(1,1)=A12*(2.*a*c)/2.
        AS1(1,1)=-S12*(2.*a*c)/sqrt(2.)
C S11(1535)
        lr=0.
        jr=1./2.
       c=-sqrt(1./3.)
        mr=1.535
        Gamma=0.15
        etta=0.485
        A12=92.5029+1.45023*Q2
        A12=A12/(1.+0.1095*Q2**2-0.000322*Q2**2*sqrt(Q2))
        XX=A12*6./COSS
        A12=A12+10.*(1.-2./(1.+cc6**2))/(1.+0.25*Q2)
        S12=-9.758811-4.231412*Q2
        S12=S12/(1.-0.7341952*Q2**2+0.5087887*Q2**2*sqrt(Q2))
        S12=S12+8.*(1.-2./(1.+cc7**2))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        AE1(1,2)=-A12*(2.*a*c)/(2.*(lr+1.))
        AS1(1,2)=-S12*(2.*a*c)/sqrt(2.)/(lr+1.)
C D13(1520) 
        lr=2.
        lr=lr-1.
        jr=3./2.
       c=-sqrt(1./3.)
        mr=1.52
        Gamma=0.12
        etta=0.5
        A12= -23.357-151.199533*Q2
        A12=A12/(1.+2.01489898*Q2**2-0.2654327*Q2**2*sqrt(Q2))
        YY=A12*6.*sqrt(2.)/COSD
        A32=3.322979*Q2-2.0339966*Q2**2+1.622563*Q2**2*sqrt(Q2)
        A32=162.458285/(1.+A32)
        ZZ=A32*2.*sqrt(6.)/COSD
        A12=A12+10.*(1.-2./(1.+cc8**2))/(1.+0.25*Q2)
        A32=A32+8.*(1.-2./(1.+cc9**2))/(1.+0.25*Q2)
        S12=1.7269289*Q2-3.25975587*Q2**2+2.905189*Q2**2*sqrt(Q2)
        S12=-67.3195419/(1.+S12)
        S12=S12+8.*(1.-2./(1.+cc10**2))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=(2.*a*c)*A12
        A32=-A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AE1(1,3)=((lr+2.)*A32-A12)/(2.*(lr+1.))
        AM1(1,3)=(A12+lr*A32)/(2.*(lr+1.))
        AS1(1,3)=-S12*(2.*a*c)/sqrt(2.)/(lr+1.)
        BB=(XX-YY)/3.
        AA=BB+(ZZ+YY)/2.
        CC=-BB+(ZZ-YY)/2.
C S11(1650)
        lr=0.
        jr=1./2.
       c=-sqrt(1./3.) 
        mr=1.65
        Gamma=0.165
        etta=0.75
        A12=(AA+BB-CC)*SINS/6.
        A12=A12+10.*(1.-2./(1.+cc11**2))/(1.+0.25*Q2)
        S12=10.*(1.-2./(1.+cc12**2))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        AE1(1,5)=-A12*(2.*a*c)/(2.*(lr+1.)) 
        AS1(1,5)=-S12*(2.*a*c)/sqrt(2.)/(lr+1.)
C S31(1620)
        lr=0.
        jr=1./2.
       c=sqrt(2./3.) 
        mr=1.63  
        Gamma=0.145       
        etta=0.25
        A12=(3.*AA-BB+CC)/18.
        A12=A12+10.*(1.-2./(1.+cc13**2))/(1.+0.25*Q2)
        S12=(-30.0+10.*(1.-2./(1.+cc14**2)))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        AE3(5)=-1.5*A12*(2.*a*c)/(2.*(lr+1.)) 
        AS3(5)=-1.5*S12*(2.*a*c)/sqrt(2.)/(lr+1.)  
C F15(1680)  The parameterization of the amplitudes corresponds to MAID global fit
        lr=3.
        lr=lr-1.
        jr=5./2.
       c=-sqrt(1./3.)
        mr=1.685
        Gamma=0.13
        etta=0.65
        A12=-(15.+70.*Q2)/(1.+1.8*Q2**2)
        A12=A12+10.*(1.-2./(1.+cc15**2))/(1.+0.25*Q2)
        A32=3.322979*Q2-2.0339966*Q2**2+1.622563*Q2**2*sqrt(Q2)
	A32=0.8*162.458285/(1.+A32)
        A32=A32+10.*(1.-2./(1.+cc16**2))/(1.+0.25*Q2)
        S12=-(10.+30.*Q2)/(1.+1.8*Q2**2)
        S12=S12+5.*(1.-2./(1.+cc17**2))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=(2.*a*c)*A12
        A32=-A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AE1(1,6)=((lr+2.)*A32-A12)/(2.*(lr+1.))
        AM1(1,6)=(A12+lr*A32)/(2.*(lr+1.))
        AS1(1,6)=-S12*(2.*a*c)/sqrt(2.)/(lr+1.)
C D33(1700) 
        lr=2.
        lr=lr-1.
        jr=3./2.
       c=sqrt(2./3.)
        mr=1.7
        Gamma=0.3
        etta=0.15
        A12=(3.*AA+2.*BB+CC)/18./sqrt(2.)
        A32=(3.*AA-CC)/6./sqrt(6.)
        A12=A12+10.*(1.-2./(1.+cc18**2))/(1.+0.25*Q2)
        A32=A32+5.*(1.-2./(1.+cc19**2))/(1.+0.25*Q2)
        S12=(10.0+10.*(1.-2./(1.+cc20**2)))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=(2.*a*c)*A12
        A32=-A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AE3(2)=1.5*((lr+2.)*A32-A12)/(2.*(lr+1.))
        AM3(2)=1.5*(A12+lr*A32)/(2.*(lr+1.))
        AS3(2)=-1.5*S12*(2.*a*c)/sqrt(2.)/(lr+1.)
C D13(1700) 
        lr=2.
        lr=lr-1.
        jr=3./2.
       c=-sqrt(1./3.)
        mr=1.7
        Gamma=0.1
        etta=0.1
        A12=-(AA-2.*BB-CC)*SIND/6./sqrt(2.)
        A32=-(AA+CC)*SIND/sqrt(6.)/2.
        A12=A12+10.*(1.-2./(1.+cc21**2))/(1.+0.25*Q2)
        A32=A32+5.*(1.-2./(1.+cc22**2))/(1.+0.25*Q2)
        S12=(10.0+10.*(1.-2./(1.+cc23**2)))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=(2.*a*c)*A12
        A32=-A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AM1(1,4)=(A12+lr*A32)/(2.*(lr+1.))
        AE1(1,4)=((lr+2.)*A32-A12)/(2.*(lr+1.))
	AS1(1,4)=-S12*(2.*a*c)/sqrt(2.)/(lr+1.)
C      P11(1710)
	mr=1.71
        Gamma=0.1
        etta=0.15
       c=-sqrt(1./3.)
        A12=0.
        S12=0.
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/(2.*pi*qr*mr*Gamma))/50.6
        AM1(1,10)=A12*(2.*a*c)/2.
        AS1(1,10)=-S12*(2.*a*c)/sqrt(2.)
C D15(1675) 
        lr=2.
        jr=5./2.
       c=-sqrt(1./3.)
        mr=1.675
        Gamma=0.15
        etta=0.4
        A12=0.
        A32=0.
        S12=0.
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=-(2.*a*c)*A12
        A32=A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AM1(1,9)=(A12-(lr+2.)*A32)/(2.*(lr+1.))
        AE1(1,9)=(lr*A32+A12)/(2.*(lr+1.))
	AS1(1,9)=-S12*(2.*a*c)/sqrt(2.)/(lr+1.)
C P13(1720) MAID global fit
        lr=1.
        jr=3./2.
       c=-sqrt(1./3.)
        mr=1.72
        Gamma=0.2
        etta=0.15
        A12=(20.+130.*Q2)/(1.+2.8*Q2**2-0.4*Q2**2*sqrt(Q2))
        A32=-(20.+220.*Q2)/(1.+4.*Q2**2)
        A12=A12+10.*(1.-2./(1.+cc24**2))/(1.+0.25*Q2)
        A32=A32+10.*(1.-2./(1.+cc25**2))/(1.+0.25*Q2)
        S12=12.*(1.-2./(1.+cc26**2))/(1.+0.25*Q2)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=-(2.*a*c)*A12
        A32=A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AM1(1,7)=(A12-(lr+2.)*A32)/(2.*(lr+1.))
        AE1(1,7)=(lr*A32+A12)/(2.*(lr+1.))
	AS1(1,7)=-S12*(2.*a*c)/sqrt(2.)/(lr+1.)
C P33(1600) 
        lr=1.
        jr=3./2.
       c=sqrt(2./3.)
        mr=1.6
        Gamma=0.3
        etta=0.2
        A12=10.*(1.-2./(1.+cc27**2))
        A32=30.*(1.-2./(1.+cc28**2))
        S12=10.*(1.-2./(1.+cc29**2))
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=-(2.*a*c)*A12
        A32=A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AM3(4)=1.5*(A12-(lr+2.)*A32)/(2.*(lr+1.))
        AE3(4)=1.5*(lr*A32+A12)/(2.*(lr+1.))
	AS3(4)=-1.5*S12*(2.*a*c)/sqrt(2.)/(lr+1.)
c-(3)         F37(1950)	****
        lr=3.
        jr=7./2.
       c=sqrt(2./3.)
        mr=1.93
        Gamma=0.3
        etta=0.4
        A12=cc30
        A32=cc31
        S12=cc32
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=-(2.*a*c)*A12
        A32=A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AM3(3)=1.5*(A12-(lr+2.)*A32)/(2.*(lr+1.))
        AE3(3)=1.5*(lr*A32+A12)/(2.*(lr+1.))
	AS3(3)=-1.5*S12*(2.*a*c)/sqrt(2.)/(lr+1.)
c-(6)         F35(1905)	****
        lr=3.
        lr=lr-1.
        jr=5./2.
       c=sqrt(2./3.)
        mr=1.89
        Gamma=0.33
        etta=0.12
        A12=cc33
        A32=cc34
        S12=cc35
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=(2.*a*c)*A12
        A32=-A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AM3(6)=1.5*(A12-(lr+2.)*A32)/(2.*(lr+1.))
        AE3(6)=1.5*(lr*A32+A12)/(2.*(lr+1.))
	AS3(6)=-1.5*S12*(2.*a*c)/sqrt(2.)/(lr+1.)
c-(7)         D35(1930)	***
        lr=2.
        jr=5./2.
       c=sqrt(2./3.)
        mr=1.960
        Gamma=0.36
        etta=0.1
        A12=cc36
        A32=cc37
        S12=cc38
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=-(2.*a*c)*A12
        A32=A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AM3(7)=1.5*(A12-(lr+2.)*A32)/(2.*(lr+1.))
        AE3(7)=1.5*(lr*A32+A12)/(2.*(lr+1.))
	AS3(7)=-1.5*S12*(2.*a*c)/sqrt(2.)/(lr+1.)
c-(10)        P33(1920)	*** + missing resonance (1975)
        lr=1.
        jr=3./2.
       c=sqrt(2./3.)
        mr=1.92
        Gamma=0.2
        etta=0.12
        A12=0.
        A32=0.
        S12=0.
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/((2.*jr+1.)*pi*qr*mr*Gamma))/50.6
        A12=-(2.*a*c)*A12
        A32=A32*(2.*a*c)/sqrt(lr*(lr+2.))
        AM3(10)=1.5*(A12-(lr+2.)*A32)/(2.*(lr+1.))
        AE3(10)=1.5*(lr*A32+A12)/(2.*(lr+1.))
	AS3(10)=-1.5*S12*(2.*a*c)/sqrt(2.)/(lr+1.)
c-(11)        P31(1910)
	mr=1.91
        Gamma=0.25
        etta=0.22
       c=sqrt(2./3.)
        kr=(mr**2-mn**2)/(2.*mr)
        Er=(mr**2+mn**2-mp**2)/(2.*mr)
        qr=sqrt(Er**2-mn**2)
        a=sqrt((kr*mn*etta)/(2.*pi*qr*mr*Gamma))/50.6
        A12=0.
        A32=0.
        S12=0.
        AM3(11)=1.5*A12*(2.*a*c)/2.
        AS3(11)=-1.5*S12*(2.*a*c)/sqrt(2.)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      MRES(1)=1.2283
      MRES(2)=1.7
      MRES(3)=1.93
      MRES(4)=1.6
      MRES(5)=1.63
      MRES(6)=1.905
      MRES(7)=1.93
      MRES(8)=1.925
      MRES(9)=1.9
      MRES(10)=1.92  
      MRES(11)=1.91
      MRES(12)=1.94  
      MRES(13)=2.
C
      MRES(21)=1.44
      MRES(22)=1.535
      MRES(23)=1.52
      MRES(24)=1.7
      MRES(25)=1.655
      MRES(26)=1.685
      MRES(27)=1.72
      MRES(28)=1.99
      MRES(29)=1.675
      MRES(30)=1.71
      MRES(31)=1.9
      MRES(32)=1.97
      MRES(33)=2.
C	
      GRES(1)=0.11
      GRES(2)=0.3
      GRES(3)=0.3
      GRES(4)=0.35
      GRES(5)=0.145
      GRES(6)=0.35
      GRES(7)=0.35
      GRES(8)=0.3
      GRES(9)=0.2
      GRES(10)=0.2
      GRES(11)=0.35
      GRES(12)=0.2
      GRES(13)=0.3
C	
      GRES(21)=0.35
      GRES(22)=0.15
      GRES(23)=0.12
      GRES(24)=0.1
      GRES(25)=0.165
      GRES(26)=0.13
      GRES(27)=0.2
      GRES(28)=0.35
      GRES(29)=0.15
      GRES(30)=0.1
      GRES(31)=0.5
      GRES(32)=0.3
      GRES(33)=0.2
C
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-(8)         P31(1925) missing resonance
c-(9)         				S31(1900) **
c-(12)        				D33(1940) *
c-(13)        F35(2000) **
Cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      DO k1=8,9
         as3(k1)=0.
         ae3(k1)=0.
         am3(k1)=0.
      ENDDO
      DO k1=12,13
         as3(k1)=0.
         ae3(k1)=0.
         am3(k1)=0.
      ENDDO
Cccccccccccccccccccccccccccccccccccc
c-(28)         F17(1990) **
c-(31)         P13(1900) **
c-(32)         				D13(2080) **
c-(33)         F15(2000) **
Cccccccccccccccccccccccccccccccccccc
c
         DO j=8,8
            am1(1,j)=0.
            ae1(1,j)=0.
            as1(1,j)=0.
         ENDDO
         DO j=11,13
            am1(1,j)=0.
            ae1(1,j)=0.
            as1(1,j)=0.
         ENDDO
c      write(6,'(''am3'',3e12.4)') (am3(i),ae3(i),as3(i),i=1,10)

      RETURN
      END
c
c------------------------------------------------------------------------------


      SUBROUTINE PHI_AMPL (S,costh,Q2,kgamma,qpion)
      IMPLICIT NONE
      SAVE
#include "multipoles.inc"
#include "highmult.inc"
#include "phiampl.inc"
#include "fit_par.inc"
      
      REAL S, Q2, kgamma, qpion ,costh
      REAL E1(2,2),E3(2)
      REAL p1,p2,p3,p4,p5
      INTEGER i, j
c
      CALL mult_ampl(S,Q2,kgamma,qpion)
      CALL high_mult(S,Q2,costh)
c
c- Calculate derivatives of Legender polinomials
c- p1=P2', p2=P3', p3=P4', p4=P3", p5=P4"
c	      
      p1=3.*costh
      p2=(15.*costh**2-3.)/2.
      p3=2.5*(7.*costh**3-3.*costh)
      p4=15.*costh
      p5=7.5*(7.*costh**2-1.)	
c      write(6,'(''ps'',10f7.3)') costh,p1,p2,p3,p4,p5
c- Calculate intermediate amplitudes for transition from 
c- multipole amplitudes to helicity amplitudes and cross sections.
c- Multipole ampl. are normalized according to CGLN Phys.Rev. 106 (1957) 1345.
c- They are equal to mult.ampl. of Devenish&Lyth/8W/pi
C       ReacM1(i,j)=0+  i=1,2 for Re,Im;        j=1,2 for pi0,pi+
C       ReacM2(i,j)=1-
C       ReacM3(i,j)=1+
C       ReacM4(i,j)=2-
C       ReacM5(i,j)=2+
C       ReacM6(i,j)=3-
C       ReacM7(i,j)=3+
      do i=1,2
c-construct the real(i=1) and imaginary(i=2) parts of intermediate amplitudes 
      do j=1,2
c-construct the amplitudes for pi0 and pi+
c	 ReacE3(i,j)=0
c	 ReacM3(i,j)=0
c	 ReacS3(i,j)=0
c	 ReacE1(i,j)=0
c	 ReacS1(i,j)=0
c	 ReacM2(i,j)=0
c	 ReacS2(i,j)=0
c	 ReacM4(i,j)=0
c	 ReacE4(i,j)=0
c	 ReacS4(i,j)=0
c	 ReacM5(i,j)=0
c	 ReacE5(i,j)=0
c	 ReacS5(i,j)=0
c	 ReacM6(i,j)=0
c	 ReacE6(i,j)=0
c	 ReacS6(i,j)=0
c	 ReacM7(i,j)=0
c	 ReacE7(i,j)=0
c	 ReacS7(i,j)=0
        E1(i,j)=ReacM6(i,j)+ReacE6(i,j)
c23456789012345678901234567890123456789012345678901234567890123456789012
        phi_amp(1,i,j)=ReacE1(i,j)+p1*(ReacM3(i,j)+ReacE3(i,j))+
     +       p3*(3.*ReacM7(i,j)+ReacE7(i,j))+
     +       3.*ReacM4(i,j)+ReacE4(i,j)+p1*(4.*ReacM6(i,j)+ReacE6(i,j))+
     +       p2*(2.*ReacM5(i,j)+ReacE5(i,j))
        phi_amp(2,i,j)=2.*ReacM3(i,j)+ReacM2(i,j)+
     +       p1*(2.*ReacM4(i,j)+3.*ReacM5(i,j))+
     +       3.*p2*ReacM6(i,j)+4.*p2*ReacM7(i,j)
        phi_amp(3,i,j)=3.*(ReacE3(i,j)-ReacM3(i,j)+E1(i,j))+
     +       p5*(ReacE7(i,j)-ReacM7(i,j))+p4*(ReacE5(i,j)-ReacM5(i,j))
        phi_amp(4,i,j)=-3.*(ReacM4(i,j)+ReacE4(i,j))+
     +       3.*(ReacM5(i,j)-ReacE5(i,j))-
     -       p4*E1(i,j)+p4*(ReacM7(i,j)-ReacE7(i,j))   
        phi_amp(5,i,j)=ReacS1(i,j)+2.*p1*ReacS3(i,j)+
     +       4.*p3*ReacS7(i,j)+3.*p2*ReacS5(i,j)-
     -       2.*ReacS4(i,j)-3.*p1*ReacS6(i,j)
        phi_amp(6,i,j)=ReacS2(i,j)+2.*p1*ReacS4(i,j)+3.*p2*ReacS6(i,j)-
     -       2.*ReacS3(i,j)-4.*p2*ReacS7(i,j)-3.*p1*ReacS5(i,j)
      end do
      end do
c
c      write(6,'(4e12.4)') ((phi_amp(1,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(2,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(3,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(4,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(5,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(6,i,j),j=1,2),i=1,2)

      do i=1,6
       do j=1,2
	phi_amp(i,1,j)=phi_amp(i,1,j)+phih(i,j)
       enddo
      enddo
c      write(6,'(4e12.4)') ((phi_amp(1,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(2,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(3,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(4,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(5,i,j),j=1,2),i=1,2)
c      write(6,'(4e12.4)') ((phi_amp(6,i,j),j=1,2),i=1,2)
c
      RETURN
      END
c
c------------------------------------------------------------------------------


      SUBROUTINE MULT_AMPL (S,Q2,kgamma,qpion)
c     ---------------------
#include "inpar.inc"
#include "gpar.inc"
#include "back_ampl.inc"
#include "multipoles.inc"
#include "amplitudes.inc"
#include "fit_par.inc"
#include "highmult.inc"

      REAL S, Q2,kgamma,qpion
      INTEGER NWRITE
      DATA NWRITE/0/
      REAL Br1(2,9),Br3(2,8)
      REAL BM33,BE23,B3M,B3Mn,B6E,B6M,B6Mn
      REAL BM11,w1,XIN1,corr
      REAL BrM1(2,2,9),BrE1(2,2,9),BrS1(2,2,9)
      REAL BrM3(2,8),BrE3(2,8),BrS3(2,8)
      REAL RMP33,IMP33,coeff
      REAL YOUTmp33r,YOUTmp33i
      REAL XIN,elab,WW,Suppress
      REAL YOUTSimr,YOUTSimi,YOUTPimr,YOUTPimi
      REAL YOUTp11r,YOUTs11r,YOUTd13r,YOUTf15r,YOUTp13r,YOUTd15r
      REAL YOUTp11i,YOUTs11i,YOUTd13i,YOUTf15i,YOUTp13i,YOUTd15i
      REAL YOUTp33r,YOUTd33r,YOUTf37r,YOUTs31r,YOUTp31r
      REAL YOUTp33i,YOUTd33i,YOUTf37i,YOUTs31i,YOUTp31i
      REAL MM1(2,2,9),MM3(2,8)
      REAL ME1(2,2,9),ME3(2,8)
      REAL MS1(2,2,9),MS3(2,8)
      REAL mpion(2),qpi(2),ratio
      REAL EV1(2),ES(2),SV1(2),SS(2),MV1(2),MS(2),SMV1(2),SMS(2)
      REAL EV3(2),SV3(2),MV3(2),SMV3(2)
      REAL Em(2),Ep(2),E00(2),Sm(2),Sp(2),S00(2)
      REAL Mm(2),Mpp(2),M0(2),SMm(2),SMp(2),SM0(2)
      INTEGER k,i,j,lun

      CALL breit_ampl(S,Q2,kgamma,qpion)
      CALL background(S,Q2)
      CALL back_high(S,Q2)
      WW=sqrt(s)
         mpion(1)=0.13498
         mpion(2)=0.13957  
	do i=1,2
      qpi(i)=sqrt(((s-mn22+mpion(i)**2)/2./WW)**2-mpion(i)**2)
	enddo
	ratio=qpi(2)/qpi(1)
      XIN=WW
      elab=(s-mn**2)/2./mn
      do i=1,2
c     -construct the real(i=1) and imaginary(i=2) parts of the
c     Breit_Wigner contributions into multipoles
         do j=1,2
c     -construct these contributions for combinations pA(1/2)(j=1),
c     nA(1/2)(j=2)
            BrM1(i,j,6)=BreitM1(i,j,6)+BreitM1(i,j,13)
            BrE1(i,j,6)=BreitE1(i,j,6)+BreitE1(i,j,13)
            BrS1(i,j,6)=BreitS1(i,j,6)+BreitS1(i,j,13)
            BrM1(i,j,7)=BreitM1(i,j,7)+BreitM1(i,j,11)    
            BrE1(i,j,7)=BreitE1(i,j,7)+BreitE1(i,j,11)    
            BrS1(i,j,7)=BreitS1(i,j,7)+BreitS1(i,j,11)    
            BrE1(i,j,2)=BreitE1(i,j,2)+BreitE1(i,j,5)
            BrS1(i,j,2)=BreitS1(i,j,2)+BreitS1(i,j,5)
            BrM1(i,j,1)=BreitM1(i,j,1)+BreitM1(i,j,10)
            BrS1(i,j,1)=BreitS1(i,j,1)+BreitS1(i,j,10)
            BrM1(i,j,3)=BreitM1(i,j,3)+BreitM1(i,j,4)+BreitM1(i,j,12)
            BrE1(i,j,3)=BreitE1(i,j,3)+BreitE1(i,j,4)+BreitE1(i,j,12)
            BrS1(i,j,3)=BreitS1(i,j,3)+BreitS1(i,j,4)+BreitS1(i,j,12)
            BrM1(i,j,9)=BreitM1(i,j,9)
            BrE1(i,j,9)=BreitE1(i,j,9)
            BrS1(i,j,9)=BreitS1(i,j,9)
            BrM1(i,j,8)=BreitM1(i,j,8)
            BrE1(i,j,8)=BreitE1(i,j,8)
            BrS1(i,j,8)=BreitS1(i,j,8)
         enddo
      enddo
c     
c 1+ P33 M1+
c	XIN1=XIN-0.004
	XIN1=XIN
      CALL interp(WSa,mp33r,51,4,XIN1,YOUTmp33r)
      RMP33=YOUTmp33r
      CALL interp(WSa,mp33i,51,4,XIN1,YOUTmp33i)
      IMP33=YOUTmp33i
c 1- P11 1(1)
      CALL interp(WSa,p11r,51,3,XIN,YOUTp11r)
      Br1(1,1)=YOUTp11r
      CALL interp(WSa,p11i,51,3,XIN,YOUTp11i)
      Br1(2,1)=YOUTp11i
c 0+ S11 1(2)
      CALL interp(WSa,s11r,51,3,XIN,YOUTs11r)
      Br1(1,2)=YOUTs11r
      CALL interp(WSa,s11i,51,3,XIN,YOUTs11i)
      Br1(2,2)=YOUTs11i
c 2- D13 1(3)
      CALL interp(WSa,d13r,51,3,XIN,YOUTd13r)
      Br1(1,3)=YOUTd13r
      CALL interp(WSa,d13i,51,3,XIN,YOUTd13i)
      Br1(2,3)=YOUTd13i
      Br1(1,4)=0.
      Br1(2,4)=0.
      Br1(1,5)=0.
      Br1(2,5)=0.
      Br1(1,8)=0.
      Br1(2,8)=0.
c 3- F15 1(6)
      CALL interp(WSa,f15r,51,3,XIN,YOUTf15r)
      Br1(1,6)=YOUTf15r
      CALL interp(WSa,f15i,51,3,XIN,YOUTf15i)
      Br1(2,6)=YOUTf15i
c 1+ P13 1(7)
      CALL interp(WSa,p13r,51,3,XIN,YOUTp13r)
      Br1(1,7)=YOUTp13r
      CALL interp(WSa,p13i,51,3,XIN,YOUTp13i)
      Br1(2,7)=YOUTp13i
c 2+ D15 1(9)
      CALL interp(WSa,d15r,51,3,XIN,YOUTd15r)
      Br1(1,9)=YOUTd15r
      CALL interp(WSa,d15i,51,3,XIN,YOUTd15i)
      Br1(2,9)=YOUTd15i

c 1+ P33 3(1)
      CALL interp(WSa,p33r,51,3,XIN,YOUTp33r)
      Br3(1,1)=YOUTp33r-0.05
      CALL interp(WSa,p33i,51,3,XIN,YOUTp33i)
      Br3(2,1)=YOUTp33i*1.017
c 2- D33 3(2)
      CALL interp(WSa,d33r,51,3,XIN,YOUTd33r)
      Br3(1,2)=YOUTd33r
      CALL interp(WSa,d33i,51,3,XIN,YOUTd33i)
      Br3(2,2)=YOUTd33i
c 3+ F37 3(3)
      CALL interp(WSa,f37r,51,3,XIN,YOUTf37r)
      Br3(1,3)=YOUTf37r
      CALL interp(WSa,f37i,51,3,XIN,YOUTf37i)
      Br3(2,3)=YOUTf37i
      Br3(1,4)=0.
      Br3(2,4)=0.
c 0+ S31 3(5)
      CALL interp(WSa,s31r,51,3,XIN,YOUTs31r)
      Br3(1,5)=YOUTs31r
      CALL interp(WSa,s31i,51,3,XIN,YOUTs31i)
      Br3(2,5)=YOUTs31i
c 3- F35 3(6)
      Br3(1,6)=0.
      Br3(2,6)=0.
c 2+ D35 3(7)
      Br3(1,7)=0.
      Br3(2,7)=0.
c 1- P31 3(8)
      CALL interp(WSa,p31r,51,3,XIN,YOUTp31r)
      Br3(1,8)=YOUTp31r
      CALL interp(WSa,p31i,51,3,XIN,YOUTp31i)
      Br3(2,8)=YOUTp31i
      do j=1,9
         do i=1,2
c     - construct Real parts of multipole amplitudes pA(1/2) and nA(1/2)
            MM1(1,i,j)=BM1(i,j)*(1.-Br1(2,j))+BrM1(1,i,j)
            ME1(1,i,j)=BE1(i,j)*(1.-Br1(2,j))+BrE1(1,i,j)
            MS1(1,i,j)=BS1(i,j)*(1.-Br1(2,j))+BrS1(1,i,j)
c     - construct Imaginary parts of multipole amplitudes pA(1/2) and nA(1/2) 
c     
            MM1(2,i,j)=BM1(i,j)*Br1(1,j)+BrM1(2,i,j)
            ME1(2,i,j)=BE1(i,j)*Br1(1,j)+BrE1(2,i,j)
            MS1(2,i,j)=BS1(i,j)*Br1(1,j)+BrS1(2,i,j)
         enddo
      enddo
c
      do i=1,2
c     -construct the real and imaginary parts of the Breit_Wigner 
c     -contributions into multipoles with I=3/2
         BrM3(i,1)=BreitM3(i,1)+BreitM3(i,4)+BreitM3(i,10)
         BrE3(i,1)=BreitE3(i,1)+BreitE3(i,4)+BreitE3(i,10)
         BrS3(i,1)=BreitS3(i,1)+BreitS3(i,4)+BreitS3(i,10)
         BrM3(i,2)=BreitM3(i,2)+BreitM3(i,12)
         BrE3(i,2)=BreitE3(i,2)+BreitE3(i,12)
         BrS3(i,2)=BreitS3(i,2)+BreitS3(i,12)
         BrM3(i,3)=BreitM3(i,3)
         BrE3(i,3)=BreitE3(i,3)
         BrS3(i,3)=BreitS3(i,3)
         BrE3(i,5)=BreitE3(i,5)+BreitE3(i,9)
         BrS3(i,5)=BreitS3(i,5)+BreitS3(i,9)
         BrM3(i,6)=BreitM3(i,6)+BreitM3(i,13)      
         BrE3(i,6)=BreitE3(i,6)+BreitE3(i,13)      
         BrS3(i,6)=BreitS3(i,6)+BreitS3(i,13)      
         BrM3(i,7)=BreitM3(i,7)
         BrE3(i,7)=BreitE3(i,7)
         BrS3(i,7)=BreitS3(i,7)
         BrM3(i,8)=BreitM3(i,8)+BreitM3(i,11)     
         BrS3(i,8)=BreitS3(i,8)+BreitS3(i,11)
      enddo     
c     - construct Real parts of multipole amplitudes A(3/2) 
      do j=1,8
         MM3(1,j)=BM3(j)*(1.-Br3(2,j))+BrM3(1,j)
         ME3(1,j)=BE3(j)*(1.-Br3(2,j))+BrE3(1,j)
         MS3(1,j)=BS3(j)*(1.-Br3(2,j))+BrS3(1,j)
c     - construct Imaginary parts of other multipole amplitudes A(3/2) 
         MM3(2,j)=BM3(j)*Br3(1,j)+BrM3(2,j)
         ME3(2,j)=BE3(j)*Br3(1,j)+BrE3(2,j)
         MS3(2,j)=BS3(j)*Br3(1,j)+BrS3(2,j)
      enddo
C----------------------------------------------------------
            MM3(1,1)=resM*RMP33
            MM3(2,1)=resM*IMP33
            do i=1,2
                MM3(i,1)=MM3(i,1)+BreitM3(i,4)
            enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C	ReacM1(i,j)=0+	i=1,2 for Re,Im;	j=1,2 for pi0,pi+	 
C	ReacM2(i,j)=1-
C	ReacM3(i,j)=1+
C	ReacM4(i,j)=2-
C	ReacM5(i,j)=2+
C	ReacM6(i,j)=3-
C	ReacM7(i,j)=3+
         do i=1,2
            ReacE1(i,1)= ME1(i,1,2)+ME3(i,5)/1.5
            ReacE1(i,2)= sqrt(2.)*(ME1(i,1,2)-ME3(i,5)/3.)
            ReacS1(i,1)= MS1(i,1,2)+MS3(i,5)/1.5
            ReacS1(i,2)= sqrt(2.)*(MS1(i,1,2)-MS3(i,5)/3.)
            ReacM2(i,1)= MM1(i,1,1)+MM3(i,8)/1.5
            ReacM2(i,2)= sqrt(2.)*(MM1(i,1,1)-MM3(i,8)/3.)*ratio
            ReacS2(i,1)= MS1(i,1,1)+MS3(i,8)/1.5
            ReacS2(i,2)= sqrt(2.)*(MS1(i,1,1)-MS3(i,8)/3.)*ratio
            ReacM3(i,1)= MM1(i,1,7)+MM3(i,1)/1.5
            ReacM3(i,2)=sqrt(2.)*(MM1(i,1,7)-MM3(i,1)/3.)*ratio
            ReacE3(i,1)= ME1(i,1,7)+ME3(i,1)/1.5
            ReacE3(i,2)= sqrt(2.)*(ME1(i,1,7)-ME3(i,1)/3.)*ratio
            ReacS3(i,1)= MS1(i,1,7)+MS3(i,1)/1.5
            ReacS3(i,2)= sqrt(2.)*(MS1(i,1,7)-MS3(i,1)/3.)*ratio
            ReacM4(i,1)= MM1(i,1,3)+MM3(i,2)/1.5
            ReacM4(i,2)= sqrt(2.)*(MM1(i,1,3)-MM3(i,2)/3.)*ratio**2
            ReacE4(i,1)= ME1(i,1,3)+ME3(i,2)/1.5
            ReacE4(i,2)= sqrt(2.)*(ME1(i,1,3)-ME3(i,2)/3.)*ratio**2
            ReacS4(i,1)= MS1(i,1,3)+MS3(i,2)/1.5
            ReacS4(i,2)= sqrt(2.)*(MS1(i,1,3)-MS3(i,2)/3.)*ratio**2
            ReacM5(i,1)= MM1(i,1,9)+MM3(i,7)/1.5
            ReacM5(i,2)= sqrt(2.)*(MM1(i,1,9)-MM3(i,7)/3.)*ratio**2
            ReacE5(i,1)= ME1(i,1,9)+ME3(i,7)/1.5
            ReacE5(i,2)= sqrt(2.)*(ME1(i,1,9)-ME3(i,7)/3.)*ratio**2
            ReacS5(i,1)= MS1(i,1,9)+MS3(i,7)/1.5
            ReacS5(i,2)= sqrt(2.)*(MS1(i,1,9)-MS3(i,7)/3.)*ratio**2
            ReacM6(i,1)= MM1(i,1,6)+MM3(i,6)/1.5
            ReacM6(i,2)= sqrt(2.)*(MM1(i,1,6)-MM3(i,6)/3.)*ratio**3
            ReacE6(i,1)= ME1(i,1,6)+ME3(i,6)/1.5
            ReacE6(i,2)= sqrt(2.)*(ME1(i,1,6)-ME3(i,6)/3.)*ratio**3
            ReacS6(i,1)= MS1(i,1,6)+MS3(i,6)/1.5
            ReacS6(i,2)= sqrt(2.)*(MS1(i,1,6)-MS3(i,6)/3.)*ratio**3
            ReacM7(i,1)= MM1(i,1,8)+MM3(i,3)/1.5
            ReacM7(i,2)= sqrt(2.)*(MM1(i,1,8)-MM3(i,3)/3.)*ratio**3
            ReacE7(i,1)= ME1(i,1,8)+ME3(i,3)/1.5
            ReacE7(i,2)= sqrt(2.)*(ME1(i,1,8)-ME3(i,3)/3.)*ratio**3
            ReacS7(i,1)= MS1(i,1,8)+MS3(i,3)/1.5
            ReacS7(i,2)= sqrt(2.)*(MS1(i,1,8)-MS3(i,3)/3.)*ratio**3
	 enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc     
        open (54,FILE='mult.dat')
      if(fitopt.eq.'GP'.or.fitopt.eq.'gp')then
        w1=ww
c        write(54,333)w1,ReacE1(1,1),ReacS1(1,1)
      endif
333   format(f4.2,16f10.5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      RETURN
      END


      SUBROUTINE BREIT_AMPL (S,Q2,kgamma, qpion)
      IMPLICIT NONE
      SAVE
#include "gpar.inc"
#include "resonances.inc"
#include "fit_par.inc"
#include "amplitudes.inc"
#include "delta_dat.inc"
      
      REAL S, Q2, kgamma, qpion
      INTEGER i, ll, j,k
      REAL M22,M2,En1,En2,kr22,qr22,kr,qr, X2, q22,k22,den
      REAL E2pion,qr2pion,r2pion,qxl,qx2l,Gpion,Ginel,Ggamma
      REAL rvirtual,rpion,WW,Wthr,Eeta,qreta,reta,Breit(2)
      REAL BreitW,Gtotal,Gtotal1,REden,IMden,qxgamma 
      REAL W2pion

c     -  Calculate Real and Imaginary parts of Breit-Wigner electroproduction amplitudes
c     -  Breit(1),Breit(2) 
c     -  Real and Imaginary parts of Breit-Wigner electroproduction multipoles
c     -  BreitM3(1,j),BreitM3(2,j) etc for I=3/2 
c     -  BreitM1(1,i,j),BreitM1(2,i,j) etc for I=1/2, 
c     -  i=1 for pM(1/2), i=2 for nM(1/2) 
        WW=sqrt(s)
      q22=qpion**2
      k22=kgamma**2
	W2pion=mn+2.*mp
      if (WW.le.W2pion) go to 166
      go to 167
166   q2pion=0.
167   continue
c     -  First do the loop for the P33(1232) resonance      
	DO i=1,1
         M22=MRES(i)**2
         M2=2.*MRES(i)
         En1=(M22+mn22+Q2)/M2  
         En2=(M22+mn22-mp22)/M2
         kr22=En1**2-mn22
         qr22=En2**2-mn22
         kr=sqrt(kr22)
         qr=sqrt(qr22)
         X2=XRES(i)**2
         rvirtual=kgamma/kr
         rpion=qpion/qr
         qxl=((X2+qr22)/(q22+X2))**Langul(i)
         Gpion=GRES(i)*rpion**(2.*Langul(i)+1.)*qxl
         Gtotal=Gpion/sqrt(s/M22)
         qxgamma=((X2+kr22)/(k22+X2))**Lprime(i)
         Ggamma=GRES(i)*rvirtual**(2.*Lprime(i)+1.)*qxgamma
         den=(M22-s)**2+(MRES(i)*Gtotal)**2
         BreitW=MRES(i)*sqrt(Gpion*Ggamma/(rvirtual*rpion))/den
         Breit(1)=BreitW*(M22-s)
         Breit(2)=BreitW*MRES(i)*Gtotal
c     -  Do the loop for real (j=1) and imaginary (j=2) parts
         do j=1,2
            BreitE3(j,i)=Breit(j)*AE3(i)
            BreitM3(j,i)=Breit(j)*AM3(i)
            BreitS3(j,i)=Breit(j)*rvirtual*AS3(i)
         enddo
      ENDDO
c     -  Do the loop for other resonances with isospin 3/2      
      DO i=2,13
         M22=MRES(i)**2	
         M2=2.*MRES(i)
         En1=(M22+mn22+Q2)/M2
         En2=(M22+mn22-mp22)/M2
         kr22=En1**2-mn22
         qr22=En2**2-mn22
         kr=sqrt(kr22)
         qr=sqrt(qr22)
         X2=XRES(i)**2
         rvirtual=kgamma/kr
         rpion=qpion/qr
         E2pion=(M22-mn22+m2pion**2)/M2
         qr2pion=E2pion**2-m2pion**2
         r2pion=(q2pion/qr2pion)**(Langul(i)+2)
         qxl=((X2+qr22)/(q22+X2))**Langul(i)
         qx2l=((X2+qr2pion)/(q2pion+X2))**(Langul(i)+2)
         Gpion=eta(i)*GRES(i)*rpion**(2.*Langul(i)+1.)*qxl         
         Ginel=(1.-eta(i))*GRES(i)*r2pion*qx2l
         Gtotal=Gpion+Ginel 
         qxgamma=((X2+kr22)/(k22+X2))**Lprime(i)
         Ggamma=GRES(i)*rvirtual**(2.*Lprime(i)+1.)*qxgamma
         den=(M22-s)**2+(MRES(i)*Gtotal)**2
         BreitW=MRES(i)*sqrt(Gpion*Ggamma/(eta(i)*rvirtual*rpion))/den
         Breit(1)=BreitW*(M22-s)
         Breit(2)=BreitW*MRES(i)*Gtotal
c     -  Do the loop for real (j=1) and imaginary (j=2) parts		 
         do j=1,2
            BreitE3(j,i)=Breit(j)*AE3(i)
            BreitM3(j,i)=Breit(j)*AM3(i)
            BreitS3(j,i)=Breit(j)*rvirtual*AS3(i)
         enddo
      ENDDO
c     
c     -  Now calculate the same for resonances with isospin 1/2.	
      DO ll=1,13
         i=ll+20
         M22=MRES(i)**2	
         M2=2.*MRES(i)
         En1=(M22+mn22+Q2)/M2
         En2=(M22+mn22-mp22)/M2
         kr22=En1**2-mn22
         qr22=En2**2-mn22
         kr=sqrt(kr22)
         qr=sqrt(qr22)
         X2=XRES(i)**2
         rvirtual=kgamma/kr
         rpion=qpion/qr
         E2pion=(M22-mn22+m2pion**2)/M2
         qr2pion=E2pion**2-m2pion**2
         r2pion=(q2pion/qr2pion)**(Langul(i)+2)
         qxl=((X2+qr22)/(q22+X2))**Langul(i)
         qx2l=((X2+qr2pion)/(q2pion+X2))**(Langul(i)+2)
         Gpion=eta(i)*GRES(i)*rpion**(2.*Langul(i)+1.)*qxl         
         Ginel=(1.-eta(i))*GRES(i)*r2pion*qx2l
         Gtotal=Gpion+Ginel
         qxgamma=((X2+kr22)/(k22+X2))**Lprime(i)
         Ggamma=GRES(i)*rvirtual**(2.*Lprime(i)+1.)*qxgamma
         den=(M22-s)**2+(MRES(i)*Gtotal)**2
         BreitW=MRES(i)*sqrt(Gpion*Ggamma/(eta(i)*rvirtual*rpion))/den
	 Breit(1)=BreitW*(M22-s)
         Breit(2)=BreitW*MRES(i)*Gtotal
c     -  Do the loop for real (j=1) and imaginary (j=2) parts	
c     -  Do the loop for  pA(1/2) (k=1) and nA(1/2)(k=2) amplitudes
         do j=1,2
            do k=1,2
               BreitE1(j,k,ll)=Breit(j)*AE1(k,ll)
               BreitM1(j,k,ll)=Breit(j)*AM1(k,ll)
               BreitS1(j,k,ll)=Breit(j)*rvirtual*AS1(k,ll)
            enddo
         enddo
      ENDDO
c     -  Now calculate the same for S11(1535) resonance		
      Wthr=mn+meta
      M22=MRES(22)**2
      M2=2.*MRES(22)
      En1=(M22+mn22+Q2)/M2
      En2=(M22+mn22-mp22)/M2
      kr22=En1**2-mn22
      qr22=En2**2-mn22
      kr=sqrt(kr22)
      qr=sqrt(qr22)
      X2=XRES(22)**2
      rvirtual=kgamma/kr 
      rpion=qpion/qr
      Eeta=(M22-mn22+meta**2)/M2
      qreta=Eeta**2-meta**2
      reta=abs(qeta/qreta)
      reta=sqrt(reta)
      Gtotal1=GRES(22)*rpion
      BreitW=MRES(22)*GRES(22)
      Gpion=0.7*GRES(22)*rpion
      Gtotal=Gpion
      REden=M22-s
      IMden=MRES(22)*Gtotal
      den=REden**2+IMden**2
      Breit(1)=BreitW*REden/den
      Breit(2)=BreitW*IMden/den
      do j=1,2
         do k=1,2
            BreitE1(j,k,2)=Breit(j)*AE1(k,2)
            BreitS1(j,k,2)=Breit(j)*rvirtual*AS1(k,2)
      enddo
      enddo
      if (WW.ge.Wthr) go to 566
      go to 567
566   IMden=MRES(22)*(Gtotal+0.3*GRES(22)*reta)
      den=REden**2+IMden**2
      Breit(1)=BreitW*REden/den
      Breit(2)=BreitW*IMden/den
      do j=1,2
         do k=1,2
            BreitE1(j,k,2)=Breit(j)*AE1(k,2)
            BreitS1(j,k,2)=Breit(j)*rvirtual*AS1(k,2)
         enddo
      enddo 
567   continue
c      write(6,'(''breite3'',10f7.3)') breite3
c      write(6,'(''breite1'',10f7.3)') breite1
c      write(6,'(''breitm3'',10f7.3)') breitm3
c      write(6,'(''breits3'',10f7.3)') breits3
c      write(6,'(''breitm1'',10f7.3)') breitm1
c      write(6,'(''breits1'',10f7.3)') breits1
c      write(6,'(''resm   '',10f7.3)') Resm
c      write(6,'(''am3    '',10f7.3)') am3 
c      write(6,'(''ae3    '',10f7.3)') ae3 
c      write(6,'(''as3    '',10f7.3)') as3 
c      write(6,'(''am1    '',10f7.3)') am1 
c      write(6,'(''ae1    '',10f7.3)') ae1 
c      write(6,'(''as1    '',10f7.3)') as1 
      RETURN
      END
c
c------------------------------------------------------------------------------

      SUBROUTINE BACKGROUND (S,Q2)
c     
      IMPLICIT NONE
      SAVE
#include "gpar.inc"
#include "back_ampl.inc"
#include "fit_par.inc"
      
      REAL S, Q2
      INTEGER NWRITE
      DATA NWRITE/0/
      INTEGER i,j
      CALL back_born(S,Q2)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	do i=1,2
		do j=1,9
       BM1(i,j)=BbM1(i,j)
       BE1(i,j)=BbE1(i,j)
       BS1(i,j)=BbS1(i,j)
c       BM1(i,j)=0.
c       BE1(i,j)=0.
c       BS1(i,j)=0.
        enddo    
	        enddo    
		do j=1,8
       BM3(j)=BbM3(j)
       BE3(j)=BbE3(j)
       BS3(j)=BbS3(j)
c       BM3(j)=0.
c       BE3(j)=0.
c       BS3(j)=0.
	        enddo   

c     
c- Background multipole amplitudes(3/2) for resonances with isospin 3/2.
c- resonances (1),(4),(10)  => (1+)
c       BM3(1),BE3(1),BS3(1)
c- resonances (2),(12)      => (2-)
c       BM3(2),BE3(2),BS3(2)
c- resonances  (3)          => (3+)
c       BM3(3),BE3(3),BS3(3)  
c- resonances  (5),(9)      => (0+)
c       BE3(5),BS3(5)
c- resonances  (6),(13)     => (3-)
c       BM3(6),BE3(6),BS3(6)
c- resonances  (7)          => (2+)  
c       BM3(7),BE3(7),BS3(7)
c- resonances  (8),(11)     => (1-)  
c       BM3(8),BE3(8),BS3(8)
           
c        do i=1,2
c- Background multipole amplitudes p(1/2) and n(1/2) 
c- for resonances with isospin 1/2.
c- resonances  (1),(10)          => (1-)
c        BM1(i,1),BS1(i,1)
c- resonances  (2),(5)           => (0+)
c        BE1(i,2),BS1(i,2)
c- resonances  (3),(4),(12)      => (2-)
c       BM1(i,3),BE1(i,3),BS1(i,3) 
c- resonances  (6),(13)          => (3-)
c       BM1(i,6),BE1(i,6),BS1(i,6) 
c- resonances  (7),(11)          => (1+)
c       BM1(i,7),BE1(i,7),BS1(i,7) 
c- resonances  (8)               => (3+)
c       BM1(i,8),BE1(i,8),BS1(i,8) 
c- resonances  (9)               => (2+)
c       BM1(i,9),BE1(i,9),BS1(i,9)   
c       enddo

      RETURN
      END

      SUBROUTINE BACK_BORN(S,Q2)
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
c BOS common block  uncomment the next line for BOS include file
#include "gpar.inc"
#include "multampl.inc"
#include "fit_par.inc"
#include "back_ampl.inc"

      REAL S, Q2
      INTEGER i,j,kpower,iw,jw,minp,maxp,mwp
      data kpower/8/
      REAL wpoint,tmult(8),wx(8),a_point,corr,cor1,cor2,qpi
      REAL bem(4,3),bep(4,3),bmm(4,3),bmp(4,3),bsm(4,3),bsp(4,3)
	REAL wdelta,del1,del2,del3,del4,del5,del6,del7,del8

c	write(6,'(''backg 1'',10f7.3)') 
c     >      BM3,BE3,BS3,
c     +      BM1,BE1,BS1,
c     +      BbM3,BbE3,BbS3,
c     +      BbM1,BbE1,BbS1,
c     +      RM3,RE3,RS3,
c     +      RM1,RE1,RS1,
c     +      E0,M1,S0,S1

      wpoint=sqrt(s)
	 qpi=sqrt(((s-mn22+mp22)/2.)**2/s-mp22)
      mwp=1
      do iw=1,maxmpoints-1
         if(wpoint.lt.wtab(iw).and.mwp.eq.1)then
            mwp=iw
            go to 101
         endif
      enddo
 101  minp=mwp-3
      if(minp.le.0)then
         minp=1
         maxp=8
      else
         maxp=mwp+4
         if(maxp.gt.maxmpoints)then
            maxp=maxmpoints
            minp=maxmpoints-7
         endif
      endif
c
      do i=0,3
         do j=1,3
c E"l"m
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=aem(iw,i+1,j)
               wx(jw)=wtab(iw)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
c            write(6,'(''tmult'',20f7.3)') (tmult(iw),iw=1,jw)
c            write(6,'(''wx'',20f7.3)') (wx(iw),iw=1,jw)
c            write(6,'(''w,a'',2f7.3)') wpoint,a_point

            bem(i+1,j)=a_point
c E"l"p
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=aep(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bep(i+1,j)=a_point
c
c M"l"m
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=amm(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bmm(i+1,j)=a_point
c M"l"p
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=amp(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bmp(i+1,j)=a_point
c S"l"m
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=asm(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bsm(i+1,j)=a_point
c S"l"p
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=asp(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bsp(i+1,j)=a_point
c
         enddo
      enddo
      CALL mixing_ampl(S,Q2)
C-constract background multipole amplitudes with I=3/2
	BbM3(1)=bmp(2,3)
	BbM3(2)=bmm(3,3)
	BbM3(3)=bmp(4,3)
	BbM3(4)=0.
	BbM3(5)=0.
	BbM3(6)=bmm(4,3)
	BbM3(7)=bmp(3,3)
	BbM3(8)=bmm(2,3)+M1(1)
c
	BbE3(1)=bep(2,3)
	BbE3(2)=bem(3,3)
	BbE3(3)=bep(4,3)
	BbE3(4)=0.
	BbE3(5)=bep(1,3)+E0(1)
	BbE3(6)=bem(4,3)
	BbE3(7)=bep(3,3)
	BbE3(8)=bem(2,3)
c
	BbS3(1)=bsp(2,3)
	BbS3(2)=bsm(3,3)
	BbS3(3)=bsp(4,3)
	BbS3(4)=0.
	BbS3(5)=bsp(1,3)+S0(1)
	BbS3(6)=bsm(4,3)
	BbS3(7)=bsp(3,3)
	BbS3(8)=bsm(2,3)+S1(1)
C-constract background multipole amplitudes with I=1/2
C-for combinations (1/2) and (0)
	BbM1(1,1)=bmm(2,2)+M1(1)
	BbM1(1,2)=0.
	BbM1(1,3)=bmm(3,2)
	BbM1(1,4)=0.
	BbM1(1,5)=0.
	BbM1(1,6)=bmm(4,2)
	BbM1(1,7)=bmp(2,2)
	BbM1(1,8)=bmp(4,2)
	BbM1(1,9)=bmp(3,2)
c
	BbM1(2,1)=bmm(2,1)+M1(2)
	BbM1(2,2)=0.
	BbM1(2,3)=bmm(3,1)
	BbM1(2,4)=0.
	BbM1(2,5)=0.
	BbM1(2,6)=bmm(4,1)
	BbM1(2,7)=bmp(2,1)
	BbM1(2,8)=bmp(4,1)
	BbM1(2,9)=bmp(3,1)
c
	BbE1(1,1)=bem(2,2)
	BbE1(1,2)=bep(1,2)+E0(1)
	BbE1(1,3)=bem(3,2)
	BbE1(1,4)=0.
	BbE1(1,5)=0.
	BbE1(1,6)=bem(4,2)
	BbE1(1,7)=bep(2,2)
	BbE1(1,8)=bep(4,2)
	BbE1(1,9)=bep(3,2)
c
	BbE1(2,1)=bem(2,1)
	BbE1(2,2)=bep(1,1)+E0(2)
	BbE1(2,3)=bem(3,1)
	BbE1(2,4)=0.
	BbE1(2,5)=0.
	BbE1(2,6)=bem(4,1)
	BbE1(2,7)=bep(2,1)
	BbE1(2,8)=bep(4,1)
	BbE1(2,9)=bep(3,1)
c
	BbS1(1,1)=bsm(2,2)+S1(1)
	BbS1(1,2)=bsp(1,2)+S0(1)
	BbS1(1,3)=bsm(3,2)
	BbS1(1,4)=0.
	BbS1(1,5)=0.
	BbS1(1,6)=bsm(4,2)
	BbS1(1,7)=bsp(2,2)
	BbS1(1,8)=bsp(4,2)
	BbS1(1,9)=bsp(3,2)
c
	BbS1(2,1)=bsm(2,1)+S1(2)
	BbS1(2,2)=bsp(1,1)+S0(2)
	BbS1(2,3)=bsm(3,1)
	BbS1(2,4)=0.
	BbS1(2,5)=0.
	BbS1(2,6)=bsm(4,1)
	BbS1(2,7)=bsp(2,1)
	BbS1(2,8)=bsp(4,1)
	BbS1(2,9)=bsp(3,1)
C-construct background multipole amplitudes with I=1/2
C-for combinations p(1/2) and n(1/2)
	do j=1,9
	BbM1(1,j)=BbM1(1,j)/3.+BbM1(2,j)
	BbM1(2,j)=2.*BbM1(2,j)-BbM1(1,j)
	BbE1(1,j)=BbE1(1,j)/3.+BbE1(2,j)
	BbE1(2,j)=2.*BbE1(2,j)-BbE1(1,j)
	BbS1(1,j)=BbS1(1,j)/3.+BbS1(2,j)
	BbS1(2,j)=2.*BbS1(2,j)-BbS1(1,j)
	enddo
       cor2=(1.1/Wpoint)*3./Q2
       BbE3(5)=BbE3(5)+0.05*cor2*(1.-2./(1.+cc41**2))
       BbS3(5)=BbS3(5)+0.05*cor2*(1.-2./(1.+cc42**2))
       BbE1(1,2)=BbE1(1,2)+0.03*cor2*(1.-2./(1.+cc43**2))
       BbS1(1,2)=BbS1(1,2)+0.03*cor2*(1.-2./(1.+cc44**2))
       cor1=35.*(3./Q2)*qpi/(1.+100.*qpi**2)
       BbM3(8)=BbM3(8)+0.06*cor1*(1.-2./(1.+cc45**2))
       BbS3(8)=BbS3(8)+0.03*cor1*(1.-2./(1.+cc46**2))
       BbM1(1,1)=BbM1(1,1)+0.1*cor1*(1.-2./(1.+cc47**2))
       BbM1(1,7)=BbM1(1,7)+0.05*cor1*(1.-2./(1.+cc48**2))
       BbS1(1,1)=BbS1(1,1)+0.05*cor1*(1.-2./(1.+cc49**2))
c
c- Background multipole amplitudes(3/2) for resonances with isospin 3/2.
c- resonances (1),(4),(10)  => (1+)
c       BM3(1),BE3(1),BS3(1)
c- resonances (2),(12)      => (2-)
c       BM3(2),BE3(2),BS3(2)
c- resonances  (3)          => (3+)
c       BM3(3),BE3(3),BS3(3)
c- resonances  (5),(9)      => (0+)
c       BE3(5),BS3(5)
c- resonances  (6),(13)     => (3-)
c       BM3(6),BE3(6),BS3(6)
c- resonances  (7)          => (2+)
c       BM3(7),BE3(7),BS3(7)
c- resonances  (8),(11)     => (1-)
c       BM3(8),BE3(8),BS3(8)
C-constract background multipole amplitudes with I=1/2
c- resonances  (1),(10)          => (1-)
c        BM1(i,1),BS1(i,1)
c- resonances  (2),(5)           => (0+)
c        BE1(i,2),BS1(i,2)
c- resonances  (3),(4),(12)      => (2-)
c       BM1(i,3),BE1(i,3),BS1(i,3) 
c- resonances  (6),(13)          => (3-)
c       BM1(i,6),BE1(i,6),BS1(i,6) 
c- resonances  (7),(11)          => (1+)
c       BM1(i,7),BE1(i,7),BS1(i,7) 
c- resonances  (8)               => (3+)
c       BM1(i,8),BE1(i,8),BS1(i,8) 
c- resonances  (9)               => (2+)
c       BM1(i,9),BE1(i,9),BS1(i,9) 
c	write(6,'(''backg 2'',10f7.3)') 
c     >      BM3,BE3,BS3,
c     +      BM1,BE1,BS1,
c     +      BbM3,BbE3,BbS3,
c     +      BbM1,BbE1,BbS1,
c     +      RM3,RE3,RS3,
c     +      RM1,RE1,RS1,
c     +      E0,M1,S0,S1
      RETURN
      END
c
c------------------------------------------------------------------------------

      SUBROUTINE MIXING_AMPL (S,Q2)
c     
      IMPLICIT NONE
      SAVE
c     
#include "gpar.inc"
#include "fit_par.inc"
#include "back_ampl.inc"

      REAL S, Q2, kgamma, qpion
      INTEGER i, ll, j,k
      REAL w,on1,on2,op1,om1,op2,om2 
      REAL B1(2),Lpv,Lmix 
c     -  here we calculate multipole amplitudes for PV-PS mixing 
        W=sqrt(s)
      on1=W-Egamma
      on2=W-Epion
      op1=on1+mn 
      om1=on1-mn
      op2=on2+mn
      om2=on2-mn  
	Lmix=0.5
        Lpv=Lmix**2/(Lmix**2+Epion**2-mp22)
        B1(1)=Lpv*3.134*(F2p-F2n)/mn/W/2.
        B1(2)=Lpv*3.134*(F2p+F2n)/mn/W/2.
C -this is loop for (+) and (0) amplitudes. (-)=0.
		DO i=1,2
	E0(i)=sqrt(op1*op2)*(W-mn)*B1(i)
	M1(i)=-sqrt(om1*om2)*(W+mn)*B1(i)
	S0(i)=sqrt(om1*op2)*op1*B1(i)
	S1(i)=-sqrt(op1*om2)*om1*B1(i)
		ENDDO

      RETURN
      END
c
c------------------------------------------------------------------------------
      
      SUBROUTINE HIGH_MULT (S,Q2,costh)
c     
      IMPLICIT NONE
      SAVE

#include "multipoles.inc"
#include "highmult.inc"
#include "amplitudes.inc"
#include "fit_par.inc"

      REAL S, Q2, costh
      INTEGER NWRITE
      DATA NWRITE/0/
      REAL p1,p2,p3,p4,p5,p6,p7,p8
      REAL qpion,mnuc,mpion
      INTEGER i,j
    
      CALL back_high(S,Q2)
         mnuc=0.93827
         mpion=0.135
      qpion=sqrt(((s-mnuc**2+mpion**2)/2.)**2/s-mpion**2)
c-P''_5
	p1=-52.5*costh+157.5*costh**3
c-P''_3
	p2=15.*costh
c-P''_4
	p3=-7.5+52.5*costh**2
c-P'_5
	p4=1.875-26.25*costh**2+39.375*costh**4
c-P'_3
	p5=-1.5+7.5*costh**2
c-P'_4
	p6=-7.5*costh+17.5*costh**3

c-P''_6
	p7=13.125-236.25*costh**2+433.125*costh**4
c-P'_6
	p8=13.125*costh-78.75*costh**3+86.625*costh**5
c-These are loops for p(1/2), n(1/2), 3/2
	  do i=1,3
	phih(1,i)=(4.*M4p(i)+E4p(i))*p4+(5.*M5p(i)+E5p(i))*p8
	phih(1,i)=phih(1,i)+(5.*M4m(i)+E4m(i))*p5
	phih(1,i)=phih(1,i)+(6.*M5m(i)+E5m(i))*p6
	phih(2,i)=(5.*M4p(i)+4.*M4m(i))*p6+(6.*M5p(i)+5.*M5m(i))*p4
	phih(3,i)=(E4p(i)-M4p(i))*p1+(E5p(i)-M5p(i))*p7
	phih(3,i)=phih(3,i)+(M4m(i)+E4m(i))*p2+(M5m(i)+E5m(i))*p3
	phih(4,i)=(M4p(i)-E4p(i)-M4m(i)-E4m(i))*p3
	phih(4,i)=phih(4,i)+(M5p(i)-E5p(i)-M5m(i)-E5m(i))*p1
	phih(5,i)=5.*S4p(i)*p4-4.*S4m(i)*p5+6.*S5p(i)*p8-5.*S5m(i)*p6
	phih(6,i)=-5.*S4p(i)*p6+4.*S4m(i)*p6-6.*S5p(i)*p4+5.*S5m(i)*p4
	  enddo
c Construct phih for reactions pi0,pi+
	do i=1,6
	 phih(i,1)=phih(i,1)+phih(i,3)/1.5
	 phih(i,2)=sqrt(2.)*(phih(i,1)-phih(i,3))
	enddo
      RETURN
      END

      SUBROUTINE BACK_HIGH(S,Q2)

      IMPLICIT NONE
      SAVE

#include "gpar.inc"
#include "multampl.inc"
#include "fit_par.inc"
#include "highmult.inc"

      REAL S, Q2
      INTEGER i,j,kpower,iw,jw,minp,maxp,mwp
      data kpower/8/
      REAL wpoint,tmult(8),wx(8),a_point
      REAL bem(2,3),bep(2,3),bmm(2,3),bmp(2,3),bsm(2,3),bsp(2,3)

      wpoint=sqrt(s)
      mwp=1
      do iw=1,maxmpoints-1
         if(wpoint.lt.wtab(iw).and.mwp.eq.1)then
            mwp=iw
            go to 101
         endif
      enddo
 101  minp=mwp-3
      if(minp.le.0)then
         minp=1
         maxp=8
      else
         maxp=mwp+4
         if(maxp.gt.maxmpoints)then
            maxp=maxmpoints
            minp=maxmpoints-7
         endif
      endif
c
      do i=4,5
         do j=1,3
c E"l"m
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=aem(iw,i+1,j)
               wx(jw)=wtab(iw)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bem(i-3,j)=a_point
c E"l"p
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=aep(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bep(i-3,j)=a_point
c
c M"l"m
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=amm(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bmm(i-3,j)=a_point
c M"l"p
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=amp(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bmp(i-3,j)=a_point
c S"l"m
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=asm(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bsm(i-3,j)=a_point
c S"l"p
            jw=0
            do iw=minp,maxp
               jw=jw+1
               tmult(jw)=asp(iw,i+1,j)
            enddo
            call polint(tmult,wx,kpower,wpoint,a_point)
            bsp(i-3,j)=a_point
c
         enddo
      enddo
C-construct background multipole amplitudes with I=1/2
C-for combinations p(1/2) and n(1/2)
	E4m(1)=bem(1,2)/3.+bem(1,1)
	E4m(2)=-bem(1,2)/3.+bem(1,1)
	E4p(1)=bep(1,2)/3.+bep(1,1)
	E4p(2)=-bep(1,2)/3.+bep(1,1)
	M4m(1)=bmm(1,2)/3.+bmm(1,1)
	M4m(2)=-bmm(1,2)/3.+bmm(1,1)
	M4p(1)=bmp(1,2)/3.+bmp(1,1)
	M4p(2)=-bmp(1,2)/3.+bmp(1,1)
	S4m(1)=bsm(1,2)/3.+bsm(1,1)
	S4m(2)=-bsm(1,2)/3.+bsm(1,1)
	S4p(1)=bsp(1,2)/3.+bsp(1,1)
	S4p(2)=-bsp(1,2)/3.+bsp(1,1)
	E5m(1)=bem(2,2)/3.+bem(2,1)
	E5m(2)=-bem(2,2)/3.+bem(2,1)
	E5p(1)=bep(2,2)/3.+bep(2,1)
	E5p(2)=-bep(2,2)/3.+bep(2,1)
	M5m(1)=bmm(2,2)/3.+bmm(2,1)
	M5m(2)=-bmm(2,2)/3.+bmm(2,1)
	M5p(1)=bmp(2,2)/3.+bmp(2,1)
	M5p(2)=-bmp(2,2)/3.+bmp(2,1)
	S5m(1)=bsm(2,2)/3.+bsm(2,1)
	S5m(2)=-bsm(2,2)/3.+bsm(2,1)
	S5p(1)=bsp(2,2)/3.+bsp(2,1)
	S5p(2)=-bsp(2,2)/3.+bsp(2,1)
C-constract background multipole amplitudes with I=3/2
	E4m(3)=bem(1,3)
	E5m(3)=bem(2,3)
	E4p(3)=bep(1,3)
	E5p(3)=bep(2,3)
	M4m(3)=bmm(1,3)
	M5m(3)=bmm(2,3)
	M4p(3)=bmp(1,3)
	M5p(3)=bmp(2,3)
	S4m(3)=bsm(1,3)
	S5m(3)=bsm(2,3)
	S4p(3)=bsp(1,3)
	S5p(3)=bsp(2,3)
c
      RETURN
      END
c
c------------------------------------------------------------------------------
        SUBROUTINE INTERP(X,Y,NPTS,NTERMS,XIN,YOUT)
C
C       interpolate between data points to evaluate a function
C       copy form R. Bevington, p.266
C
C
C       THIS PROGRAM INTERPOLATES BETWEEN DATA POINTS.
C       (X,Y) ARE DATA POINTS,
C       NPTS= NUMBER OF SUCH DATA POINT PAIRS,
C       NTERMS= ORDER OF THE POLYNOMIAL USED + 1,
C       XIN= INPUT X,
C       YOUT= INTERPOLATED VALUE.
C       WE ARE OK FOR UP TO 10 TERMS.
C
        implicit none
        INTEGER I1,I2,K,J,NPTS,NTERMS,i,ix,imax,ixmax
        REAL DELTAX,SUM,PROD,DENOM
        REAL DELTA(10),A(10)
        REAL  X(51),Y(51),XIN,YOUT
        DO 19 I=1,NPTS
        IF (XIN-X(I)) 13,17,19
13      I1=I-NTERMS/2
        IF (I1) 15,15,21
15      I1=1
        GO TO 21
17      YOUT=Y(I)
        GO TO 61
19      CONTINUE
        I1=NPTS-NTERMS+1
21      I2=I1+NTERMS-1
        IF (NPTS-I2) 23,31,31
23      I2=NPTS
        I1=I2-NTERMS+1
        IF (I1) 26,26,31
26      I1=1
        NTERMS=I2-I1+1
31      DENOM=X(I1+1)-X(I1)
        DELTAX=(XIN-X(I1))/DENOM
        DO 35 I=1,NTERMS
        IX=I1+I-1
35      DELTA(I)=(X(IX)-X(I1))/DENOM
        A(1)=Y(I1)
        DO 50 K=2,NTERMS
        PROD=1.
        SUM=0.
        IMAX=K-1
        IXMAX=I1+IMAX
        DO 49 I=1,IMAX
           J=K-I
           PROD=PROD*(DELTA(K)-DELTA(J))
 49        SUM=SUM-A(J)/PROD
 50        A(K)=SUM+Y(IXMAX)/PROD
           SUM=A(1)
           DO 57 J=2,NTERMS
              PROD=1.
              IMAX=J-1
              DO 56 I=1,IMAX
 56              PROD=PROD*(DELTAX-DELTA(I))
 57              SUM=SUM+A(J)*PROD
                 YOUT=SUM
 61              CONTINUE
        RETURN
        END     




      SUBROUTINE HEL_AMPL(S,Q2,THETA,kgamma)
      IMPLICIT NONE
      SAVE
#include "gpar.inc"
#include "inpar.inc"
#include "phiampl.inc"
#include "helicity.inc"

      REAL S,Q2,theta,kgamma
      CHARACTER*10 opt
      INTEGER NWRITE
      INTEGER i,j
      REAL teta,qq,sx,cx,sx2
      REAL PHreac(6,3,2)

      teta=theta*pi/180.
      qq=sqrt(Q2)/kgamma
      sx=sin(teta)
      cx=cos(teta)
      sx2=sx**2
c- Calculate Intermediate amplitudes for different reactions. 
c
c- Loop by "i" counts for amplitudes. Loop by "j" stands for 
c- real(j=1) and imaginary(j=2)  parts.
      do i=1,6
         do j=1,2
c  - Reaction g + p => p + pi0
            PHreac(i,1,j)=phi_amp(i,j,1)
c  - Reaction g + p => n + pi+
            PHreac(i,2,j)=phi_amp(i,j,2) 
         end do
      end do	
c
c- Now get the helicity amplitudes.
c- Here as before "j" stands for real and imaginary parts.
c- Loop for "i" done for the reactions:
c- i=1 - g + p => p + pi0
c- i=2 - g + p => n + pi+
c- h1,h2,...h6 are helicity amplitudes f_m2_m1_l,   
c- where m1 and m2 are helicities of initial and final nucleons, l is 
c- helicity of incomming photon. 
      do i=1,2
               do j=1,2
c - h1 = f++- = f--+	
            h1(i,j)=-sx*(PHreac(3,i,j)+cx*PHreac(4,i,j))/a2
c - h2 = f-++ = -f+--	
	h2(i,j)=2.*PHreac(1,i,j)-2.*cx*PHreac(2,i,j)+sx2*PHreac(4,i,j)
	h2(i,j)=-h2(i,j)/a2
c- h3 = f+-+ = f-+-	
            h3(i,j)=-sx2*PHreac(4,i,j)/a2
c- h4 = f+++ = f---	
	h4(i,j)=2.*PHreac(2,i,j)+PHreac(3,i,j)+cx*PHreac(4,i,j)
	h4(i,j)=sx*h4(i,j)/a2
c- h5 = f--0 = -f++0	
            h5(i,j)=qq*(PHreac(5,i,j)+cx*PHreac(6,i,j))
c- h6 = f+-0 = f-+0	
            h6(i,j)=qq*sx*PHreac(6,i,j)
         end do
      end do
      RETURN
      END
c
c------------------------------------------------------------------------------

      SUBROUTINE response_functions (s)
      
      IMPLICIT NONE
      SAVE
      
#include "helicity.inc"
#include "res_func.inc"

#include "janrtab.inc"

      REAL s,mpion(2),mnuc,qpi,qkratio,sinth
      INTEGER NWRITE
      INTEGER i,j
      REAL w
      real eps,csf,cs2f,snf,sn2f,sige,sigz,sigez
      real p_target/1./, p_beam/1./
      real azimut/0./, epsilon/0.9/

      mnuc=0.93827
      mpion(1)=0.13498
      mpion(2)=0.13957
      
      DO i=1,2
      
        qpi = sqrt(((s-mnuc**2+mpion(i)**2)/2.)**2/s-mpion(i)**2)
        qkratio = qpi/((s-mnuc**2)/2./sqrt(s))  
         
        sigmaT(i)=0.
        sigmaL(i)=0.
        sigmaTT(i)=0.
        sigmaTL(i)=0.
        sigmaTLP(i)=0.
        sigmaTy(i)=0.
        sigmaLy(i)=0.
        sigmaTTx(i)=0.
        sigmaTTy(i)=0.
        sigmaTTz(i)=0.
        sigmaTLx(i)=0.
        sigmaTLy(i)=0.
        sigmaTLz(i)=0.
        sigmaTTPx(i)=0.
        sigmaTTPz(i)=0.
        sigmaTLPx(i)=0.
        sigmaTLPy(i)=0.
        sigmaTLPz(i)=0.
        
        Do j=1,2
c j=1 for real parts
c j=2 for imaginary parts
          sigmaT(i) = sigmaT(i) + 
     +     H1(i,j)**2.+ H2(i,j)**2.+ H3(i,j)**2.+ H4(i,j)**2.
          sigmaL(i) = sigmaL(i) + H5(i,j)**2. + H6(i,j)**2.
          sigmaTT(i) = sigmaTT(i)-H4(i,j)*H1(i,j)+H3(i,j)*H2(i,j)
          sigmaTL(i) = sigmaTL(i)+H5(i,j)*(H1(i,j)-H4(i,j))+ 
     -     H6(i,j) * ( H2(i,j) + H3(i,j) )
          sigmaTLPx(i)= sigmaTLPx(i)+
     1     H5(i,j)*(-H3(i,j)+H2(i,j))+H6(i,j)*(H1(i,j)+H4(i,j))
          sigmaTLPy(i)= sigmaTLPy(i)-
     1     H5(i,j)*( H3(i,j)+H2(i,j))+H6(i,j)*(H1(i,j)-H4(i,j))
          sigmaTLPz(i)= sigmaTLPz(i)+
     1     H5(i,j)*( H1(i,j)+H4(i,j))+H6(i,j)*(H3(i,j)-H2(i,j))
          sigmaTTPx(i)= sigmaTTPx(i)+
     1     H1(i,j)*H2(i,j)+H3(i,j)*H4(i,j)
          sigmaTTPz(i)= sigmaTTPz(i)+
     1     H1(i,j)**2.- H2(i,j)**2.+ H3(i,j)**2.- H4(i,j)**2.
        EndDo
        
        sigmaTLP(i) = sigmaTLP(i)+
     1   H5(i,2)*(H1(i,1)-H4(i,1))- H5(i,1)*(H1(i,2)-H4(i,2))+
     1   H6(i,2)*(H2(i,1)+H3(i,1))- H6(i,1)*(H2(i,2)+H3(i,2))
        sigmaTy(i)  = sigmaTy(i)
     1    +H1(i,2)*H2(i,1)-H1(i,1)*H2(i,2)
     1    +H3(i,2)*H4(i,1)-H3(i,1)*H4(i,2)
        sigmaLy(i)  = sigmaLy(i)
     1    +H5(i,2)*H6(i,1)-H5(i,1)*H6(i,2)
        sigmaTLx(i) = sigmaTLx(i)+
     1   H5(i,1)*(H3(i,2)-H2(i,2))-H5(i,2)*(H3(i,1)-H2(i,1))
     1   -H6(i,1)*(H1(i,2)+H4(i,2))+H6(i,2)*(H1(i,1)+H4(i,1))
        sigmaTLy(i) = sigmaTLy(i)
     1   -H5(i,1)*(H3(i,2)+H2(i,2))+H5(i,2)*(H3(i,1)+H2(i,1))
     1   +H6(i,1)*(H1(i,2)-H4(i,2))-H6(i,2)*(H1(i,1)-H4(i,1))
        sigmaTLz(i) = sigmaTLz(i)
     1   -H5(i,1)*(H1(i,2)+H4(i,2))+H5(i,2)*(H1(i,1)+H4(i,1))
     1   +H6(i,1)*(H2(i,2)-H3(i,2))-H6(i,2)*(H2(i,1)-H3(i,1))
        sigmaTTx(i) = sigmaTTx(i)
     1   +H1(i,1)*H3(i,2)-H1(i,2)*H3(i,1)
     1   -H2(i,1)*H4(i,2)+H2(i,2)*H4(i,1)
        sigmaTTy(i) = sigmaTTy(i)
     1   -H1(i,1)*H3(i,2)+H1(i,2)*H3(i,1)
     1   -H2(i,1)*H4(i,2)+H2(i,2)*H4(i,1)
        sigmaTTz(i) = sigmaTTz(i)
     1   -H1(i,1)*H4(i,2)+H1(i,2)*H4(i,1)
     1   -H2(i,1)*H3(i,2)+H2(i,2)*H3(i,1)

c FROM COLE:
c As for the definitions, I use complex variables:
c      sig_t  = (cabs(hh1)**2+cabs(hh2)**2
c     1         +cabs(hh3)**2+cabs(hh4)**2)/2.
c      sig_l  =  cabs(hh5)**2+cabs(hh6)**2
c      sig_tt =  real(-conjg(hh1)*hh4+conjg(hh2)*hh3)ccccccc
c
c      sig_lt  =   real(conjg(hh5)*(hh1-hh4)  +
c     1                 conjg(hh6)*(hh2+hh3))/sqrt(2.0)
c      sig_ltp = -aimag(conjg(hh5)*(hh1-hh4)  +
c     1                 conjg(hh6)*(hh2+hh3))/sqrt(2.0)c
c
c Target polarization (see Knochlein et al., nucl-th/9506029)c
c
c      sigt_t_y   = -aimag(conjg(hh1)*hh2+conjg(hh3)*hh4)
c      sigt_l_y   = -2*aimag(conjg(hh5)*hh6)
c
c      sigt_lt_x  = aimag(conjg(hh5)*(hh3-hh2) -
c     1                   conjg(hh6)*(hh1+hh4))/sqrt(2.)
c      sigt_ltp_x =  real(conjg(hh5)*(hh2-hh3) +
c     1                   conjg(hh6)*(hh1+hh4))/sqrt(2.)
c
c      sigt_lt_y  = aimag(-conjg(hh5)*(hh2+hh3) +
c     1                    conjg(hh6)*(hh1-hh4))/sqrt(2.)
c      sigt_ltp_y =  real(-conjg(hh5)*(hh2+hh3) +
c     1                    conjg(hh6)*(hh1-hh4))/sqrt(2.)
c
c      sigt_lt_z  = aimag(-conjg(hh5)*(hh1+hh4) +
c     1                    conjg(hh6)*(hh2-hh3))/sqrt(2.)
c      sigt_ltp_z = -real(-conjg(hh5)*(hh1+hh4) +
c     1                    conjg(hh6)*(hh2-hh3))/sqrt(2.)
c
c      sigt_tt_x  =  aimag(conjg(hh1)*hh3-conjg(hh2)*hh4)
c      sigt_tt_y  = -aimag(conjg(hh1)*hh3+conjg(hh2)*hh4)
c      sigt_tt_z  = -aimag(conjg(hh1)*hh4+conjg(hh2)*hh3)
c
c      sigt_ttp_x = real(conjg(hh1)*hh2+conjg(hh3)*hh4)
c      sigt_ttp_z = (cabs(hh1)**2-cabs(hh2)**2
c     1             +cabs(hh3)**2-cabs(hh4)**2)/2.
c
c Recoil polarization
c
c      sigr_lt_y = -aimag(-conjg(hh5)*(hh2+hh3) +
c     1                    conjg(hh6)*(hh1-hh4))/sqrt(2.)
c      sigr_ltp_x =  real(conjg(hh5)*(hh3-hh2) +
c     1                   conjg(hh6)*(hh1+hh4))/sqrt(2.)
c      sigr_ttp_z = (cabs(hh1)**2+cabs(hh2)**2
c     1             -cabs(hh3)**2-cabs(hh4)**2)/2.
c
c
c     Calculate direction of target polarization vector with respect
c     to r.h. coordinate system x,y,z with z along q vector, y normal to
c     hadron decay plane.
c
c      pt  = kf_mag*sin(the_e)/sqrt(q2+nu**2)
c      pl  = sqrt(1-pt**2)
c      pn  = -pt * sin(phicm)
c      pt  =  pt * cos(phicm)
c
c      pol_elec = pol_e
c      pol_targ = pol_t
c
c      cph  = cos(phicm)
c      sph  = sin(phicm)
c      c2ph = cos(2*phicm)
c      s2ph = sin(2*phicm)c
c
c      vl   = epsilon
c      vtt  = epsilon
c      vlt  = sqrt(2*epsilon*(1+epsilon))
c      vltp = sqrt(2*epsilon*(1-epsilon))
c      vttp = sqrt(1-epsilon**2)
c
c      ekin = sqrt(Q2)/qv_mag_cm
c
c      sigma_t_y   = fkt  * vt   * pn * sigt_t_y
c      sigma_l_y   = fkt  * vl   * pn * sigt_l_y         *ekin**2
c      sigma_lt_x  = fkt  * vlt  * pt * sigt_lt_x  * sph *ekin
c      sigma_lt_y  = fkt  * vlt  * pn * sigt_lt_y  * cph *ekin
c      sigma_lt_z  = fkt  * vlt  * pl * sigt_lt_z  * sph *ekin
c      sigma_tt_x  = fkt  * vtt  * pt * sigt_tt_x  * s2ph
c      sigma_tt_y  = fkt  * vtt  * pn * sigt_tt_y  * c2ph
c      sigma_tt_z  = fkt  * vtt  * pl * sigt_tt_z  * s2ph
c      sigma_ltp_x = fkt  * vltp * pt * sigt_ltp_x * cph *ekin
c      sigma_ltp_y = fkt  * vltp * pn * sigt_ltp_y * sph *ekin
c      sigma_ltp_z = fkt  * vltp * pl * sigt_ltp_z * cph *ekin
c      sigma_ttp_x = fkt  * vttp * pt * sigt_ttp_x
c      sigma_ttp_z = fkt  * vttp * pl * sigt_ttp_z
c
c      sigpt  = sigma_t_y + sigma_l_y
c    >  sigma_lt_x + sigma_lt_y + sigma_lt_z
c    >  sigma_tt_x + sigma_tt_y + sigma_tt_z
c
c      siget  = sigma_ltp_x+sigma_ltp_y+sigma_ltp_z
c      siget  = siget+sigma_ttp_x+sigma_ttp_z
c
c      sige   =  pol_elec*sigma_ltp
c      sigpt  =  pol_targ*sigpt
c      siget  =  pol_elec*pol_targ*siget
c
c
c END of Cole 
 
      sigmaT(i)    = sigmaT(i)    * qkratio / 2.
      sigmaL(i)    = sigmaL(i)    * qkratio
      sigmaTT(i)   = sigmaTT(i)   * qkratio
      sigmaTL(i)   = sigmaTL(i)   * qkratio / sqrt(2.)
      sigmaTLP(i)  = sigmaTLP(i)  * qkratio / sqrt(2.)
      sigmaTy(i)   = sigmaTy(i)   * qkratio
      sigmaLy(i)   = sigmaLy(i)   * qkratio * 2.
      sigmaTLx(i)  = sigmaTLx(i)  * qkratio / sqrt(2.)
      sigmaTLy(i)  = sigmaTLy(i)  * qkratio / sqrt(2.)
      sigmaTLz(i)  = sigmaTLz(i)  * qkratio / sqrt(2.)
      sigmaTTx(i)  = sigmaTTx(i)  * qkratio 
      sigmaTTy(i)  = sigmaTTy(i)  * qkratio 
      sigmaTTz(i)  = sigmaTTz(i)  * qkratio 
      sigmaTLPx(i) = sigmaTLPx(i) * qkratio / sqrt(2.)
      sigmaTLPy(i) = sigmaTLPy(i) * qkratio / sqrt(2.)
      sigmaTLPz(i) = sigmaTLPz(i) * qkratio / sqrt(2.)
      sigmaTTPx(i) = sigmaTTPx(i) * qkratio 
      sigmaTTPz(i) = sigmaTTPz(i) * qkratio / 2. 

c MAID variables:
      EPS = epsilon
      CSF =  cos(azimut)      
      SNF =  sin(azimut)      
      CS2F =  cos(2.*azimut)      
      SN2F =  sin(2.*azimut)      

c MAID polarizaed beam and target cross section:
c target spin
      sigz = SQRT(2.*EPS*(1.+EPS)) * SNF * sigmaTLZ(i) + 
     >      EPS * SN2F * sigmaTTZ(i)
c beam and target spin
      sigez = SQRT(2.*EPS*(1.-EPS)) * CSF * sigmaTLPZ(i)+
     >  SQRT(1.-EPS**2) * sigmaTTPZ(i)

c beam spin
      sige = SQRT(2.*EPS*(1.-EPS)) * SNF * sigmaTLP(i)

      sigma(i) = sigmaT(i) + 
     >  epsilon * sigmaL(i) + 
     >  epsilon * sigmaTT(i) * cos(2. * azimut) +
     >  sqrt(2.*epsilon*(1.+epsilon))*sigmaTL(i)*cos(azimut )
c pol. beam and target terms
     >  + P_beam * sige 
     >  + P_beam * P_target * sigez 
     >  + p_target * sigz     

      ENDDO
      RETURN
      END


c------------------------------------------------------------------------------

      Subroutine mult_ana(qsq)
   
      IMPLICIT NONE
      SAVE

#include "inpar.inc"
#include "janrtab.inc"
#include "fit_par.inc"
#include "remult.inc"
#include "multampl.inc"
     
      Integer Ierr
   
      INTEGER iw,iq2,icos,iphi,i,j
      INTEGER Nw,Nq2,Ncos,Nphi
      REAL Q2, W, CosTh, Phi, Epsilon, nu, ep, mp, El,qsq
      REAL Fint,Finthigh, gauss, W_min, W_step
      external gauss,fint,finthigh
      integer lm,is,ii,jj,ifile,lun
      character*10 fnamr(24)
      data fnamr/'rm0p.dat','rm0m.dat','re0p.dat','re0m.dat',
     +'rs0p.dat','rs0m.dat',
     +     'rm1p.dat','rm1m.dat','re1p.dat','re1m.dat','rs1p.dat',
     +'rs1m.dat',
     +     'rm2p.dat','rm2m.dat','re2p.dat','re2m.dat','rs2p.dat',
     +'rs2m.dat',
     +     'rm3p.dat','rm3m.dat','re3p.dat','re3m.dat','rs3p.dat',
     +'rs3m.dat'/
      character*10 fnamrr(24)
      character*10 fnamhigh(12)
      data fnamhigh/'rm4p.dat','rm4m.dat','re4p.dat','re4m.dat',
     +'rs4p.dat','rs4m.dat',
     +     'rm5p.dat','rm5m.dat','re5p.dat','re5m.dat','rs5p.dat',
     +'rs5m.dat'/
      common/ifint/lm,is,ii,jj,w,q2
      DATA mp/0.938272/
c     _end_var
c     
c     executable code for routine Mult_ana:
c---- 6----------------------------------------------------------------72
      ierr=0
         cc1=start_value(1)
         cc2=start_value(2)
         cc3=start_value(3)
         cc4=start_value(4)
         cc5=start_value(5)
         cc6=start_value(6)
         cc7=start_value(7)
         cc8=start_value(8)
         cc9=start_value(9)
         cc10=start_value(10)
         cc11=start_value(11)
         cc12=start_value(12)
         cc13=start_value(13)
         cc14=start_value(14)
         cc15=start_value(15)
         cc16=start_value(16)
         cc17=start_value(17)
         cc18=start_value(18)
         cc19=start_value(19)
         cc20=start_value(20)
         cc21=start_value(21)
         cc22=start_value(22)
         cc23=start_value(23)
         cc24=start_value(24)
         cc25=start_value(25)
         cc26=start_value(26)
         cc27=start_value(27)
         cc28=start_value(28)
         cc29=start_value(29)
         cc30=start_value(30)
         cc31=start_value(31)
         cc32=start_value(32)
         cc33=start_value(33)
         cc34=start_value(34)
         cc35=start_value(35)
         cc36=start_value(36)
         cc37=start_value(37)
         cc38=start_value(38)
         cc39=start_value(39)
         cc40=start_value(40)
         cc41=start_value(41)
         cc42=start_value(42)
         cc43=start_value(43)
         cc44=start_value(44)
         cc45=start_value(45)
         cc46=start_value(46)
         cc47=start_value(47)
         cc48=start_value(48)
         cc49=start_value(49)

c       write(6,'(10f7.2)') 
c     >cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10,
c     +cc11,cc12,cc13,cc14,cc15,cc16,cc17,cc18,cc19,cc20,
c     +cc21,cc22,cc23,cc24,cc25,cc26,cc27,cc28,cc29,cc30,
c     +cc31,cc32,cc33,cc34,cc35,cc36,cc37,cc38,cc39,cc40,
c     +cc41,cc42,cc43,cc44,cc45,cc46,cc47,cc48,cc49
C
C     
c
      IF( fitopt.NE.'M' .AND. fitopt.NE.'m' ) THEN
         Nw = 101
         W_min=1.1
         W_step=0.01
      ELSE
         If ( Wstep.GT.0 .AND. Wmin+Wstep.LT.Wmax )Then
            Nw = INT( (Wmax - Wmin) / Wstep ) + 1
            W_min=Wmin
            W_step=Wstep
         Else
            if( Wmin.Eq.Wmax )then
               Nw=1
               W_min=Wmin
               W_step=Wstep
            endif
         EndIf
      ENDIF
      write(6,*)'Number of W points = ',nw
c
         do jj=1,1
      do lm=0,3
         do is=1,6
            lun=(lm+2)*10+is
            ifile=lm*6+is
            open(unit=lun,file=fnamr(ifile),status='UNKNOWN')
            write(lun,*)fnamr(ifile)
            write(lun,*)
         enddo
      enddo   
c     
c change to pass as argument
c      Q2 = (Q2max+Q2min)/2. 
      q2 = qsq

      Do iw=1,Nw
c     Loop over W
         W = W_min + float(iw-1)*W_step
         el=(w**2-mp**2)/2./mp
         do lm=0,3
            do is=1,3
               do ii=1,6
                  Rint(lm,is,ii,jj)=gauss(Fint,-1.,1.,IntAccur)
               enddo
            enddo
            do is=1,3
               call re_mult()
            enddo
            lun=(lm+2)*10
            write(lun+1,201)W,El,(Mlp(lm,is,jj),is=1,3)
            write(lun+2,201)W,El,(Mlm(lm,is,jj),is=1,3)
            write(lun+3,201)W,El,(Elp(lm,is,jj),is=1,3)
            write(lun+4,201)W,El,(Elm(lm,is,jj),is=1,3)
            write(lun+5,201)W,El,(Slp(lm,is,jj),is=1,3)
            write(lun+6,201)W,El,(Slm(lm,is,jj),is=1,3)
         enddo
      EndDo
      do lm=0,3
         do is=1,6
            lun=(lm+2)*10+is
            close(lun)
         enddo
      enddo   
      do lm=4,5
         do is=1,6
            lun=(lm+2)*10+is
            ifile=(lm-4)*6+is
            open(unit=lun,file=fnamhigh(ifile),status='UNKNOWN')
            write(lun,*)fnamhigh(ifile)
            write(lun,*)
         enddo
      enddo   
c     
c change to pass as argument
c      Q2 = (Q2max+Q2min)/2. 
      q2 = qsq

      Do iw=1,Nw
c     Loop over W
         W = W_min + float(iw-1)*W_step
         el=(w**2-mp**2)/2./mp
         do lm=4,5
            do is=1,3
               do ii=1,6
                  Rint(lm,is,ii,jj)=gauss(Finthigh,-1.,1.,IntAccur)
               enddo
            enddo
            do is=1,3
               call re_mult()
            enddo
            lun=(lm+2)*10
            write(lun+1,201)W,El,(Mlp(lm,is,jj),is=1,3)
            write(lun+2,201)W,El,(Mlm(lm,is,jj),is=1,3)
            write(lun+3,201)W,El,(Elp(lm,is,jj),is=1,3)
            write(lun+4,201)W,El,(Elm(lm,is,jj),is=1,3)
            write(lun+5,201)W,El,(Slp(lm,is,jj),is=1,3)
            write(lun+6,201)W,El,(Slm(lm,is,jj),is=1,3)
         enddo
      EndDo
      do lm=4,5
         do is=1,6
            lun=(lm+2)*10+is
            close(lun)
         enddo
      enddo   
      enddo   
 201  format(2x,2F7.3,2x,3F12.5)
c
      return
      end
      
      Subroutine Re_mult()

      IMPLICIT NONE
      SAVE

#include "remult.inc"
      
      common/ifint/lm,is,ii,jj,w,q2
      integer lm,is,ii,jj
      real w,q2 
c ii, W and Q2 are not used here 
c
      if(lm.ge.0)then
         Elp(lm,is,jj)=(-2.*lm*Rint(lm,is,1,jj)+Rint(lm,is,2,jj)+
     +2.*lm*Rint(lm,is,3,jj)-
     -lm*(lm-1)*Rint(lm,is,4,jj))/(2.*(lm+1.)*(2.*lm+1.))
         Slp(lm,is,jj)=(Rint(lm,is,5,jj)-
     -lm*Rint(lm,is,6,jj))/(lm+1.)/(2.*lm+1.)
         Mlp(lm,is,jj)=0.
         Mlm(lm,is,jj)=0.
         Elm(lm,is,jj)=0.
         Slm(lm,is,jj)=0.
      endif
      if(lm.ge.1)then
         Mlp(lm,is,jj)=(2.*Rint(lm,is,1,jj)+
     +Rint(lm,is,2,jj)+2.*lm*Rint(lm,is,3,jj)+
     +(lm-1.)*(lm+2.)*Rint(lm,is,4,jj))/(2.*(lm+1.)*(2.*lm+1.))
         Mlm(lm,is,jj)=(-2.*Rint(lm,is,1,jj)-
     - Rint(lm,is,2,jj)+2.*(lm+1)*Rint(lm,is,3,jj)-
     -(lm-1.)*(lm+2.)*Rint(lm,is,4,jj))/(2.*lm*(2.*lm+1.))
         Slm(lm,is,jj)=(Rint(lm,is,5,jj)+
     +(lm+1.)*Rint(lm,is,6,jj))/(lm*(2.*lm+1.))
      endif
      if(lm.ge.2)then
	Elm(lm,is,jj)=(2.*(lm+1.)*Rint(lm,is,1,jj)+Rint(lm,is,2,jj)-
     -2.*(lm+1.)*Rint(lm,is,3,jj)-
     -(lm+1.)*(lm+2.)*Rint(lm,is,4,jj))/(2.*lm*(2.*lm+1.))
      endif

      return
      end


      Function Finthigh(x)

      IMPLICIT NONE
      SAVE
      
#include "rephi.inc"
      
      common/ifint/lm,is,i,j,w,q2

      Real Finthigh,x
      Real W,Q2
      Integer i,j,is,lm,ii,ir
      real p0,p1,p2,p3,p4,p5,p6,rim(6,3),PHis(6,3)
      real pim(6,2),pip(6,5,3),PHisotop(6,3)
      real extB,extR
      integer ierr

      call janr_rephi(W,Q2,x,ierr)
	extB=1.
	extR=0.
	do ii=1,6
	do ir=1,3
	 Phisotop(ii,ir)=Phisot(1,ii,ir)*extB+Phisot(2,ii,ir)*extR
	enddo
	enddo
c
c- is=1 for (3/2), is=2 for (1/2), is=3 for (0)
c
c- (0)
      if(is.eq.3)then
         do ii=1,6
            PHis(ii,is)=Phisotop(ii,3)
         enddo
c- (1/2)
      elseif(is.eq.2)then
         do ii=1,6
            PHis(ii,is)=Phisotop(ii,1)+2.*Phisotop(ii,2)
         enddo
c- (3/2)
      elseif(is.eq.1)then
         do ii=1,6
            PHis(ii,is)=Phisotop(ii,1)-Phisotop(ii,2)
         enddo
      endif
c
	p0=1.-x**2
c-P_4
	p1=(3.-30.*x**2+35.*x**4)/8.
c-P_5
	p2=(15.*x-70.*x**3+63.*x**5)/8.
c-P'_4
        p3=-7.5*x+17.5*x**3
c-P'_5
        p4=1.875-26.25*x**2+39.375*x**4
c-P''_4
        p5=-7.5+52.5*x**2  
c-P''_5
        p6=-52.5*x+157.5*x**3

      rim(1,is)=PHis(3,is)+x*PHis(4,is)
      rim(2,is)=2.*PHis(1,is)-2.*x*PHis(2,is)+p0*PHis(4,is)
      rim(3,is)=PHis(2,is)
      rim(4,is)=PHis(4,is)
      rim(5,is)=PHis(5,is)+x*PHis(6,is)
      rim(6,is)=PHis(6,is)
c
	pim(1,1)=-9.*p0*p3/45.
	pim(1,2)=-11.*p0*p4/60.
	pim(2,1)=4.5*p1
	pim(2,2)=5.5*p2
	pim(3,1)=-pim(1,1)
	pim(3,2)=-pim(1,2)
	pim(4,1)=p0**2*p5/80.
	pim(4,2)=11.*p0**2*p6/1680.
	pim(5,1)=pim(2,1)
	pim(5,2)=pim(2,2)
	pim(6,1)=pim(3,1)
	pim(6,2)=pim(3,2)
c
      if(lm.eq.5)then
         do ii=1,6
		pip(ii,lm,is)=pim(ii,2)*rim(ii,is)
         enddo
c
      elseif(lm.eq.4)then  
         do ii=1,6   
		pip(ii,lm,is)=pim(ii,1)*rim(ii,is)
         enddo
      endif  
c
      finthigh=0.0 
      goto(11,12,13,14,15,16)i
c
 11   continue
c - Integeral 1 (I1)   
      finthigh=pip(1,lm,is)
      go to 113
c
 12   continue
c - Integeral 2 (I2)
      finthigh=pip(2,lm,is)
      go to 113
c
 13   continue
c - Integeral 3 (I3)   
      finthigh=pip(3,lm,is)
      go to 113
c
 14   continue
c - Integeral 4 (I4)
      finthigh=pip(4,lm,is)
      go to 113
c
 15   continue
c - Integeral 5 (I5)
      finthigh=pip(5,lm,is)
      go to 113
c
 16   continue
c - Integeral 6 (I6)
      finthigh=pip(6,lm,is)
c
113   continue
      return
      end

      SUBROUTINE JANR_RePhi(Wi,Q2i,costhi,ierr)
      
      IMPLICIT NONE
      SAVE

#include "gpar.inc"
#include "rephi.inc"

      REAL Wi, Q2i, epsilon, costhi, phidi
      INTEGER ierr
      INTEGER i,j
      REAL sinth,theta,im_phi_save(6,3),phi
      REAL kgamma,qpion,k22,q22
      REAL W, costh, phid
      REAL Q2, S
      
      w=wi
      q2=q2i
      costh=costhi
      phid=phidi

      CALL JANR_ini_point(Q2)
c      
      theta=acos(costh)
      theta=theta*180./pi
      sinth=sqrt(1.-costh**2)
      phi=phidi*pi/180.
c
      s=w**2
      w2=2.*w
      E=s-mn22
      Egamma=(E-Q2)/w2
      Epion=(E+mp22)/w2
      q22=Epion**2-mp22
      k22=(w-Egamma)**2-mn22
      qpion=sqrt(q22)
      kgamma=sqrt(k22)
      qk=qpion*kgamma*costh-Epion*Egamma
      t=mp22-Q2+2.*qk
      u=2.*mn22+mp22-Q2-s-t
      CALL etta
      CALL born(S,Q2)
      CALL RE_PHI(S,Q2,COSTH)
      ierr = 0
      RETURN
      END

c------------------------------------------------------------------------------
      SUBROUTINE ETTA
      IMPLICIT NONE
      SAVE
      
#include "ettapar.inc"

      INTEGER ETampl(6)
      INTEGER ETisot(3)
      INTEGER i,k,m 

c
c      IF (NWRITE .LT. 1) THEN
c        NWRITE = NWRITE + 1
c        CRMESS='This is a DUMMY routine, this message written once'
c        CALL RECMES(CRNAME,'I',CRMESS)
c      ENDIF
c
c- Definition of coefficents that denote the symmetry of Invariant 
c- Amplitudes under the interchange "s<==>u". 
c
      ETampl(1) = 1
      ETampl(2) = 1
      ETampl(6) = 1
      do i=3,5
         ETampl(i) = -1
      end do
c
c- Coefficients which determine the symmetry of Isospin Amplitudes 
c- (+),(-),(0) under the interchange "s<==>u".
c
      ETisot(1) = 1
      ETisot(2) = -1
      ETisot(3) = 1
c
c Products of two symmetry coefficients.
c
      DO k=1,6
        do m=1,3
           ett(k,m) = ETampl(k)*ETisot(m)
        enddo
      ENDDO
      RETURN
      END

c------------------------------------------------------------------------------

      SUBROUTINE BORN(S,Q2)

      IMPLICIT NONE
      SAVE

#include "ettapar.inc"
#include "gpar.inc"
#include "bornterms.inc"
#include "fit_par.inc"
      
      REAL S,Q2
      REAL GAMMF
      REAL d,d1,tau,f_1(2),f_2(2)
      REAL q1,z,Gmp,Gep,Gmn,Gen,F1p,F1n
      REAL EF1(3),EF2(3),Re(6,3),Rm(6,3)
      REAL Rro(6,3),Rom(6,3)
      REAL Lambda_pi,Lambda_ro
      REAL Rpi_re(6,3)
      REAL Rro_re(6,3)
      REAL Rom_re(6,3),RRe(6,3),RRm(6,3)
      REAL F_pion, F_dipole,gvec,gro,gro2,gom,gom2
      REAL Lamdom,Lamdro,gvecom,gvecro
      REAL ae1,be1,be2,be3,am1,bm1,bm2,bm3
      REAL ce1,ce2,de1,de2,de3,de4,cm1,dm1,dm2,dm3
      INTEGER i,j,k

c   Calculate Born term contribution (nucleon and pion poles)
c     into invariant amplitudes.
      
	d=1./(s-mn22)
	d1=u-mn22
        tau=Q2/4./mn22
	q1=sqrt(Q2)
c- Pion Formfactor
c- Dipole formula
        F_pion=1./(1.+Q2/0.54)
        F_dipole=1./(1.+Q2/0.71)**2
c-  Parameterizations from NP B (Proc.Suppl.) 159(2006)127
        ae1=-0.0578
        be1=11.1
        be2=13.6
        be3=33.0
        am1=0.15
        bm1=11.1
        bm2=19.6
        bm3=7.54
        ce1=1.25
        ce2=1.3
        de1=-9.86
        de2=305.
        de3=-758.
        de4=802.
        cm1=1.81
        dm1=14.1
        dm2=20.7
        dm3=68.7
        Gep=(1.+ae1*tau)/(1.+be1*tau+be2*tau**2+be3*tau**3)
        Gmp=2.793*(1.+am1*tau)/(1.+bm1*tau+bm2*tau**2+bm3*tau**3)
        Gen=ce1*tau+ce2*tau**2
        Gen=Gen/(1.+de1*tau+de2*tau**2+de3*tau**3+de4*tau**4)
c- this is my parametrization for Gen
	Gen=Gen*(1.+1./(1.+2.*(Q2-0.4)**2))
        Gmn=-1.913*(1.+cm1*tau)/(1.+dm1*tau+dm2*tau**2+dm3*tau**3)
c- Neutron Pauli Formfactors
	F1p=(Gep+tau*Gmp)/(1.+tau)
	F2p=(Gmp-Gep)/(1.+tau)/mn2
	F1n=(Gen+tau*Gmn)/(1.+tau)
	F2n=(Gmn-Gen)/(1.+tau)/mn2
c- Isovecotor and Isoscalar Pauli Formfactors (1 is V, 2 is S). 
	F_1(1)=F1p-F1n
	F_1(2)=F1p+F1n
	F_2(1)=F2p-F2n
	F_2(2)=F2p+F2n
	do k=1,2
c-  This loop counts for isotopic amplitudes (+) and (-),
c-  which are determined by isovector form factors	
	  EF1(k)=F_1(1)
	  EF2(k)=F_2(1)
	end do
c-  The following are isotopic form factors (0) 
c-  which are equal to isoscalar form factors
	EF1(3)=F_1(2)
	EF2(3)=F_2(2)	
	do k=1,3
c-  This loop counts for isotopic states (+),(-),(0)	
	  RRe(1,k)=EF1(k)
	  RRe(2,k)=-EF1(k)
	  RRe(3,k)=-EF1(k)/2.
	  RRe(4,k)=0.
	  RRe(5,k)=0.
	  RRe(6,k)=0.
	  RRm(1,k)=mn2*EF2(k)
	  RRm(2,k)=0.
	  RRm(3,k)=0.
	  RRm(4,k)=EF2(k)
	  RRm(5,k)=-qk*EF2(k)
	  RRm(6,k)=2.*EF2(k)
	end do
	do i=1,6
	  do j=1,3
	    Re(i,j)=3.134*RRe(i,j)*(d+ett(i,j)/d1)
	    Rm(i,j)=3.134*RRm(i,j)*(d+ett(i,j)/d1)
	  end do
	end do
c- The coefficient 3.134 (microbarn)**(1/2) comes from
c- g*e/(8piGeV), where g**2/4pi=13.8, e**2/4pi=1/137.	
	Re(3,2)=Re(3,2)-6.268*f_pion/(t-mp22)
c	The omega exchange contribution: Rw
	gvec=F_dipole
        Lamdom=1.2
        Lamdro=1.5
        gvecom=gvec*Lamdom**2/(Lamdom**2-t)/(t-0.612)
        gvecro=gvec*Lamdro**2/(Lamdro**2-t)/(t-0.591)
	gro=0.75*gvecro*1.1*1.2
	gro2=2.2*gvecro*0.9*1.2
	gom=22.6*gvecom*1.1*1.2
	gom2=-6.44*gvecom*1.2*1.2
        Rro(1,3)=mn2*gro+t*gro2
        Rro(2,3)=0.5*(Q2+mp22-t)*gro2
        Rro(3,3)=(u-s)*gro2/4.
        Rro(4,3)=0.
        Rro(5,3)=0.
        Rro(6,3)=2.*gro
        Rom(1,1)=mn2*gom+t*gom2
        Rom(2,1)=0.5*(Q2+mp22-t)*gom2
        Rom(3,1)=(u-s)*gom2/4.
        Rom(4,1)=0.
        Rom(5,1)=0.
        Rom(6,1)=2.*gom
c
	do i=1,6
          do j=1,2
        k=j+1.
            Rro(i,j)=0.
            Rom(i,k)=0.
          end do
        end do 
c 	
	do i=1,6
	  do j=1,3
	br(1,i,j)=Rom(i,j)+Rro(i,j)+Re(i,j)+Rm(i,j)
	  end do
	end do
      RETURN
      END

      SUBROUTINE RE_PHI(S,Q2,COSTH)

      IMPLICIT NONE
      SAVE

#include "gpar.inc"
#include "bornterms.inc"
#include "rephi.inc"

      REAL S, Q2, COSTH
      INTEGER i,j
      REAL on1,on2,op1,op2,om1,om2,wm,wp,xx
      REAL br1,br2,br3,br4,br5,br6
      REAL AK(6),W

c
c      IF (NWRITE .LT. 1) THEN
c        NWRITE = NWRITE + 1
c        CRMESS='This is a DUMMY routine, this message written once'
c        CALL RECMES(CRNAME,'I',CRMESS)
c      ENDIF
c
c- Determination of Real parts of Intermediate Amplitudes through 
c- Born Invariant Amplitudes  
c
      W = sqrt(S)
      on1=W-Egamma
      on2=W-Epion
      op1=on1+mn
      om1=on1-mn
      op2=on2+mn
      om2=on2-mn
      Wm=W-mn
      Wp=W+mn
      XX=Egamma*qk-Epion*Q2
      AK(1)=sqrt(op1*op2)/w2
      AK(2)=sqrt(om1*om2)/w2
      AK(3)=AK(2)*op2
      AK(4)=AK(1)*om2
      AK(5)=sqrt(om1*op2)/(Q2*w2)
      AK(6)=sqrt(op1*om2)/(Q2*w2)
c     - Loop by Isospin states (+), (-) & (0)
	 do i=1,3
	 do j=1,2
            br5=br(j,5,i)+(E+qk)*br(j,6,i)/2.
c	  
c	  
         PHisot(j,1,i)=AK(1)*(Wm*br(j,1,i)-br5)
c     
         PHisot(j,2,i)=-AK(2)*(Wp*br(j,1,i)+br5)
c     
         br1=2.*br(j,3,i)-br(j,2,i)
         br2=br(j,6,i)/2.-br(j,4,i)
         br3=Q2*br(j,1,i)
c     
         PHisot(j,3,i)=AK(3)*(br1+Wp*br2)
c     
         PHisot(j,4,i)=AK(4)*(-br1+Wm*br2)
c     
         br4=op1*(br3+Wm*br5+W2*om1*(br(j,2,i)-Wp*br(j,6,i)/2.))
         PHisot(j,5,i)=AK(5)*(br4-XX*(br1+Wp*br2))
c     
         br6=-om1*(br3-Wp*br5+W2*op1*(br(j,2,i)+Wm*br(j,6,i)/2.))
c     
         PHisot(j,6,i)=AK(6)*(br6+XX*(br1-Wm*br2))
      end do
      end do
      RETURN
      END

      Function Fint(x)

      IMPLICIT NONE
      SAVE
      
#include "rephi.inc"

      common/ifint/lm,is,i,j,w,q2

      Real Fint,x
      Real W,Q2
      Integer i,j,is,lm,ii
      real p,p1,p2,rim,PHis(6,3,2)
      integer ierr

      call janr_rephi(W,Q2,x,ierr)

c- (0)
      if(is.eq.3)then
         do ii=1,6
            PHis(ii,is,j)=Phisot(j,ii,3)
         enddo
c- (1/2)
      elseif(is.eq.2)then
         do ii=1,6
            PHis(ii,is,j)=Phisot(j,ii,1)+2.*Phisot(j,ii,2)
         enddo
c- (3/2)
      elseif(is.eq.1)then
         do ii=1,6
            PHis(ii,is,j)=Phisot(j,ii,1)-Phisot(j,ii,2)
         enddo
      endif
c
      goto(21,22,23,24)lm+1
 21   continue
      p=1.
      p1=0.
      p2=0.
      go to 111
c
 22   continue
      p=x
      p1=-sqrt(1.-x**2)
      p2=0.
      go to 111
c
 23   continue
      p=3.*x**2/2.-0.5
      p1=-3.*x*sqrt(1-x**2)
      p2=3.*(1.-x**2)
      go to 111
c
 24   continue
      p=(5.*x**2-3.)*x/2.
      p1=-sqrt(1.-x**2)*(15.*x**2-3.)/2.
      p2=15.*x*(1-x**2)
c
 111  continue
      fint=0.0
      goto(11,12,13,14,15,16)i
c
 11   continue
c - Integeral 1 (I1)
      if(lm.lt.1)return
      rim=p1*sqrt(1.-x**2)*(PHis(3,is,j)+x*PHis(4,is,j))
      fint=rim*(2.*lm+1.)/(2.*lm*(lm+1.))
      go to 113
c
 12   continue 
c - Integeral 2 (I2)
      rim=p*(2.*PHis(1,is,j)-2.*x*PHis(2,is,j)+(1.-x**2)*PHis(4,is,j))
      fint=rim*(2.*lm+1.)/2.
      go to 113
c
 13   continue
c - Integeral 3 (I3)
      if(lm.lt.1)return
      rim=p1*sqrt(1.-x**2)*PHis(2,is,j)
      fint=-rim*(2.*lm+1.)/(2.*lm*(lm+1.))
      go to 113
c
 14   continue
c - Integeral 4 (I4)
      if(lm.lt.2)return
      rim=p2*(1.-x**2)*PHis(4,is,j)
      fint=rim*(2.*lm+1.)/(2.*lm*(lm+1.)*(lm-1)*(lm+2))
      go to 113
c
 15   continue 
c - Integeral 5 (I5)
      rim=p*(PHis(5,is,j)+x*PHis(6,is,j))
      fint=rim*(2.*lm+1.)/2.
      go to 113
c
 16   continue
c - Integeral 6 (I6)
      if(lm.lt.1)return
      rim=p1*sqrt(1.-x**2)*PHis(6,is,j)
      fint=-rim*(2.*lm+1.)/(2.*lm*(lm+1.))
c
113   continue
      return
      end

      Subroutine read_multables
c
c - Read pregenerated tables for multipole amplitudes. 
c     March 2000
      IMPLICIT NONE
      SAVE

#include "multampl.inc"

      integer ifile, lun
      character*10 fnamr(24)
      data fnamr/'rm0p.dat','rm0m.dat','re0p.dat','re0m.dat','rs0p.dat',
     +'rs0m.dat',
     +'rm1p.dat','rm1m.dat','re1p.dat','re1m.dat','rs1p.dat','rs1m.dat',
     +'rm2p.dat','rm2m.dat','re2p.dat','re2m.dat','rs2p.dat','rs2m.dat',
     +'rm3p.dat','rm3m.dat','re3p.dat','re3m.dat','rs3p.dat','rs3m.dat'/
c
      real w,e,m23,m12,ms
      integer i,j,iw,ii
c
      do i=0,3
         do j=1,6
            lun=(i+3)*10+j
            ifile=i*6+j
            open(unit=lun,file=fnamr(ifile),status="old")
            do ii=1,2
               read(lun,*)
            enddo
         enddo
c
         lun=(i+3)*10
         do iw=1,maxmpoints
c Read rm"L"p.dat
            read(lun+1,*,end=101)w,e,m23,m12,ms
            amp(iw,i+1,1)=ms
            amp(iw,i+1,2)=m12
            amp(iw,i+1,3)=m23
c Read rm"L"m.dat
            read(lun+2,*,end=101)w,e,m23,m12,ms
            amm(iw,i+1,1)=ms
            amm(iw,i+1,2)=m12
            amm(iw,i+1,3)=m23
c Read re"L"p.dat
            read(lun+3,*,end=101)w,e,m23,m12,ms
            aep(iw,i+1,1)=ms
            aep(iw,i+1,2)=m12
            aep(iw,i+1,3)=m23
c Read re"L"m.dat
            read(lun+4,*,end=101)w,e,m23,m12,ms
            aem(iw,i+1,1)=ms
            aem(iw,i+1,2)=m12
            aem(iw,i+1,3)=m23
c Read rs"L"p.dat
            read(lun+5,*,end=101)w,e,m23,m12,ms
            asp(iw,i+1,1)=ms
            asp(iw,i+1,2)=m12
            asp(iw,i+1,3)=m23
c Read rs"L"m.dat
            read(lun+6,*,end=101)w,e,m23,m12,ms
            asm(iw,i+1,1)=ms
            asm(iw,i+1,2)=m12
            asm(iw,i+1,3)=m23
            wtab(iw)=w
         enddo
 101     continue
         do j=1,6
            lun=(i+3)*10+j
            close(lun)
         enddo
      enddo
c
      return
      end


      Subroutine high_multables

c - Read pregenerated tables for multipole amplitudes. 
c     March 2000

      IMPLICIT NONE
      SAVE

#include "multampl.inc"
      
      integer ifile, lun
      character*10 fnamhigh(12)
c----6----------------------------------------------------------------72
      data fnamhigh/'rm4p.dat','rm4m.dat','re4p.dat','re4m.dat',
     +'rs4p.dat','rs4m.dat','rm5p.dat',
     +'rm5m.dat','re5p.dat','re5m.dat','rs5p.dat','rs5m.dat'/
c
      real w,e,m23,m12,ms
      integer i,j,iw,ii
c
      do i=4,5
         do j=1,6
            lun=(i+3)*10+j
            ifile=(i-4)*6+j
            open(unit=lun,file=fnamhigh(ifile),status="old")
            do ii=1,2
               read(lun,*)
            enddo
         enddo
c
         lun=(i+3)*10
         do iw=1,maxmpoints
c Read rm"L"p.dat
            read(lun+1,*,end=101)w,e,m23,m12,ms
            amp(iw,i+1,1)=ms
            amp(iw,i+1,2)=m12
            amp(iw,i+1,3)=m23
c Read rm"L"m.dat
            read(lun+2,*,end=101)w,e,m23,m12,ms
            amm(iw,i+1,1)=ms
            amm(iw,i+1,2)=m12
            amm(iw,i+1,3)=m23
c Read re"L"p.dat
            read(lun+3,*,end=101)w,e,m23,m12,ms
            aep(iw,i+1,1)=ms
            aep(iw,i+1,2)=m12
            aep(iw,i+1,3)=m23
c Read re"L"m.dat
            read(lun+4,*,end=101)w,e,m23,m12,ms
            aem(iw,i+1,1)=ms
            aem(iw,i+1,2)=m12
            aem(iw,i+1,3)=m23
c Read rs"L"p.dat
            read(lun+5,*,end=101)w,e,m23,m12,ms
            asp(iw,i+1,1)=ms
            asp(iw,i+1,2)=m12
            asp(iw,i+1,3)=m23
c Read rs"L"m.dat
            read(lun+6,*,end=101)w,e,m23,m12,ms
            asm(iw,i+1,1)=ms
            asm(iw,i+1,2)=m12
            asm(iw,i+1,3)=m23
            wtab(iw)=w
         enddo
 101     continue
         do j=1,6
            lun=(i+3)*10+j
            close(lun)
         enddo
      enddo
c
      return
      end
