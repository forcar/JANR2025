      SUBROUTINE HIGH_MULT (S,Q2,costh)
c     
      IMPLICIT NONE
      SAVE
c     
c     _begin_inc
c     include files :
c     ---------------------
#include "multipoles.inc"
#include "highmult.inc"
#include "amplitudes.inc"
#include "fit_par.inc"
c     _end_inc
c     
c     _begin_var
c     input/output variables:
c     -----------------------
      REAL S, Q2, costh
c     
c     Local pre-defined variables:
c     ---------------------------
c     RCS information: 
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$RCSfile$')
      PARAMETER (CREVIS = '$Revision$')
      PARAMETER (CSTATE = '$State$')
      PARAMETER (CDATE = '$Date$')
      PARAMETER (CAUTHO = '$Author$')
      DATA CRCSID/   
     1     '$Id$'   
     2     /   
c     Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
c     CHARACTER*100  CRMESS
      PARAMETER (CRNAME='HIGH_MULT')
      PARAMETER (CRAUTH='')
c     
c     Local User defined variables:
c     -----------------------------
      INTEGER NWRITE
      DATA NWRITE/0/
      REAL p1,p2,p3,p4,p5,p6,p7,p8
      REAL M4p(3),E4p(3),S4p(3)
      REAL M4m(3),E4m(3),S4m(3)
      REAL M5p(3),E5p(3),S5p(3)
      REAL M5m(3),E5m(3),S5m(3)
      REAL ww,highmult
      INTEGER i,j
c     _end_var
c     
c     executable code for routine MULT_AMPL:
c---- 6---------------------------------------------------------------
c     
	ww=sqrt(s)
	highmult=0.33*((ww-1.076)/s)**2
	E4m(1)=0.25*highmult
	M4m(1)=0.4*highmult
	E4m(3)=-highmult
	M4m(3)=0.2*highmult
	E4p(1)=0.33*highmult
	M4p(1)=-0.03*highmult
	E4p(3)=-0.55*highmult
	M4p(3)=0.2*highmult
	E5m(1)=0.14*highmult
	M5m(1)=0.14*highmult
	E5m(3)=-0.07*highmult
	M5m(3)=-0.05*highmult
	E5p(1)=0.16*highmult
	M5p(1)=-0.05*highmult
	E5p(3)=-0.25*highmult
	M5p(3)=0.3*highmult
	S4m(1)=-0.14*highmult
	S4m(3)=0.16*highmult
	S4p(1)=0.5*highmult
	S4p(3)=-0.8*highmult
	S5m(1)=-0.14*highmult
	S5m(3)=0.18*highmult
	S5p(1)=0.26*highmult
	S5p(3)=-0.4*highmult

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
	phih(1,i)=phih(1,i)+(5.*M4m(i)+E4m(i))*p5+(6.*M5m(i)+E5m(i))*p6
	phih(2,i)=(5.*M4p(i)+4.*M4m(i))*p6+(6.*M5p(i)+5.*M5m(i))*p4
	phih(3,i)=(E4p(i)-M4p(i))*p1+(E5p(i)-M5p(i))*p7
	phih(3,i)=phih(3,i)+(M4m(i)+E4m(i))*p2+(M5m(i)+E5m(i))*p3
	phih(4,i)=(M4p(i)-E4p(i)-M4m(i)-E4m(i))*p3
	phih(4,i)=phih(4,i)+(M5p(i)-E5p(i)-M5m(i)-E5m(i))*p1
	phih(5,i)=5.*S4p(i)*p4-4.*S4m(i)*p5+6.*S5p(i)*p8-5.*S5m(i)*p6
	phih(6,i)=-5.*S4p(i)*p6+4.*S4m(i)*p6-6.*S5p(i)*p4+5.*S5m(i)*p4
	phih(2,i)=0.
	  enddo

      RETURN
      END
