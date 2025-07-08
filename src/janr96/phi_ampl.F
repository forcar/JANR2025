      SUBROUTINE PHI_AMPL (S,costh,Q2,kgamma,qpion)
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for subroutine PHI_AMPL
c
c  Purpose:
c  --------
c
c  Input Parameters:  (Name - Type - Meaning)
c  ----------------
c
c  Output Parameters:  (Name - Type - Meaning)
c  -----------------
c
c  Other routines:
c  ---------------
c
c  Notes:
c  ------
c
c  Author:         Created:  Wed Dec 11 16:01:19 EST 1996
c  -------
c
c  Major revisions:
c  ----------------
c     
c
c_end_doc
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
c BOS common block  uncomment the next line for BOS include file
#include "multipoles.inc"
#include "highmult.inc"
#include "phiampl.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
      REAL S, Q2, kgamma, qpion ,costh
c
c  Local pre-defined variables:
c  ---------------------------
c  RCS information: 
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$RCSfile$')
      PARAMETER (CREVIS = '$Revision$')
      PARAMETER (CSTATE = '$State$')
      PARAMETER (CDATE = '$Date$')
      PARAMETER (CAUTHO = '$Author$')
      DATA CRCSID/   
     1'$Id$'   
     2/   
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
c      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='PHI_AMPL')
      PARAMETER (CRAUTH='I. Aznauryan')
c
c  Local User defined variables:
c  -----------------------------
c      INTEGER NWRITE
      REAL E1(2,2),E3(2)
      REAL p1,p2,p3,p4,p5
      INTEGER i, j
c_end_var
c
c  executable code for routine PHI_AMPL:
c----6----------------------------------------------------------------72
c
c      IF (NWRITE .LT. 1) THEN
c        NWRITE = NWRITE + 1
c        CRMESS='This is a DUMMY routine, this message written once'
c        CALL RECMES(CRNAME,'I',CRMESS)
c      ENDIF
c
c- 
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
c
c- Calculate intermediate amplitudes for transition from 
c- multipole amplitudes to helicity amplitudes and cross sections.
c- Multipole ampl. are normalized according to CGLN Phys.Rev. 106 (1957) 1345.
c- They are equal to mult.ampl. of Devenish&Lyth/8W/pi
c- First do the loop for I=1/2.
c- mult.ampl. (1) => (1-)
c- mult.ampl. (2) => (0+)
c- mult.ampl. (3) => (2-)
c- mult.ampl. (6) => (3-)
c- mult.ampl. (7) => (1+)
c- mult.ampl. (8) => (3+)
c- mult.ampl. (9) => (2+)
c
      do i=1,2
c-construct the real(i=1) and imaginary(i=2) parts of intermediate amplitudes 
      do j=1,2
c-construct the amplitudes for combinations pA(1/2)(j=1), nA(1/2)(j=2)
        E1(i,j)=MM1(i,j,6)+ME1(i,j,6)
c23456789012345678901234567890123456789012345678901234567890123456789012
        phi_amp(1,i,j)=ME1(i,j,2)+p1*(MM1(i,j,7)+ME1(i,j,7))+
     +       p3*(3.*MM1(i,j,8)+ME1(i,j,8))+
     +       3.*MM1(i,j,3)+ME1(i,j,3)+p1*(4.*MM1(i,j,6)+ME1(i,j,6))+
     +       p2*(2.*MM1(i,j,9)+ME1(i,j,9))
        phi_amp(2,i,j)=2.*MM1(i,j,7)+MM1(i,j,1)+
     +       p1*(2.*MM1(i,j,3)+3.*MM1(i,j,9))+
     +       3.*p2*MM1(i,j,6)+4.*p2*MM1(i,j,8)
        phi_amp(3,i,j)=3.*(ME1(i,j,7)-MM1(i,j,7)+E1(i,j))+
     +       p5*(ME1(i,j,8)-MM1(i,j,8))+p4*(ME1(i,j,9)-MM1(i,j,9))
        phi_amp(4,i,j)=-3.*(MM1(i,j,3)+ME1(i,j,3))+
     +       3.*(MM1(i,j,9)-ME1(i,j,9))-
     -       p4*E1(i,j)+p4*(MM1(i,j,8)-ME1(i,j,8))   
        phi_amp(5,i,j)=MS1(i,j,2)+2.*p1*MS1(i,j,7)+
     +       4.*p3*MS1(i,j,8)+3.*p2*MS1(i,j,9)-
     -       2.*MS1(i,j,3)-3.*p1*MS1(i,j,6)
        phi_amp(6,i,j)=MS1(i,j,1)+2.*p1*MS1(i,j,3)+3.*p2*MS1(i,j,6)-
     -       2.*MS1(i,j,7)-4.*p2*MS1(i,j,8)-3.*p1*MS1(i,j,9)
      end do
      end do
c- Now calculate amplitudes for I=3/2.
c- mult.ampl. (1) => (1+)
c- mult.ampl. (2) => (2-)
c- mult.ampl. (3) => (3+)
c- mult.ampl. (5) => (0+)
c- mult.ampl. (6) => (3-)
c- mult.ampl. (7) => (2+)
c- mult.ampl. (8) => (1-)
      do i=1,2
c-construct the real(i=1) and imaginary(i=2) parts of intermediate amplitudes 
c
        E3(i)=MM3(i,6)+ME3(i,6)
c23456789123456789123456789123456789123456789123456789123456789123456789
        phi_amp(1,i,3)=ME3(i,5)+p1*(MM3(i,1)+ME3(i,1))+
     +          p3*(3.*MM3(i,3)+ME3(i,3))+
     +		3.*MM3(i,2)+ME3(i,2)+p1*(4.*MM3(i,6)+ME3(i,6))+
     + 		p2*(2.*MM3(i,7)+ME3(i,7))
        phi_amp(2,i,3)=MM3(i,8)+2.*MM3(i,1)+p1*(2.*MM3(i,2)+3.*MM3(i,7))+
     +		4.*p2*MM3(i,3)+3.*p2*MM3(i,6)
        phi_amp(3,i,3)=3*E3(i)+3.*(ME3(i,1)-MM3(i,1))+
     +          p5*(ME3(i,3)-MM3(i,3))+  
     +		p4*(ME3(i,7)-MM3(i,7))
        phi_amp(4,i,3)=-3.*(MM3(i,2)+ME3(i,2))+3.*(MM3(i,7)-ME3(i,7))-
     -		p4*E3(i)+p4*(MM3(i,3)-ME3(i,3)) 
        phi_amp(5,i,3)=MS3(i,5)+2.*p1*MS3(i,1)+4.*p3*MS3(i,3)+
     +          3.*p2*MS3(i,7)-
     -		2.*MS3(i,2)-3.*p1*MS3(i,6)
        phi_amp(6,i,3)=MS3(i,8)+2.*p1*MS3(i,2)+3.*p2*MS3(i,6)-
     -		2.*MS3(i,1)-4.*p2*MS3(i,3)-3.*p1*MS3(i,7)
c
      end do
	do i=1,6
	do j=1,3
c	phih(i,j)=0.
	phi_amp(i,1,j)=phi_amp(i,1,j)+phih(i,j)
	enddo
	enddo
c
      RETURN
      END
c
c------------------------------------------------------------------------------



