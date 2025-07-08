      SUBROUTINE RE_PHI(S,Q2,COSTH)
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for subroutine RE_PHI
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
c  Author:         Created:  Mon Dec 16 16:42:54 EST 1996
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
#include "gpar.inc"
#include "bornterms.inc"
#include "rephi.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
      REAL S, Q2, COSTH
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
      PARAMETER (CRNAME='RE_PHI')
      PARAMETER (CRAUTH='')
c
c  Local User defined variables:
c  -----------------------------
c      INTEGER NWRITE
      INTEGER i,j
      REAL on1,on2,op1,op2,om1,om2,wm,wp,xx
      REAL br1,br2,br3,br4,br5,br6
      REAL AK(6),W
c_end_var
c
c  executable code for routine RE_PHI:
c----6----------------------------------------------------------------72
c
c      IF (NWRITE .LT. 1) THEN
c        NWRITE = NWRITE + 1
c        CRMESS='This is a DUMMY routine, this message written once'
c        CALL RECMES(CRNAME,'I',CRMESS)
c      ENDIF
c
c- Determination of Real parts of Intermediate Amplitudes through 
c- Real parts of Invariant Amplitudes determined by Dispersion Relations. 
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
c
c------------------------------------------------------------------------------




