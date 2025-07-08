      SUBROUTINE CROSS_SECTIONS (azimut,qkratio,epsilon,sinth)
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for subroutine CROSS_SECTIONS
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
c  Author:  H. Aznauryan       Created:  Mon Dec 16 19:31:55 EST 1996
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
#include "helicity.inc"
#include "cross_sec.inc"
#include "fit_par.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
	REAL azimut,qkratio,epsilon,sinth,eps,az
	REAL sigmaTLP(2)
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
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='CROSS_SECTIONS')
      PARAMETER (CRAUTH='')
c
c  Local User defined variables:
c  -----------------------------
      INTEGER NWRITE
      INTEGER i,j
c_end_var
c
c  executable code for routine CROSS_SECTIONS:
c----6----------------------------------------------------------------72
c-loops for reactions p+pi0, n+pi+

      DO i=1,2
        sigmaT(i)=0.
        sigmaL(i)=0.
        sigmaTT(i)=0.
        sigmaTL(i)=0.
        Do j=1,2
          sigmat(i) = sigmat(i) + 
     +                H1(i,j)**2.+ H2(i,j)**2.+ H3(i,j)**2.+ H4(i,j)**2.
          sigmaL(i) = sigmaL(i) + H5(i,j)**2. + H6(i,j)**2.
          sigmaTT(i) = sigmaTT(i)-H4(i,j)*H1(i,j)+H3(i,j)*H2(i,j)
          sigmaTL(i) = sigmaTL(i)+H5(i,j)*(H1(i,j)-H4(i,j))+ 
     -                 H6(i,j) * ( H2(i,j) + H3(i,j) )
        EndDo
      sigmaT(i)=sigmaT(i) * qkratio / 2.
      sigmaL(i)=sigmaL(i) * qkratio
      sigmaTT(i)=sigmaTT(i) * qkratio
      sigmaTL(i)=sigmaTL(i) * qkratio / sqrt(2.)
      
c      sigmaP(i)=sigmaTT(i)/(sinth**2.)
c      sigmaI(i)=sigmaTL(i)/sinth/qkratio
c      sigmaI(i)=sigmaTL(i)/sinth

      sigma(i) = sigmaT(i) + epsilon * sigmaL(i) + 
     + epsilon * sigmaTT(i) * cos(2. * azimut) + sqrt( 2. * epsilon *
     *(1. + epsilon)) * sigmaTL(i) * cos( azimut )

c        sigmaTLP(i)=0.
c          sigmaTLP(i) = sigmaTLP(i)+
c     1        H5(i,2)*(H1(i,1)-H4(i,1))- H5(i,1)*(H1(i,2)-H4(i,2))+
c     1        H6(i,2)*(H2(i,1)+H3(i,1))- H6(i,1)*(H2(i,2)+H3(i,2))
c      sigmaTLP(i)  = sigmaTLP(i)  * qkratio / sqrt(2.)/sinth
c	 sigma(i)= sigmaTLP(i)
      ENDDO
      RETURN
      END
c
c------------------------------------------------------------------------------




