      SUBROUTINE response_functions (azimut,qkratio,epsilon,sinth,w)
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for subroutine response_functions
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
#include "res_func.inc"
#include "janrtab.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
	REAL azimut,qkratio,epsilon,sinth
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
      REAL w
c_end_var
c
c  executable code for routine CROSS_SECTIONS:
c----6----------------------------------------------------------------72
c-loops for reactions p+pi0, n+pi+

      DO i=1,2
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
          sigmaTLPx(i)= sigmaTLPx(i)+
     1                  H5(i,j)*(-H3(i,j)+H2(i,j))+H6(i,j)*(H1(i,j)+H4(i,j))
          sigmaTLPy(i)= sigmaTLPy(i)-
     1                  H5(i,j)*( H3(i,j)+H2(i,j))+H6(i,j)*(H1(i,j)-H4(i,j))
          sigmaTLPz(i)= sigmaTLPz(i)+
     1                  H5(i,j)*( H1(i,j)+H4(i,j))+H6(i,j)*(H3(i,j)-H2(i,j))
          sigmaTTPx(i)= sigmaTTPx(i)+
     1                  H1(i,j)*H2(i,j)+H3(i,j)*H4(i,j)
          sigmaTTPz(i)= sigmaTTPz(i)+
     1                  H1(i,j)**2.- H2(i,j)**2.+ H3(i,j)**2.- H4(i,j)**2.
        EndDo
          sigmaTLP(i) = sigmaTLP(i)+
     1        H5(i,2)*(H1(i,1)-H4(i,1))- H5(i,1)*(H1(i,2)-H4(i,2))+
     1        H6(i,2)*(H2(i,1)+H3(i,1))- H6(i,1)*(H2(i,2)+H3(i,2))
          sigmaTy(i)  = sigmaTy(i)
     1         +H1(i,2)*H2(i,1)-H1(i,1)*H2(i,2)
     1         +H3(i,2)*H4(i,1)-H3(i,1)*H4(i,2)
          sigmaLy(i)  = sigmaLy(i)
     1                  +H5(i,2)*H6(i,1)-H5(i,1)*H6(i,2)
          sigmaTLx(i) = sigmaTLx(i)+
     1         H5(i,1)*(H3(i,2)-H2(i,2))-H5(i,2)*(H3(i,1)-H2(i,1))
     1         -H6(i,1)*(H1(i,2)+H4(i,2))+H6(i,2)*(H1(i,1)+H4(i,1))
          sigmaTLy(i) = sigmaTLy(i)
     1                 -H5(i,1)*(H3(i,2)+H2(i,2))+H5(i,2)*(H3(i,1)+H2(i,1))
     1                 +H6(i,1)*(H1(i,2)-H4(i,2))-H6(i,2)*(H1(i,1)-H4(i,1))
          sigmaTLz(i) = sigmaTLz(i)
     1                 -H5(i,1)*(H1(i,2)+H4(i,2))+H5(i,2)*(H1(i,1)+H4(i,1))
     1                 +H6(i,1)*(H2(i,2)-H3(i,2))-H6(i,2)*(H2(i,1)-H3(i,1))
          sigmaTTx(i) = sigmaTTx(i)
     1                 +H1(i,1)*H3(i,2)-H1(i,2)*H3(i,1)
     1                 -H2(i,1)*H4(i,2)+H2(i,2)*H4(i,1)
          sigmaTTy(i) = sigmaTTy(i)
     1                 -H1(i,1)*H3(i,2)+H1(i,2)*H3(i,1)
     1                 -H2(i,1)*H4(i,2)+H2(i,2)*H4(i,1)
          sigmaTTz(i) = sigmaTTz(i)
     1                 -H1(i,1)*H4(i,2)+H1(i,2)*H4(i,1)
     1                 -H2(i,1)*H3(i,2)+H2(i,2)*H3(i,1)



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

c
      ENDDO
      RETURN
      END
c
c------------------------------------------------------------------------------




