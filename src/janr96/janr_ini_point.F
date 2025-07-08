      SUBROUTINE JANR_INI_point(Q2)
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for subroutine JANR_INI_POINT
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
c  Author:         Created:  Wed Dec 11 17:16:27 EST 1996
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
#include "resonances.inc"
#include "amplitudes.inc"
#include "fit_par.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
c
      REAL Q2
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
      PARAMETER (CRNAME='JANR_INI_POINT')
      PARAMETER (CRAUTH='')
c
c  Local User defined variables:
c  -----------------------------
c      INTEGER NWRITE
      INTEGER i,j,k1
	
c_end_var
c
c  executable code for routine JANR_INI_POINT:
c----6----------------------------------------------------------------72
c
      pi=2.*acos(0.E0)
      a2=sqrt(2.E0)
      mn=0.93827
      if(pi_type.ge.2)then
         mp=0.13957
      else
         mp=0.13498
      endif
      meta=0.547
      m2pion=2.*mp
      mn22=mn**2
      mp22=mp**2
      mn2=2.*mn
C-Resonances with (****),(***),(**) are included
c-(1)	P33(1232)
c-(2)	D33(1700)
c-(3)	F37(1950)
c-(4)	P33(1600)
c-(5)	S31(1620)
c-(6)	F35(1905)
c-(7)	D35(1930)
c-(8)	P31(1925) missing resonance
c-(9)	S31(1900)
c-(10)	P33(1920)+ missing resonance (1975)
c-(11)	P31(1910)
c-(12)	D33(1940)
c-(13)	F35(2000)
c-(21)	P11(1440)
c-(22) 	S11(1535)
c-(23) 	D13(1520)
c-(24) 	D13(1700)
c-(25) 	S11(1650)
c-(26) 	F15(1680)
c-(27) 	P13(1720)
c-(28) 	F17(1990)
c-(29) 	D15(1675)
c-(30)	P11(1710)
c-(31)	P13(1900)+3 missing resonances (1870),(1980),(1955)
c-(32)	D13(2080)
c-(33)	F15(2000)+ missing resonance (1955)
      MRES(1)=1.228
      MRES(2)=1.65
      MRES(3)=1.93
      MRES(4)=1.6
      MRES(5)=1.62
      MRES(6)=1.905
      MRES(7)=1.93
      MRES(8)=1.925
      MRES(9)=1.9
      MRES(10)=1.92
      MRES(11)=1.91
      MRES(12)=1.94
      MRES(13)=2.
C      
      MRES(21)=1.4255
      MRES(22)=1.52544
      MRES(23)=1.51455
      MRES(24)=1.68
      MRES(25)=1.65
      MRES(26)=1.682
      MRES(27)=1.8
      MRES(28)=1.99
      MRES(29)=1.675
      MRES(30)=1.71
      MRES(31)=1.9
      MRES(32)=1.97
      MRES(33)=2.
C	
      GRES(1)=0.11
      GRES(2)=0.25
      GRES(3)=0.3
      GRES(4)=0.3
      GRES(5)=0.14
      GRES(6)=0.35
      GRES(7)=0.35
      GRES(8)=0.3
      GRES(9)=0.2
      GRES(10)=0.2
      GRES(11)=0.25
      GRES(12)=0.2
      GRES(13)=0.3
C	
      GRES(21)=0.3396
      GRES(22)=0.11808
      GRES(23)=0.1086
      GRES(24)=0.1
      GRES(25)=0.11
      GRES(26)=0.13
      GRES(27)=0.32
      GRES(28)=0.35
      GRES(29)=0.16
      GRES(30)=0.11
      GRES(31)=0.5
      GRES(32)=0.3
      GRES(33)=0.2
C
C	
      do  i=2,40
         XRES(i)=0.5
      end do
	XRES(1)=0.2
	XRES(21)=0.375
	XRES(22)=0.5
	XRES(23)=1.2651
	XRES(27)=0.4
	XRES(2)=0.2
c	For I=3/2 AM3=M(3/2)
c	For I=1/2 AAM1(1,i)=M(1/2)/3.,	AAM1(2,i)=M(0)
C	cMi,cEi,cSi are fitting parameters for multipoles with I=3/2
C	cMPi,cEPi,cSPi are fitting parameters for multipoles pA(1/2)
C	P33(1232)
      AM3(1)=4.7276
      AE3(1)=-0.02505
      AS3(1)=-0.15265
C      D33(1700)
      AM3(2)=0.019
      AE3(2)=-0.27
      AS3(2)=0.065
C	F37(1950)
      AM3(3)=0.17
      AE3(3)=-0.006
      AS3(3)=-0.065
C       P33(1600)  
      AM3(4)=0.    
      AE3(4)=0.
      AS3(4)=0.
C       S31(1620)
      AE3(5)=0.
      AS3(5)=0.
C       F35(1905)
      AM3(6)=0.04
      AE3(6)=0.04
      AS3(6)=0.04   
C       D35(1930)
      AM3(7)=0.
      AE3(7)=0.04
      AS3(7)=0.
C       P31(1925)
      AM3(8)=0.
      AE3(8)=-0.04   
      AS3(8)=0.   
C	
      DO k1=9,13
         as3(k1)=0.
         ae3(k1)=0.
         am3(k1)=0.
      ENDDO
C	P11(1440)
      AM1(1,1)=0.0507
      AS1(1,1)=0.14324
C	S11(1535)
      AE1(1,2)=0.65466
      AS1(1,2)=-0.257668
C	D13(1520)
      AM1(1,3)=0.2458
      AE1(1,3)=0.17973
      AS1(1,3)=-0.162

C	D13(1700)
      AM1(1,4)=-0.06
      AE1(1,4)=-0.06
      AS1(1,4)=-0.036

C      S11(1650)
      AE1(1,5)=0.148
      AS1(1,5)=0.0

C	F15(1680)
      AM1(1,6)=0.112
      AE1(1,6)=0.143
      AS1(1,6)=-0.033          
C      P13(1720)
      AM1(1,7)=0.06
      AE1(1,7)=0.007
      AS1(1,7)=-0.06
C	
      DO i=1,2
         DO j=8,13
            am1(i,j)=0.
            ae1(i,j)=0.
         ENDDO
      ENDDO
      RETURN
      END
c
c------------------------------------------------------------------------------




