      SUBROUTINE LANG_DEF
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for subroutine LANG_DEF
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
#include "resonances.inc"
#include "fit_par.inc"
c_end_inc
c
c_begin_var
c  input/output variables:
c  -----------------------
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
      PARAMETER (CRNAME='LANG_DEF')
      PARAMETER (CRAUTH='')
c
c  Local User defined variables:
c  -----------------------------
c      INTEGER NWRITE
      INTEGER i,k
c_end_var
c
c  executable code for routine Lang_def:
c----6----------------------------------------------------------------72
c
c- Langul(i) is the degree of photon momenta in Breit-Wigner formula for 
c- M(l+),M(l-): Langul(i)=l.
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
      eta(2)=0.16
      eta(3)=0.42
      eta(4)=0.1
      eta(5)=0.29
      eta(6)=0.1
      eta(7)=0.12
      eta(8)=0.01
      eta(9)=0.13
      eta(10)=0.16
      eta(11)=0.18
      eta(12)=0.01
      eta(13)=0.01


      eta(21)=0.6
      eta(22)=0.4
      eta(23)=0.54
      eta(23)=1.0
      eta(24)=0.1
      eta(25)=0.6
      eta(26)=0.62
      eta(27)=0.15
      eta(28)=0.06
      eta(29)=0.38
      eta(30)=0.18
      eta(31)=0.2
      eta(32)=0.1
      eta(33)=0.01

      RETURN
      END
c
c------------------------------------------------------------------------------
