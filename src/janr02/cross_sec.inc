c      include file CROSS_SEC
c
c_begin_doc
c  RCS ID string
c  $Id: cross_sec.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $
c
c  Documentation for include file CROSS_SEC
c
c  Purpose:
c  --------
c
c  Used in routines:
c  ----------------
c
c  Notes:
c  ------
c
c  Author:         Created:  Mon Dec 16 23:08:39 EST 1996
c  -------
c
c  Major revisions:
c  ----------------
c
c_end_doc
c
c
c  Remember NO NESTED INCLUDES
c
c  COMMON BLOCK NAME should be the same as include file name!
c
c_begin_var
c  variables for include file CROSS_SEC:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CICROSS_SEC
      DATA CICROSS_SEC/
     1'$Id: cross_sec.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in CROSS_SEC common block
c  --------------------------------------------
c
c_end_var
c
       REAL sigmaT(3),sigmaL(3),sigmaTT(3),sigmaTL(3),sigma(3)
       REAL sigmaP(3),sigmaI(3)
c      Cross section of transverse photons. 
      COMMON/CROSS_SEC/sigmaT,sigmaL,sigmaTT,sigmaTL,sigma,sigmaP,
     +                 sigmaI
c----6----------------------------------------------------------------72




