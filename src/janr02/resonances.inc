c      include file RESONANCES
c
c_begin_doc
c  RCS ID string
c  $Id: resonances.inc,v 1.1 2003/11/20 05:16:09 lcsmith Exp $
c
c  Documentation for include file RESONANCES
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
c  Author:         Created:  Wed Dec 11 16:45:51 EST 1996
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
c  variables for include file RESONANCES:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIRESONANCES
      DATA CIRESONANCES/
     1'$Id: resonances.inc,v 1.1 2003/11/20 05:16:09 lcsmith Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in RESONANCES common block
c  --------------------------------------------
c
c_end_var
c
       INTEGER MAXres
       PARAMETER (MAXres=40)
       REAL MRES(MAXres), GRES(MAXres), XRES(MAXres), eta(MAXres)
       INTEGER Langul(MAXres), Lprime(MAXres) 
       COMMON/RESONANCES/MRES,GRES,XRES,Langul,Lprime,eta
c----6----------------------------------------------------------------72




