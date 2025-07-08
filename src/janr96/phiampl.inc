c      include file PHIAMPL
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for include file PHIAMPL
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
c  Author:         Created:  Wed Dec 11 16:56:50 EST 1996
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
c  variables for include file PHIAMPL:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIPHIAMPL
      DATA CIPHIAMPL/
     1'$Id$'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in PHIAMPL common block
c  --------------------------------------------
c
c_end_var
c
       REAL phi_amp  ! 
       COMMON/PHIAMPL/phi_amp(6,2,3)
c----6----------------------------------------------------------------72




