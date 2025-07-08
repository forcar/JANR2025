c      include file INT_FUNC
c
c_begin_doc
c  RCS ID string
c  $Id: int_func.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $
c
c  Documentation for include file INT_FUNC
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
c  Author:         Created:  Mon Dec 16 11:41:50 EST 1996
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
c  variables for include file INT_FUNC:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIINT_FUNC
      DATA CIINT_FUNC/
     1'$Id: int_func.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in INT_FUNC common block
c  --------------------------------------------
c
c_end_var
c
       REAL H(6,3)
       INTEGER jm,km
       COMMON/INT_FUNC/H,jm,km
c----6----------------------------------------------------------------72




