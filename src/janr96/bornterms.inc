c      include file BORN
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for include file BORN
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
c  Author:         Created:  Mon Dec 16 11:13:05 EST 1996
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
c  variables for include file BORN:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIBORN
      DATA CIBORN/
     1'$Id$'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in BORN common block
c  --------------------------------------------
c
c_end_var
c
      REAL br(2,6,3) 
       COMMON/BORNTERMS/br
c----6----------------------------------------------------------------72




