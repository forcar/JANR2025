      SUBROUTINE SAFE_EXIT
c
c_begin_doc
c  RCS ID string
c  $Id: safe_exit.F,v 1.1 2003/11/20 05:16:09 lcsmith Exp $
c
c  Documentation for subroutine SAFE_EXIT
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
c  Author:  H. Aznauryan       Created:  Wed Dec 11 17:43:38 EST 1996
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
      PARAMETER (CFILE = '$RCSfile: safe_exit.F,v $')
      PARAMETER (CREVIS = '$Revision: 1.1 $')
      PARAMETER (CSTATE = '$State: Exp $')
      PARAMETER (CDATE = '$Date: 2003/11/20 05:16:09 $')
      PARAMETER (CAUTHO = '$Author: lcsmith $')
      DATA CRCSID/   
     1'$Id: safe_exit.F,v 1.1 2003/11/20 05:16:09 lcsmith Exp $'   
     2/   
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='SAFE_EXIT')
      PARAMETER (CRAUTH='')
c
c  Local User defined variables:
c  -----------------------------
c_end_var
c
c  executable code for routine safe_exit:
c----6----------------------------------------------------------------72
       write(f_out,*)'Program is interupted by user'
       write(6,*)'Program is interupted by user'
       Close(f_out)
       END
c------------------------------------------------------------------------------




