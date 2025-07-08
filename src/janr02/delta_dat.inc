c      include file DELTA_DAT
c
c_begin_doc
c  RCS ID string
c  $Id: delta_dat.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $
c
c  Documentation for include file DELTA_DAT
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
c  Author:         Created:  Mon Dec 16 19:21:20 EST 1996
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
c  variables for include file KINPAR:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIDELTA_DAT
      DATA CIDELTA_DAT/
     1'$Id: delta_dat.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in KINPAR common block
c  --------------------------------------------
c
c_end_var
c
      REAL xM(74),xE(74),xS(74)
      REAL xOM(74),xOE(74),xOS(74),xH(74)
      REAL yM(74),yE(74),yS(74)
      REAL yOM(74),yOE(74),yOS(74),yH(74)
      COMMON/DELTA_DAT/xM,xE,xS,xOM,xOE,xOS,xH,yM,yE,yS,yOM,yOE,yOS,yH
c----6----------------------------------------------------------------72




