c      include file BACK_AMPL
c
c_begin_doc
c  RCS ID string
c  $Id: back_ampl.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $
c
c  Documentation for include file BACK_AMPL
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
c  Author:         Created:  Thu Dec 12 19:51:49 EST 1996
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
c  variables for include file BACK_AMPL:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIBACK_AMPL
      DATA CIBACK_AMPL/
     1'$Id: back_ampl.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in ADD_AMPL common block
c  --------------------------------------------
c
c_end_var
c
       REAL BM3,BE3,BS3
       REAL BM1,BE1,BS1
       REAL BbM3,BbE3,BbS3
       REAL BbM1,BbE1,BbS1
       REAL RM3,RE3,RS3
       REAL RM1,RE1,RS1
       REAL E0,M1,S0,S1


	COMMON/BACK_AMPL/BM3(8),BE3(8),BS3(8),
     +      BM1(2,9),BE1(2,9),BS1(2,9),
     +      BbM3(8),BbE3(8),BbS3(8),
     +      BbM1(2,9),BbE1(2,9),BbS1(2,9),
     +      RM3(8),RE3(8),RS3(8),
     +      RM1(2,9),RE1(2,9),RS1(2,9),
     +      E0(2),M1(2),S0(2),S1(2)
c----6----------------------------------------------------------------72




