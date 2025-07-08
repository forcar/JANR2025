c      include file AMPLITUDES
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for include file AMPLITUDES
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
c  variables for include file AMPLITUDES:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIAMPLITUDES
      DATA CIAMPLITUDES/
     1'$Id$'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in AMPLITUDES common block
c  --------------------------------------------
c
c_end_var
c
       INTEGER MAX1res
       PARAMETER (MAX1res=13)
       INTEGER MAX3res
       PARAMETER (MAX3res=13)
       INTEGER MAX13res
       PARAMETER (MAX13res=40)
       REAL AM3,AE3,AS3  ! 
       REAL AM1,AE1,AS1,ratio  ! 
       REAL BreitE3,BreitM3,BreitS3
       REAL BreitE1,BreitM1,BreitS1
       COMMON/AMPLITUDES/BreitE3(2,MAX3res),BreitE1(2,2,MAX1res),
     +      BreitM3(2,MAX3res),BreitS3(2,MAX3res),
     +      BreitM1(2,2,MAX1res),BreitS1(2,2,MAX1res),
     +      AM3(MAX3res),AE3(MAX3res),AS3(MAX3res),
     +      AM1(2,MAX1res),AE1(2,MAX1res),AS1(2,MAX1res),ratio
c----6----------------------------------------------------------------72




