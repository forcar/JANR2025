c      include file MULTIPOLES
c
c_begin_doc
c  RCS ID string
c  $Id: multipoles.inc,v 1.1 2003/11/20 05:16:09 lcsmith Exp $
c
c  Documentation for include file MULTIPOLES
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
c  variables for include file MULTIPOLES:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIMULTIPOLES
      DATA CIMULTIPOLES/
     1'$Id: multipoles.inc,v 1.1 2003/11/20 05:16:09 lcsmith Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in MULTIPOLES common block
c  --------------------------------------------
c
c_end_var
c
       REAL ReacM2,ReacM3,ReacM4,ReacM5,ReacM6,ReacM7  ! 
       REAL ReacE1,ReacE3,ReacE4,ReacE5,ReacE6,ReacE7  ! 
       REAL ReacS1,ReacS2,ReacS3,ReacS4,ReacS5,ReacS6,ReacS7  ! 
       REAL MM3,ME3,MS3  ! 
       REAL MM1,ME1,MS1  ! 
       REAL WSa(51)
       REAL mp33r(51),mp33i(51)
       REAL s11r(51),s11i(51)
       REAL s31r(51),s31i(51)
       REAL p11r(51),p11i(51)
       REAL p13r(51),p13i(51)
       REAL p31r(51),p31i(51)
       REAL p33r(51),p33i(51)
       REAL d13r(51),d13i(51)
       REAL d15r(51),d15i(51)
       REAL d33r(51),d33i(51)
       REAL f15r(51),f15i(51)
       REAL f37r(51),f37i(51)
       REAL simr(13),simi(13)
       REAL pimr(13),pimi(13)
       COMMON/MULTIPOLES/ ReacM2(2,2),ReacM3(2,2),ReacM4(2,2),
     +      ReacM5(2,2),ReacM6(2,2),ReacM7(2,2),
     +      ReacE1(2,2),ReacE3(2,2),ReacE4(2,2),
     +      ReacE5(2,2),ReacE6(2,2),ReacE7(2,2),
     +      ReacS1(2,2),ReacS2(2,2),ReacS3(2,2),ReacS4(2,2),
     +      ReacS5(2,2),ReacS6(2,2),ReacS7(2,2),
     +      MM3(2,8),ME3(2,8),MS3(2,8),
     +      MM1(2,2,9),ME1(2,2,9),MS1(2,2,9),
     +      mp33r,mp33i,
     +      s11r,s11i,s31r,s31i,
     +      p11r,p11i,p13r,p13i,p31r,p31i,p33r,p33i,
     +      d13r,d13i,d15r,d15i,d33r,d33i,
     +      f15r,f15i,f37r,f37i,
     +      simr,simi,pimr,pimi,WSa
c----6----------------------------------------------------------------72




