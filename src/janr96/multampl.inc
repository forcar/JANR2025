c      include file MULTAMPL
c
c_begin_doc
c  RCS ID string
c  $Id$
c
c  Documentation for include file MULTAMPL
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
c----6----------------------------------------------------------------72
c  
c  Variables contained in MULTAMPL common block
c  --------------------------------------------
c
c_end_var
c
c
       INTEGER Maxmpoints
       PARAMETER (Maxmpoints=101)
       REAL Wtab(Maxmpoints)
       REAL AEm(Maxmpoints,4,3),AEp(Maxmpoints,4,3)   
       REAL ASm(Maxmpoints,4,3),ASp(Maxmpoints,4,3)   
       REAL AMm(Maxmpoints,4,3),AMp(Maxmpoints,4,3)  
       REAL RAEm(Maxmpoints,4,3),RAEp(Maxmpoints,4,3)   
       REAL RASm(Maxmpoints,4,3),RASp(Maxmpoints,4,3)   
       REAL RAMm(Maxmpoints,4,3),RAMp(Maxmpoints,4,3)  
       COMMON/MULTampl/Wtab,AEm,AEp,ASm,ASp,AMm,AMp,
     +     RAEm,RAEp,RASm,RASp,RAMm,RAMp
       save
c----6----------------------------------------------------------------72




