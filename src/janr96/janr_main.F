      PROGRAM JANR_MAIN
c     
c     _begin_doc
c     RCS ID string
c     $Id: JANR_main.F 1.1 1996/12/11 20:41:20 stepanya Exp $
c     
c     Documentation for subroutine JANR_MAIN
c     
c     Purpose:
c     --------
c     
c     Input Parameters:  (Name - Type - Meaning)
c     ----------------
c     
c     Output Parameters:  (Name - Type - Meaning)
c     -----------------
c     
c     Other routines:
c     ---------------
c     
c     Notes:
c     ------
c     
c     Author: S.Stepanyan  Created:  Wed Jun 23 15:40:59 EST 1996
c     -------
c     
c     Major revisions:
c     ----------------
c     
c     
c     _end_doc`
c     
      IMPLICIT NONE
      SAVE
c     
c     _begin_inc
c     include files :
c     ---------------------
#include "inpar.inc"
c     _end_inc
c     
c     _begin_var
c     input/output variables:
c     -----------------------
c     
c     Local pre-defined variables:
c     ---------------------------
c     RCS information: 
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$RCSfile: JANR_main.F $')
      PARAMETER (CREVIS = '$Revision: 1.1 $')
      PARAMETER (CSTATE = '$State: Exp $')
      PARAMETER (CDATE = '$Date: 1996/12/11 20:41:20 $')
      PARAMETER (CAUTHO = '$Author: stepanya $')
      DATA CRCSID/   
     1     '$Id: JANR_main.F 1.1 1996/12/11 20:41:20 stepanya Exp $'   
     2     /   
c     Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='JANR_MAIN')
      PARAMETER (CRAUTH='S.Stepanyan')
c     
c     Local User defined variables:
c     -----------------------------
      INTEGER ierr
c     _end_var
c     
c     executable code for routine JANR_MAIN:
c---- 6----------------------------------------------------------------72
c     Initialisation of interrupt procedure
      CALL assign_sig
c     End of initialisation of interrupt procedure
      CALL janr_init(ierr)
      call mult_ana
      IF( fitopt.NE.'M' .AND. fitopt.NE.'m' ) THEN
         call read_multables
         call regge_multables
         If( fitopt.EQ.'F' .OR. fitopt.EQ.'f' .OR. 
     +        fitopt.EQ.'MF' .OR. fitopt.EQ.'mf' ) Then
            CALL janr_FIT
         ElseIf( fitopt.EQ.'G' .OR. fitopt.EQ.'g'.or. 
     +           fitopt.EQ.'GP' .OR. fitopt.EQ.'gp'.or.
     +           fitopt.EQ.'GR' .OR. fitopt.EQ.'gr') Then
            CALL janr_EVGen
         EndIf
c     ELSEIF( fitopt.EQ.'M' .OR. fitopt.EQ.'m' ) THEN
c     call mult_ana
c     
      ENDIF
      END
c
c------------------------------------------------------------------------------




