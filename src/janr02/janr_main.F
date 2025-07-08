      PROGRAM JANR_MAIN

      IMPLICIT NONE
      SAVE

#include "inpar.inc"

      INTEGER ierr

      CALL assign_sig    
      CALL janr_init(ierr)
      call mult_ana
      
      IF( fitopt.NE.'M' .AND. fitopt.NE.'m' ) THEN
      
         call read_multables
         call regge_multables
         call high_multables
         
         If( fitopt.EQ.'F' .OR. fitopt.EQ.'f' .OR. 
     +        fitopt.EQ.'MF' .OR. fitopt.EQ.'mf' ) Then
            CALL janr_FIT
         ElseIf( fitopt.EQ.'G' .OR. fitopt.EQ.'g'.or. 
     +           fitopt.EQ.'GP' .OR. fitopt.EQ.'gp'.or.
     +           fitopt.EQ.'GR' .OR. fitopt.EQ.'gr') Then
            CALL janr_EVGen
         EndIf
         
c     ELSEIF( fitopt.EQ.'M' .OR. fitopt.EQ.'m' ) THEN
c      call mult_ana
     
      ENDIF
      
      END



