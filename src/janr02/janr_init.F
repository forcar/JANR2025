      SUBROUTINE JANR_INIT (ierr)

      IMPLICIT NONE
      SAVE

#include "janr_data.inc"
#include "delta_dat.inc"
#include "inpar.inc"
#include "fit_par.inc"
#include "janrtab.inc"
#include "multipoles.inc"

      INTEGER ierr
      REAL at2,at3,at4,at5,at8,at9     

      CHARACTER*10 dummych
      INTEGER i,LINPUT
      data LINPUT/21/
      integer ipow, itest
      
      CALL lang_def
      
      write(6,*)'Reading the input file'
      
      OPEN (UNIT=LINPUT,FILE='janr_input',STATUS='OLD')
      READ (LINPUT,51)dummych,fitopt 
      Write(6,*)'Action is    ',fitopt
 51   FORMAT(A10,A10)
      READ (LINPUT,*)
      READ (LINPUT,*) Pi_type2fit, obsv2fit
     
      IF((Pi_type2fit.LE.0).OR.(Pi_type2fit.GE.4)) THEN
         WRITE(*,*)'You choice of type of Pion is not correct'
         WRITE(*,*)'Possible choices are [1/2/3] for [Pi0/Pi+/Pi-]'
         RETURN
      ENDIF

      do ipow = 1, 4 
         itest = amod( float(obsv2fit), 10.**(ipow) ) / 10**(ipow-1) 
         if ( itest .GT. 9 ) then
            write(*,*)'Your choice of observables is incorrect'
            write(*,*)'Possible choices are 1, 10, 11, 100, 101 etc '
            return
         endif
      enddo

      READ (LINPUT,*) Ebeam
      READ (LINPUT,*) WMIN,WMAX,WSTEP
      READ (LINPUT,*) Q2MIN,Q2MAX,Q2STEP
      READ (LINPUT,*) COSMIN,COSMAX,COSSTEP
      READ (LINPUT,*) PHIMIN,PHIMAX,PHISTEP
      READ (LINPUT,*) IntAccur
      READ (LINPUT,*) tabfile
      READ (LINPUT,*) inffile
 52   FORMAT(A10,A60)
      READ (LINPUT,*)
      READ (LINPUT,*) MAXEVENTS
      READ (LINPUT,*) MAXITER
      READ (LINPUT,*) Randomseed
c     
      READ (LINPUT,*)
      READ (LINPUT,54)maxcalls
      READ (LINPUT,54)tolorance 
      READ (LINPUT,54)errorup
 54   format(f8.2)
 
      READ (LINPUT,*) parms
      
      If( parms.Le.MAXpar .And. parms.GT.0 ) Then
         READ (LINPUT,*)
         Do i=1,parms
            read(LINPUT,53)plist(i),pname(i), start_value(i)
     +           ,step_size(i),par_stat(i),low_lim(i),up_lim(i)
 53         FORMAT(i2,1x,A10,2f12.2,i10,2f10.4)
         EndDo
      Else
         Write(6,*)'Wrong Number of parameters    ',parms
         Return
      EndIf
      
      READ (LINPUT,*)
      READ (LINPUT,52)datform,Datafile
      READ (LINPUT,*)
      READ (LINPUT,*)lmini,lmino,lminf
      READ (LINPUT,*)chimax
      READ (LINPUT,*)Nbad
      
      If( Nbad.gt.0 )Then
         Do i=1,Nbad
            READ (LINPUT,*)ibad(1,i),ibad(2,i)
         EndDo
      EndIf
      
      call janr_input_said()

      RETURN
      END





