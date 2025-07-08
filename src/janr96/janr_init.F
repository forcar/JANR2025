      SUBROUTINE JANR_INIT (ierr)
c     
c     _begin_doc
c     RCS ID string
c     $Id$
c     
c     Documentation for subroutine JANR_INIT
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
c     Author:         Created:  Wed Dec 11 17:16:27 EST 1996
c     -------
c     
c     Major revisions:
c     ----------------
c     
c     
c     _end_doc
c     
      IMPLICIT NONE
      SAVE
c     
c     _begin_inc
c     include files :
c     ---------------------
c     BOS common block  uncomment the next line for BOS include file
#include "janr_data.inc"
#include "delta_dat.inc"
#include "inpar.inc"
#include "fit_par.inc"
#include "janrtab.inc"
#include "multipoles.inc"
c     _end_inc
c     
c     _begin_var
c     input/output variables:
c     -----------------------
      INTEGER ierr
       REAL at2,at3,at4,at5,at8,at9     
c     Local pre-defined variables:
c     ---------------------------
c     RCS information: 
      CHARACTER*132  CFILE, CREVIS, CSTATE, CDATE, CAUTHO, CRCSID
      PARAMETER (CFILE = '$RCSfile$')
      PARAMETER (CREVIS = '$Revision$')
      PARAMETER (CSTATE = '$State$')
      PARAMETER (CDATE = '$Date$')
      PARAMETER (CAUTHO = '$Author$')
      DATA CRCSID/   
     1     '$Id$'   
     2     /   
c     Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      PARAMETER (CRNAME='JANR_INIT')
      PARAMETER (CRAUTH='stepanyan')
c     
c     Local User defined variables:
c     -----------------------------
      CHARACTER*10 dummych
      INTEGER i,LINPUT
      data LINPUT/21/
      integer ipow, itest
c     _end_var
c     
c     executable code for routine JANR_INIT:
c---- 6----------------------------------------------------------------72
c     Begin of constants input
      CALL lang_def
c     End of constants input

      write(6,*)'Reading the input file'
      OPEN (UNIT=LINPUT,FILE='janr_input',STATUS='OLD')
      READ (LINPUT,51)dummych,fitopt 
      Write(6,*)'Action is    ',fitopt
 51   FORMAT(A10,A10)
      READ (LINPUT,*)
      READ (LINPUT,*) Pi_type2fit, obsv2fit
C     
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
c     Read SAID tables
      open(UNIT=67,file='/home/ampana/janr/multipols/mp33.dat',status='old')
      do i=1,51
      read(67,*,end=110)wsa(i),mp33r(i),mp33i(i),at2,at3,at4
	 wsa(i)=wsa(i)/1000.
	 mp33r(i)=mp33r(i)/52.437
	 mp33i(i)=mp33i(i)/52.437
      enddo
 110  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/s11.dat',status='old')
      do i=1,51
      read(67,*,end=111)wsa(i),at2,at3,at4,at5,s11r(i),s11i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 111  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/s31.dat',status='old')
      do i=1,51
      read(67,*,end=112)wsa(i),at2,at3,at4,at5,s31r(i),s31i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 112  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/p11.dat',status='old')
      do i=1,51
      read(67,*,end=113)wsa(i),at2,at3,at4,at5,p11r(i),p11i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 113  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/p13.dat',status='old')
      do i=1,51
      read(67,*,end=114)wsa(i),at2,at3,at4,at5,p13r(i),p13i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 114  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/p31.dat',status='old')
      do i=1,51
      read(67,*,end=115)wsa(i),at2,at3,at4,at5,p31r(i),p31i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 115  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/p33.dat',status='old')
      do i=1,51
      read(67,*,end=116)wsa(i),at2,at3,at4,at5,p33r(i),p33i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 116  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/d13.dat',status='old')
      do i=1,51
      read(67,*,end=117)wsa(i),at2,at3,at4,at5,d13r(i),d13i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 117  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/d15.dat',status='old')
      do i=1,51
      read(67,*,end=118)wsa(i),at2,at3,at4,at5,d15r(i),d15i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 118  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/d33.dat',status='old')
      do i=1,51
      read(67,*,end=119)wsa(i),at2,at3,at4,at5,d33r(i),d33i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 119  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/f15.dat',status='old')
      do i=1,51
      read(67,*,end=120)wsa(i),at2,at3,at4,at5,f15r(i),f15i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 120  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/f37.dat',status='old')
      do i=1,51
      read(67,*,end=121)wsa(i),at2,at3,at4,at5,f37r(i),f37i(i),at8,at9
	 wsa(i)=wsa(i)/1000.
      enddo
 121  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/pim.dat',status='old')
      do i=1,13
         read(67,*,end=122)wsa(i),pimr(i),pimi(i)
         wsa(i)=wsa(i)/1000.
      enddo
 122  continue
      close(67)
      open(UNIT=67,file='/home/ampana/janr/multipols/sim.dat',status='old')
      do i=1,13
         read(67,*,end=123)wsa(i),simr(i),simi(i)
         wsa(i)=wsa(i)/1000.
      enddo
 123  continue
      close(67)

      RETURN
      END
c
c------------------------------------------------------------------------------




