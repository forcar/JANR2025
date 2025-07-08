      SUBROUTINE JANR_FIT
c
c_begin_doc
c  RCS ID string
c  $Id: JANR_main.F 1.1 1996/12/11 20:41:20 stepanya Exp $
c
c  Documentation for subroutine JANR_MAIN
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
c  Author: S.Stepanyan   Created:  Wed Jun 23 15:40:59 EST 1996
c  -------
c
c  Major revisions:
c  ----------------
c     
c
c_end_doc`
c
      IMPLICIT NONE
      SAVE
c
c_begin_inc
c  include files :
c  ---------------------
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
      PARAMETER (CFILE = '$RCSfile: JANR_main.F $')
      PARAMETER (CREVIS = '$Revision: 1.1 $')
      PARAMETER (CSTATE = '$State: Exp $')
      PARAMETER (CDATE = '$Date: 1996/12/11 20:41:20 $')
      PARAMETER (CAUTHO = '$Author: stepanya $')
      DATA CRCSID/   
     1'$Id: JANR_main.F 1.1 1996/12/11 20:41:20 stepanya Exp $'   
     2/   
c  Module information:
      CHARACTER*(*)  CRNAME, CRAUTH
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='JANR_MAIN')
      PARAMETER (CRAUTH='stepanya')
c
c  Local User defined variables:
c  -----------------------------
       INTEGER I,ERRORflag
       double precision stv,stsz,dzero,ar_list(35),llim,ulim
       data dzero/0.D00/
       EXTERNAL JANR_fmin
       real fputil
c_end_var
c
c  executable code for routine JANR_MAIN:
c----6----------------------------------------------------------------72
c Initialize Minuit, Define I/O
          CALL mninit(lmini,lmino,lminf)
c     Define the title of the Job
          CALL mnseti('Dispersion Analyses of Single pion Production')
c     Define parameters, Set the initial values
          fix_param=0.
          DO i=1,parms
             stv=start_value(i)
             stsz=step_size(i)
             llim=low_lim(i)
             ulim=up_lim(i)
             call mnparm(Plist(i),Pname(i),stv,stsz,llim,ulim,errorflag)
             if(errorflag.ne.0)then
                write(6,*)'Unable to define parameter  ',pname(i),
     +               i,' Error',errorflag
             endif
          ENDDO
          ar_list(1)=1.
          CALL mnexcm(JANR_fmin,'call ufc',ar_list,1,errorflag,fputil)
          DO i=1,parms
             ar_list(1)=i
c             ar_list(2)=100
             ar_list(2)=10
             if(abs(par_stat(i)).eq.1)then
                print *,'Scaning parameter ',i
                call mnexcm(JANR_fmin,'scan',ar_list,2,errorflag,fputil)
             endif
             if(par_stat(i).ge.1)then
                fix_param=fix_param+1.
                ar_list(1)=float(Plist(i))
                CALL mnexcm(JANR_fmin,'FIX',ar_list,1,errorflag)
                if(errorflag.ne.0)then
                   write(6,*)'Unable to Fix the parameter  ',pname(i),
     +                  i,' Error',errorflag
                endif
             endif
          ENDDO
c   
          print *,'Number of parameters = ',parms,'  Number of fixed parameters = ',fix_param
          fputil=0.
          ar_list(1)=1.
          CALL mnexcm(JANR_fmin,'set print',ar_list,1,errorflag,fputil)
c     
          ar_list(1)=errorup
          CALL mnexcm(JANR_fmin,'set errordef',ar_list,1,errorflag,fputil)
c     
          ar_list(1)=2.
          CALL mnexcm(JANR_fmin,'call ufc',ar_list,1,errorflag,fputil)
c     
          ar_list(1)=maxcalls
          ar_list(2)=tolorance
          write(6,*)'Starting fit', ar_list(1), ar_list(2)
          CALL mnexcm(JANR_fmin,'MIG',ar_list,2,fputil,errorflag)
          write(6,*),'MIGRAD',errorflag
          write(lmino,*)'MIGRAD',errorflag
c     
	if(errorflag.ne.0)then
          ar_list(1)=tolorance
          ar_list(2)=1
          ar_list(3)=2
          ar_list(4)=3
          CALL mnexcm(JANR_fmin,'MINOS',ar_list,1,errorflag,fputil)
          write(6,*),'MINOS',errorflag
          write(lmino,*),'MINOS',errorflag
	endif
C    
        print *,'Sending STOP to MINUI '
        CALL mnexcm(JANR_fmin,'STOP',0,0,fputil,errorflag)
       END
c
c------------------------------------------------------------------------------
      Subroutine prterr
c
      Implicit none
#include "fit_par.inc"
      Double precision fmin,fedm,errdef,eplus,eminus,eparab,globcc,parv,parerr,bn1,bn2
      Integer npari,nparx,istat,i,iv
      CHARACTER*10 pnamei(35)
c
      call mnstat(fmin,fedm,errdef,npari,nparx,istat)
      write(6,*)
      write(6,*)fmin,fmin/ndf,fedm,errdef,npari,nparx,istat
      write(6,*)
      write(6,*)'Parameter error'
c
      do i=1,parms
         call mnpout(i,pnamei(i),parv,parerr,bn1,bn2,iv)
         if(par_stat(i).le.0)then
            call mnerrs(-i,eplus,eminus,eparab,globcc)
            write(6,*)i,'  ',pnamei(i),parv,parerr,bn1,bn2
            write(6,*)i,eplus,eminus,eparab,globcc
            write(f_out,*)i,'  ',pnamei(i),parv,parerr
         else
c            write(f_out,*)i,'  ',pnamei(i),parv
         endif
      enddo
c
      return
      end









