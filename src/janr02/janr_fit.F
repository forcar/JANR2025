      SUBROUTINE JANR_FIT

      IMPLICIT NONE
      SAVE

#include "fit_par.inc"

      INTEGER I,ERRORflag
      double precision stv,stsz,dzero,ar_list(35),llim,ulim
      data dzero/0.D00/
      EXTERNAL JANR_fmin
      real fputil
       
c     Initialize Minuit, Define I/O
      CALL mninit(lmini,lmino,lminf)
          
c     Define the title of the Job
      CALL mnseti('JANR Analysis of Single Pion Production')
          
c     Define parameters, Set the initial values
      fix_param=0.
      do i=1,parms
        stv  = start_value(i)
        stsz = step_size(i)
        llim = low_lim(i)
        ulim = up_lim(i)
        call mnparm(plist(i),pname(i),stv,stsz,llim,ulim,errorflag)
        if(errorflag.ne.0)then
        write(6,*)'Unable to define parameter  ',pname(i),
     +  i,' Error',errorflag
        endif
      enddo
      
      ar_list(1)=1.
      CALL mnexcm(JANR_fmin,'CALL UFC',ar_list,1,errorflag,fputil)
      
      do i=1,parms
        ar_list(1)=i
c       ar_list(2)=100
        ar_list(2)=10
        if(abs(par_stat(i)).eq.1)then
          print *,'Scaning parameter ',i
          call mnexcm(JANR_fmin,'SCAN',ar_list,2,errorflag,fputil)
        endif
        if(par_stat(i).ge.1)then
          fix_param=fix_param+1.
          ar_list(1)=float(Plist(i))
          CALL mnexcm(JANR_fmin,'FIX',ar_list,1,errorflag)
          if(errorflag.ne.0)then
            write(6,*)'Unable to Fix the parameter  ',pname(i),
     +      i,' Error',errorflag
          endif
        endif
      enddo
   
      print *,'Number of parameters = ',parms,'  Number of fixed parameters = ',fix_param

      fputil=0.
      ar_list(1)=1.
      CALL mnexcm(JANR_fmin,'SET PRINT',ar_list,1,errorflag,fputil)
     
      ar_list(1)=errorup
      CALL mnexcm(JANR_fmin,'SET ERRORDEF',ar_list,1,errorflag,fputil)
     
      ar_list(1)=2.
      CALL mnexcm(JANR_fmin,'CALL UFC',ar_list,1,errorflag,fputil)
     
      ar_list(1)=maxcalls
      ar_list(2)=tolorance
      write(6,*)'Starting fit', ar_list(1), ar_list(2)
      CALL mnexcm(JANR_fmin,'MIG',ar_list,2,fputil,errorflag)
      
      write(6,*),'MIGRAD',errorflag
      write(lmino,*)'MIGRAD',errorflag
     
      if(errorflag.ne.0)then
        ar_list(1)=tolorance
        ar_list(2)=1
        ar_list(3)=2
        ar_list(4)=3
        CALL mnexcm(JANR_fmin,'MINOS',ar_list,1,errorflag,fputil)
        write(6,*),'MINOS',errorflag
        write(lmino,*),'MINOS',errorflag
      endif
  
      print *,'Sending STOP to MINUI '
      CALL mnexcm(JANR_fmin,'STOP',0,0,fputil,errorflag)
      
      end
      
      subroutine prterr

      implicit none
#include "fit_par.inc"
      double precision fmin,fedm,errdef,eplus,eminus,eparab,globcc,parv,parerr,bn1,bn2
      integer npari,nparx,istat,i,iv
      CHARACTER*10 pnamei(35)

      call mnstat(fmin,fedm,errdef,npari,nparx,istat)
      write(6,*)
      write(6,*)fmin,fmin/ndf,fedm,errdef,npari,nparx,istat
      write(6,*)
      write(6,*)'Parameter error'

      do i=1,parms
         call mnpout(i,pnamei(i),parv,parerr,bn1,bn2,iv)
         if(par_stat(i).le.0)then
            call mnerrs(-i,eplus,eminus,eparab,globcc)
            write(6,*)i,'  ',pnamei(i),parv,parerr,bn1,bn2
            write(6,*)i,eplus,eminus,eparab,globcc
            write(f_out,*)i,'  ',pnamei(i),parv,parerr
         else
c           write(f_out,*)i,'  ',pnamei(i),parv
         endif
      enddo

      return
      end









