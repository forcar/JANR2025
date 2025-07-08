      Subroutine mult_ana
c     
c     Documentation for subroutine Mult_ana
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
c     -------------
c     
c     Notes:
c     ------
c     
c     Author:  S.Stepanyan     Created:  Wed Jun 23 17:16:39 EST 1996
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
#include "inpar.inc"
#include "janrtab.inc"
#include "fit_par.inc"
#include "remult.inc"
#include "multampl.inc"
c     
      Integer Ierr
c     
c     Local User defined variables:
c     -----------------------------
c     
      INTEGER iw,iq2,icos,iphi,i,j
      INTEGER Nw,Nq2,Ncos,Nphi
      REAL Q2, W, CosTh, Phi, Epsilon, nu, ep, mp, El
      REAL Fint, gauss, W_min, W_step
      external gauss,fint
      integer lm,is,ii,jj,ifile,lun
      character*10 fnamr(24)
      data fnamr/'rm0p.dat','rm0m.dat','re0p.dat','re0m.dat','rs0p.dat','rs0m.dat',
     +     'rm1p.dat','rm1m.dat','re1p.dat','re1m.dat','rs1p.dat','rs1m.dat',
     +     'rm2p.dat','rm2m.dat','re2p.dat','re2m.dat','rs2p.dat','rs2m.dat',
     +     'rm3p.dat','rm3m.dat','re3p.dat','re3m.dat','rs3p.dat','rs3m.dat'/
      character*10 fnamrr(24)
      data fnamrr/'rrm0p.dat','rrm0m.dat','rre0p.dat','rre0m.dat',
     +     'rrs0p.dat','rrs0m.dat',
     +     'rrm1p.dat','rrm1m.dat','rre1p.dat','rre1m.dat','rrs1p.dat','rrs1m.dat',
     +     'rrm2p.dat','rrm2m.dat','rre2p.dat','rre2m.dat','rrs2p.dat','rrs2m.dat',
     +     'rrm3p.dat','rrm3m.dat','rre3p.dat','rre3m.dat','rrs3p.dat','rrs3m.dat'/
      common/ifint/lm,is,ii,jj,w,q2
      DATA mp/0.938/
c     _end_var
c     
c     executable code for routine Mult_ana:
c---- 6----------------------------------------------------------------72
      ierr=0
      cM1=start_value(1)
      cE1=start_value(2)
      cS1=start_value(3)
      cM2=start_value(4)
      cE2=start_value(5)
      cS2=start_value(6)
      cM3=start_value(7)
      cE3=start_value(8)
      cS3=start_value(9)
      cMP1=start_value(10)
      cSP1=start_value(11)
      cEP2=start_value(12)
      cSP2=start_value(13)
      cMP3=start_value(14)
      cEP3=start_value(15)
      cSP3=start_value(16)
      cMP4=start_value(17)
      cEP4=start_value(18)
      cSP4=start_value(19)
      cEP5=start_value(20)
      cSP5=start_value(21)
      cMP6=start_value(22)
      cEP6=start_value(23)
      cSP6=start_value(24)
      cMP7=start_value(25)
      cEP7=start_value(26)
      cSP7=start_value(27)
      cMP12=start_value(28)
      cEP12=start_value(29)
      cSP12=start_value(30)
      write(6,111)(start_value(i),i=1,parms)
 111  format(25x,'Value of parameters => 1-30 ',/,
     +     '   cM1 =',f9.5,'   cE1 =',f9.5,'   cS1 =',f9.5,/,
     +     '   cM2 =',f9.5,'   cE2 =',f9.5,'   cS2 =',f9.5,/,
     +     '   cM3 =',f9.5,'   cE3 =',f9.5,'   cS3 =',f9.5,/,
     +     '   cMP1 =',f9.5,'   cSP1 =',f9.5,'   cEP2 =',f9.5,'   cSP2 =',f9.5,/,
     +     '   cMP3 =',f9.5,'   cEP3 =',f9.5,'   cSP3 =',f9.5,/,
     +     '   cMP4 =',f9.5,'   cEP4 =',f9.5,'   cSP4 =',f9.5,/,
     +     '   cEP5 =',f9.5,'   cSP5 =',f9.5,/,
     +     '   cMP6 =',f9.5,'   cEP6 =',f9.5,'   cSP6 =',f9.5,/,
     +     '   cMP7 =',f9.5,'   cEP7 =',f9.5,'   cSP7 =',f9.5,/,
     +     '   cMP12 =',f9.5,'   cEP12 =',f9.5,'   cSP12 =',f9.5,/)
C     
      IF( fitopt.NE.'M' .AND. fitopt.NE.'m' ) THEN
         Nw = 101
         W_min=1.1
         W_step=0.01
      ELSE
         If ( Wstep.GT.0 .AND. Wmin+Wstep.LT.Wmax )Then
            Nw = INT( (Wmax - Wmin) / Wstep ) + 1
            W_min=Wmin
            W_step=Wstep
         Else
            if( Wmin.Eq.Wmax )then
               Nw=1
               W_min=Wmin
               W_step=Wstep
            endif
         EndIf
      ENDIF
      write(6,*)'Number of W points = ',nw
c
         do jj=1,1
      do lm=0,3
         do is=1,6
            lun=(lm+2)*10+is
            ifile=lm*6+is
            open(unit=lun,file=fnamr(ifile),status='UNKNOWN')
            write(lun,*)fnamr(ifile)
            write(lun,*)
         enddo
      enddo   
c     
      Q2 = (Q2max+Q2min)/2. 
c
      Do iw=1,Nw
c     Loop over W
         W = W_min + float(iw-1)*W_step
         el=(w**2-mp**2)/2./mp
         do lm=0,3
            do is=1,3
               do ii=1,6
                  Rint(lm,is,ii,jj)=gauss(Fint,-1.,1.,IntAccur)
               enddo
            enddo
            do is=1,3
               call re_mult()
            enddo
            lun=(lm+2)*10
            write(lun+1,201)W,El,(Mlp(lm,is,jj),is=1,3)
            write(lun+2,201)W,El,(Mlm(lm,is,jj),is=1,3)
            write(lun+3,201)W,El,(Elp(lm,is,jj),is=1,3)
            write(lun+4,201)W,El,(Elm(lm,is,jj),is=1,3)
            write(lun+5,201)W,El,(Slp(lm,is,jj),is=1,3)
            write(lun+6,201)W,El,(Slm(lm,is,jj),is=1,3)
         enddo
      EndDo
      do lm=0,3
         do is=1,6
            lun=(lm+2)*10+is
            close(lun)
         enddo
      enddo   
      enddo   
         do jj=2,2
      do lm=0,3
         do is=1,6
            lun=(lm+2)*10+is
            ifile=lm*6+is
            open(unit=lun,file=fnamrr(ifile),status='UNKNOWN')
            write(lun,*)fnamrr(ifile)
            write(lun,*)
         enddo
      enddo   
c     
      Q2 = (Q2max+Q2min)/2. 
c
      Do iw=1,Nw
c     Loop over W
         W = W_min + float(iw-1)*W_step
         el=(w**2-mp**2)/2./mp
         do lm=0,3
            do is=1,3
               do ii=1,6
                  Rint(lm,is,ii,jj)=gauss(Fint,-1.,1.,IntAccur)
               enddo
            enddo
            do is=1,3
               call re_mult()
            enddo
            lun=(lm+2)*10
            write(lun+1,201)W,El,(Mlp(lm,is,jj),is=1,3)
            write(lun+2,201)W,El,(Mlm(lm,is,jj),is=1,3)
            write(lun+3,201)W,El,(Elp(lm,is,jj),is=1,3)
            write(lun+4,201)W,El,(Elm(lm,is,jj),is=1,3)
            write(lun+5,201)W,El,(Slp(lm,is,jj),is=1,3)
            write(lun+6,201)W,El,(Slm(lm,is,jj),is=1,3)
         enddo
      EndDo
 201  format(2x,2F7.3,2x,3F12.5)
      do lm=0,3
         do is=1,6
            lun=(lm+2)*10+is
            close(lun)
         enddo
      enddo   
      enddo   

      return
      end
      
