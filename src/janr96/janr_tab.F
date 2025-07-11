      SUBROUTINE JANR_TAB(Ierr)
c     
c     _begin_doc
c     RCS ID string
c     $Id$
c     
c     Documentation for subroutine JANR_TAB
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
c     Author:  S.Stepanyan   Created:  Wed Jun 23 17:16:39 EST 1999
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
#include "inpar.inc"
#include "cross_sec.inc"
#include "janrtab.inc"
#include "fit_par.inc"
#include "res_func.inc"
c     _end_inc
c     
c     _begin_var
c     input/output variables:
c     -----------------------
      Integer Ierr
c     
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
      PARAMETER (CRNAME='JANR_TAB')
      PARAMETER (CRAUTH='Stepanyan')
c     
c     Local User defined variables:
c     -----------------------------
      character*10 fnamm3(38)
      data fnamm3/'jre0p.dat','jrs0p.dat','jrm1m.dat','jrs1m.dat',
     +     'jrm1p.dat','jre1p.dat','jrs1p.dat',
     +     'jrm2m.dat','jre2m.dat','jrs2m.dat',
     +     'jrm2p.dat','jre2p.dat','jrs2p.dat',
     +     'jrm3m.dat','jre3m.dat','jrs3m.dat',
     +     'jrm3p.dat','jre3p.dat','jrs3p.dat',
     +     'jie0p.dat','jis0p.dat','jim1m.dat','jis1m.dat',
     +     'jim1p.dat','jie1p.dat','jis1p.dat',
     +     'jim2m.dat','jie2m.dat','jis2m.dat',
     +     'jim2p.dat','jie2p.dat','jis2p.dat',
     +     'jim3m.dat','jie3m.dat','jis3m.dat',
     +     'jim3p.dat','jie3p.dat','jis3p.dat'/
c      data fnamm3/'Re0p.dat','Rs0p.dat','Rm1m.dat','Rs1m.dat',
c     +     'Rm1p.dat','Re1p.dat','Rs1p.dat',
c     +     'Rm2m.dat','Re2m.dat','Rs2m.dat',
c     +     'Rm2p.dat','Re2p.dat','Rs2p.dat',
c     +     'Rm3m.dat','Re3m.dat','Rs3m.dat',
c     +     'Rm3p.dat','Re3p.dat','Rs3p.dat',
c     +     'Ie0p.dat','Is0p.dat','Im1m.dat','Is1m.dat',
c     +     'Im1p.dat','Ie1p.dat','Is1p.dat',
c     +     'Im2m.dat','Ie2m.dat','Is2m.dat',
c     +     'Im2p.dat','Ie2p.dat','Is2p.dat',
c     +     'Im3m.dat','Ie3m.dat','Is3m.dat',
c     +     'Im3p.dat','Ie3p.dat','Is3p.dat'/
c     
      INTEGER iw,iq2,icos,iphi,i,j,ifile,lun
      INTEGER Nw,Nq2,Ncos,Nphi
      REAL Q2, W, CosTh, Phi, Epsilon, nu, ep, mp
      REAL mpion,Egamma,Epion,qpion,kgamma,transfer
      DATA mp/0.938/
c     _end_var
c     
c     executable code for routine JANR_RUN:
c---- 6----------------------------------------------------------------72
c     Get the number of data points in the table
      Ierr = 0
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
      IF ( Wstep.GT.0 .AND. Wmin+Wstep.LT.Wmax )THEN
         Nw = INT( (Wmax - Wmin) / Wstep ) + 1
      ELSE
         If( Wmin.Eq.Wmax .or. Wstep.Eq.0. )Then
            Nw=1
         Else
            Nw = INT( (Wmax - Wmin) / abs(Wstep) )
         EndIf
      ENDIF
c     
      IF ( cosstep.GT.0 .AND. cosmin+cosstep.LT.cosmax )THEN
         Ncos = INT( (cosmax - cosmin) / cosstep ) + 1
      ELSE
         If( cosmin.Eq.cosmax )Then
            Ncos=1
         EndIf
      ENDIF
c     
      IF ( phistep.GT.0 .AND. phimin+phistep.LT.phimax )THEN
         Nphi = INT( (phimax - phimin) / phistep ) + 1
      ELSE
         If( phimin.Eq.phimax )Then
            Nphi=1
         EndIf
      ENDIF
c-
c-
c     
      OPEN (UNIT=57,FILE=tabfile,STATUS='unknown')
c     Make additional printouts if flag is set, this set of printouts can be changed.
      if(fitopt.eq.'GP'.or.fitopt.eq.'gp')then
c     open files for MM3
         do i=1,38
            ifile=i
            Lun=60+ifile
            open(unit=lun,file=fnamm3(ifile),status='unknown')
         enddo
      elseif(fitopt.eq.'GR'.or.fitopt.eq.'gr')then
         do i=1,2
            open(unit=60+i,file='resp_f1.dat',status='unknown')
            write(60+i,*)'Q2,w,costh,phi,
     +sigmaTLP,sigmaTy,sigmaLy,sigmaTTx,sigmaTTy,
     +sigmaTTz,sigmaTLx,sigmaTLy,sigmaTLz,
     +sigmaTTPx,sigmaTTPz,sigmaTLPx,sigmaTLPy,sigmaTLPz'
         enddo
      endif
c     
         Q2 = (Q2min + Q2max)/2.
         write(6,*)nw,ncos,nphi,'  at Q2=',Q2
c     
         Do iw=1,Nw
c     Loop over W
            W = Wmin + float(iw-1)*Wstep
c     Get the epsilon
            nu=(w**2+q2-mp**2)/2./mp
            ep=ebeam-nu
            epsilon=1./(1.+2.*(q2+nu**2)/(4.*ebeam*ep-q2))
            do icos=1,Ncos
c     Loop over cos(Theta*)
               costh=cosmin + float(icos-1)*cosstep
               dO iphi=1,Nphi
c     Loop over phi*
                  phi=phimin + float(iphi-1)*phistep
c     Calculate cross section
c                  write(6,*) 'before', phi , phistep
                  CALL janr_RUN(W,Q2,epsilon,costh,phi,ierr)                  
c                  write(6,*) 'after', phi , phistep
	mpion=0.139
      Egamma=(w**2-mp**2-Q2)/w/2.
      Epion=(w**2-mp**2+mpion**2)/w/2.
      qpion=sqrt(Epion**2-mpion**2)
      kgamma=sqrt((w-Egamma)**2-mp**2)
      transfer=mpion**2-Q2+2.*(qpion*kgamma*costh-Epion*Egamma)
      
c      cs_tab(iq2,iw,icos,iphi) =3.14*sigma(pi_type2fit)/(qpion*kgamma)
      cs_tab(iq2,iw,icos,iphi) =sigma(pi_type2fit)
                  
c                  write(57,*)q2,w,epsilon,costh,phi,cs_tab(iq2,iw,icos,iphi)
       write(57,201)w,costh,cs_tab(iq2,iw,icos,iphi)
c     write(6,*)q2,w,epsilon,costh,phi,cs_tab(iq2,iw,icos,iphi)
                  if(fitopt.eq.'GR'.or.fitopt.eq.'gr')then
                     do i=1,2
                        write(60+i,*)q2,w,costh,phi,
     +                        sigmaTLP(i),sigmaTy(i),sigmaLy(i),sigmaTTx(i),sigmaTTy(i),
     +                       sigmaTTz(i),sigmaTLx(i),sigmaTLy(i),sigmaTLz(i),
     +                       sigmaTTPx(i),sigmaTTPz(i),sigmaTLPx(i),sigmaTLPy(i),sigmaTLPz(i)
                     enddo
                  endif
               enDdO
            enddo
         EndDo
         close(57)
 201  format(2x,F7.3,2x,2F12.4,2x,F12.3)
         RETURN
 999     continue
         write(6,*)'Error while reading the data file'
         RETURN
         END
c
c------------------------------------------------------------------------------




