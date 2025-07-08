      SUBROUTINE MULT_AMPL (S,Q2,kgamma,qpion)
c     
c     _begin_doc
c     RCS ID string
c     $Id$
c     
c     Documentation for subroutine MULT_AMPL
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
c     Author:         Created:  Mon Dec 16 10:31:34 EST 1996
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
#include "gpar.inc"
#include "back_ampl.inc"
#include "multipoles.inc"
#include "amplitudes.inc"
#include "fit_par.inc"
c     _end_inc
c     
c     _begin_var
c     input/output variables:
c     -----------------------
      REAL S, Q2,kgamma,qpion
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
c     CHARACTER*100  CRMESS
      PARAMETER (CRNAME='MULT_AMPL')
      PARAMETER (CRAUTH='')
c     
c     Local User defined variables:
c     -----------------------------
      INTEGER NWRITE
      DATA NWRITE/0/
      REAL Br1(2,9),Br3(2,8)
      REAL BM33,BE23,B3M,B3Mn,B6E,B6M,B6Mn
      REAL BM11,w1
      REAL BrM1(2,2,9),BrE1(2,2,9),BrS1(2,2,9)
      REAL BrM3(2,8),BrE3(2,8),BrS3(2,8)
      REAL RMP33,IMP33,coeff
      REAL YOUTmp33r,YOUTmp33i
      REAL XIN,elab,WW,Suppress
      REAL YOUTSimr,YOUTSimi,YOUTPimr,YOUTPimi
      REAL YOUTp11r,YOUTs11r,YOUTd13r,YOUTf15r,YOUTp13r,YOUTd15r
      REAL YOUTp11i,YOUTs11i,YOUTd13i,YOUTf15i,YOUTp13i,YOUTd15i
      REAL YOUTp33r,YOUTd33r,YOUTf37r,YOUTs31r,YOUTp31r
      REAL YOUTp33i,YOUTd33i,YOUTf37i,YOUTs31i,YOUTp31i
      INTEGER k,i,j,lun
c     _end_var
c     
c     executable code for routine MULT_AMPL:
c---- 6---------------------------------------------------------------
c     -72
c     
      CALL breit_ampl(S,Q2,kgamma,qpion)
      CALL background(S,Q2)
      WW=sqrt(s)
      XIN=WW
      elab=(s-mn**2)/2./mn
c      if (WW.ge.1.34) go to 366
      if (WW.le.1.34) go to 366
	go to 367
366      CALL interp1(WSa,simr,13,3,XIN,YOUTSimr)
      CALL interp1(WSa,simi,13,3,XIN,YOUTSimi)
c      CALL interp1(WSa,pimr,13,3,XIN,YOUTPimr)
c      CALL interp1(WSa,pimi,13,3,XIN,YOUTPimi)
      do j=1,2
         BreitE1(2,j,2)=BreitE1(1,j,2)*YOUTSimr/(1.-YOUTSimi)
         BreitS1(2,j,2)=BreitS1(1,j,2)*YOUTSimr/(1.-YOUTSimi)
c         BreitM1(2,j,1)=BreitM1(1,j,1)*YOUTPimr/(1.-YOUTPimi)
c         BreitS1(2,j,1)=BreitS1(1,j,1)*YOUTPimr/(1.-YOUTPimi)
      enddo
367   continue
      do i=1,2
c     -construct the real(i=1) and imaginary(i=2) parts of the
c     Breit_Wigner contributions into multipoles
         do j=1,2
c     -construct these contributions for combinations pA(1/2)(j=1),
c     nA(1/2)(j=2)
            BrM1(i,j,6)=BreitM1(i,j,6)+BreitM1(i,j,13)
            BrE1(i,j,6)=BreitE1(i,j,6)+BreitE1(i,j,13)
            BrS1(i,j,6)=BreitS1(i,j,6)+BreitS1(i,j,13)
            BrM1(i,j,7)=BreitM1(i,j,7)+BreitM1(i,j,11)    
            BrE1(i,j,7)=BreitE1(i,j,7)+BreitE1(i,j,11)    
            BrS1(i,j,7)=BreitS1(i,j,7)+BreitS1(i,j,11)    
            BrE1(i,j,2)=BreitE1(i,j,2)+BreitE1(i,j,5)
            BrS1(i,j,2)=BreitS1(i,j,2)+BreitS1(i,j,5)
            BrM1(i,j,1)=BreitM1(i,j,1)+BreitM1(i,j,10)
            BrS1(i,j,1)=BreitS1(i,j,1)+BreitS1(i,j,10)
            BrM1(i,j,3)=BreitM1(i,j,3)+BreitM1(i,j,4)+BreitM1(i,j,12)
            BrE1(i,j,3)=BreitE1(i,j,3)+BreitE1(i,j,4)+BreitE1(i,j,12)
            BrS1(i,j,3)=BreitS1(i,j,3)+BreitS1(i,j,4)+BreitS1(i,j,12)
            BrM1(i,j,9)=BreitM1(i,j,9)
            BrE1(i,j,9)=BreitE1(i,j,9)
            BrS1(i,j,9)=BreitS1(i,j,9)
            BrM1(i,j,8)=BreitM1(i,j,8)
            BrE1(i,j,8)=BreitE1(i,j,8)
            BrS1(i,j,8)=BreitS1(i,j,8)
         enddo
      enddo
c     
c 1+ P33 M1+
      CALL interp(WSa,mp33r,51,3,XIN,YOUTmp33r)
      RMP33=YOUTmp33r
      CALL interp(WSa,mp33i,51,3,XIN,YOUTmp33i)
      IMP33=YOUTmp33i
c 1- P11 1(1)
      CALL interp(WSa,p11r,51,3,XIN,YOUTp11r)
      Br1(1,1)=YOUTp11r
      CALL interp(WSa,p11i,51,3,XIN,YOUTp11i)
      Br1(2,1)=YOUTp11i
c 0+ S11 1(2)
      CALL interp(WSa,s11r,51,3,XIN,YOUTs11r)
      Br1(1,2)=YOUTs11r
      CALL interp(WSa,s11i,51,3,XIN,YOUTs11i)
      Br1(2,2)=YOUTs11i
c 2- D13 1(3)
      CALL interp(WSa,d13r,51,3,XIN,YOUTd13r)
      Br1(1,3)=YOUTd13r
      CALL interp(WSa,d13i,51,3,XIN,YOUTd13i)
      Br1(2,3)=YOUTd13i
      Br1(1,4)=0.
      Br1(2,4)=0.
      Br1(1,5)=0.
      Br1(2,5)=0.
      Br1(1,8)=0.
      Br1(2,8)=0.
c 3- F15 1(6)
      CALL interp(WSa,f15r,51,3,XIN,YOUTf15r)
      Br1(1,6)=YOUTf15r
      CALL interp(WSa,f15i,51,3,XIN,YOUTf15i)
      Br1(2,6)=YOUTf15i
c 1+ P13 1(7)
      CALL interp(WSa,p13r,51,3,XIN,YOUTp13r)
      Br1(1,7)=YOUTp13r
      CALL interp(WSa,p13i,51,3,XIN,YOUTp13i)
      Br1(2,7)=YOUTp13i
c 2+ D15 1(9)
      CALL interp(WSa,d15r,51,3,XIN,YOUTd15r)
      Br1(1,9)=YOUTd15r
      CALL interp(WSa,d15i,51,3,XIN,YOUTd15i)
      Br1(2,9)=YOUTd15i

c 1+ P33 3(1)
      CALL interp(WSa,p33r,51,3,XIN,YOUTp33r)
      Br3(1,1)=YOUTp33r
      CALL interp(WSa,p33i,51,3,XIN,YOUTp33i)
      Br3(2,1)=YOUTp33i
c 2- D33 3(2)
      CALL interp(WSa,d33r,51,3,XIN,YOUTd33r)
      Br3(1,2)=YOUTd33r
      CALL interp(WSa,d33i,51,3,XIN,YOUTd33i)
      Br3(2,2)=YOUTd33i
c 3+ F37 3(3)
      CALL interp(WSa,f37r,51,3,XIN,YOUTf37r)
      Br3(1,3)=YOUTf37r
      CALL interp(WSa,f37i,51,3,XIN,YOUTf37i)
      Br3(2,3)=YOUTf37i
      Br3(1,4)=0.
      Br3(2,4)=0.
c 0+ S31 3(5)
      CALL interp(WSa,s31r,51,3,XIN,YOUTs31r)
      Br3(1,5)=YOUTs31r
      CALL interp(WSa,s31i,51,3,XIN,YOUTs31i)
      Br3(2,5)=YOUTs31i
c 3- F35 3(6)
      Br3(1,6)=0.
      Br3(2,6)=0.
c 2+ D35 3(7)
      Br3(1,7)=0.
      Br3(2,7)=0.
c 1- P31 3(8)
      CALL interp(WSa,p31r,51,3,XIN,YOUTp31r)
      Br3(1,8)=YOUTp31r
      CALL interp(WSa,p31i,51,3,XIN,YOUTp31i)
      Br3(2,8)=YOUTp31i
      do j=1,9
         do i=1,2
c     - construct Real parts of multipole amplitudes pA(1/2) and nA(1/2)
            MM1(1,i,j)=BM1(i,j)*(1.-Br1(2,j))+BrM1(1,i,j)
            ME1(1,i,j)=BE1(i,j)*(1.-Br1(2,j))+BrE1(1,i,j)
            MS1(1,i,j)=BS1(i,j)*(1.-Br1(2,j))+BrS1(1,i,j)
c     - construct Imaginary parts of multipole amplitudes pA(1/2) and nA(1/2) 
c     
            MM1(2,i,j)=BM1(i,j)*Br1(1,j)+BrM1(2,i,j)
            ME1(2,i,j)=BE1(i,j)*Br1(1,j)+BrE1(2,i,j)
            MS1(2,i,j)=BS1(i,j)*Br1(1,j)+BrS1(2,i,j)
         enddo
      enddo
c
      do i=1,2
c     -construct the real and imaginary parts of the Breit_Wigner 
c     -contributions into multipoles with I=3/2
         BrM3(i,1)=BreitM3(i,1)+BreitM3(i,4)+BreitM3(i,10)
         BrE3(i,1)=BreitE3(i,1)+BreitE3(i,4)+BreitE3(i,10)
         BrS3(i,1)=BreitS3(i,1)+BreitS3(i,4)+BreitS3(i,10)
         BrM3(i,2)=BreitM3(i,2)+BreitM3(i,12)
         BrE3(i,2)=BreitE3(i,2)+BreitE3(i,12)
         BrS3(i,2)=BreitS3(i,2)+BreitS3(i,12)
         BrM3(i,3)=BreitM3(i,3)
         BrE3(i,3)=BreitE3(i,3)
         BrS3(i,3)=BreitS3(i,3)
         BrE3(i,5)=BreitE3(i,5)+BreitE3(i,9)
         BrS3(i,5)=BreitS3(i,5)+BreitS3(i,9)
         BrM3(i,6)=BreitM3(i,6)+BreitM3(i,13)      
         BrE3(i,6)=BreitE3(i,6)+BreitE3(i,13)      
         BrS3(i,6)=BreitS3(i,6)+BreitS3(i,13)      
         BrM3(i,7)=BreitM3(i,7)
         BrE3(i,7)=BreitE3(i,7)
         BrS3(i,7)=BreitS3(i,7)
         BrM3(i,8)=BreitM3(i,8)+BreitM3(i,11)     
         BrS3(i,8)=BreitS3(i,8)+BreitS3(i,11)
      enddo     
    	   Suppress=(WW/1.2286)**7.2
c     - construct Real parts of (1+) multipole amplitudes A(3/2)
      MM3(1,1)=BM3(1)*(1.-Br3(2,1))+BrM3(1,1)*Suppress
      ME3(1,1)=BE3(1)*(1.-Br3(2,1))+BrE3(1,1)
      MS3(1,1)=BS3(1)*(1.-Br3(2,1))+BrS3(1,1)
c     - construct Imaginary parts of (1+) multipole amplitudes A(3/2)
      MM3(2,1)=BM3(1)*Br3(1,1)+BrM3(2,1)*Suppress  
      ME3(2,1)=BE3(1)*Br3(1,1)+BrE3(2,1)
      MS3(2,1)=BS3(1)*Br3(1,1)+BrS3(2,1)
      if (WW.ge.1.2283) go to 666
	go to 667
666    	    MM3(1,1)=cm1*RMP33
	    MM3(2,1)=cm1*IMP33
667   continue

c     - construct Real parts of other multipole amplitudes A(3/2) 
      do j=2,8
         MM3(1,j)=BM3(j)*(1.-Br3(2,j))+BrM3(1,j)
         ME3(1,j)=BE3(j)*(1.-Br3(2,j))+BrE3(1,j)
         MS3(1,j)=BS3(j)*(1.-Br3(2,j))+BrS3(1,j)
c     - construct Imaginary parts of other multipole amplitudes A(3/2) 
         MM3(2,j)=BM3(j)*Br3(1,j)+BrM3(2,j)
         ME3(2,j)=BE3(j)*Br3(1,j)+BrE3(2,j)
         MS3(2,j)=BS3(j)*Br3(1,j)+BrS3(2,j)
      enddo
CCCC for calculating multipoles in mFm
c	coeff=10.
CCCCC for calculating observables
	coeff=1.
      do j=1,9
         do i=1,2
         do k=1,2
            MM1(k,i,j)=coeff*MM1(k,i,j)
            ME1(k,i,j)=coeff*ME1(k,i,j)
            MS1(k,i,j)=coeff*MS1(k,i,j)
         enddo
      enddo
         enddo
      do j=1,8
         do k=1,2
            MM3(k,j)=coeff*MM3(k,j)
            ME3(k,j)=coeff*ME3(k,j)
            MS3(k,j)=coeff*MS3(k,j)
         enddo
      enddo
c      do j=1,2
c	ME3(j,5)=ME3(j,5)/1.5+ME1(j,1,2)
c	MS3(j,5)=MS3(j,5)/1.5+MS1(j,1,2)
c	MM3(j,8)=MM3(j,8)/1.5+MM1(j,1,1)
c	MS3(j,8)=MS3(j,8)/1.5+MS1(j,1,1)
c	MM3(j,1)=MM3(j,1)/1.5+MM1(j,1,7)
c	ME3(j,1)=ME3(j,1)/1.5+ME1(j,1,7)
c	MS3(j,1)=MS3(j,1)/1.5+MS1(j,1,7)
c	MM3(j,2)=MM3(j,2)/1.5+MM1(j,1,3)
c	ME3(j,2)=ME3(j,2)/1.5+ME1(j,1,3)
c	MS3(j,2)=MS3(j,2)/1.5+MS1(j,1,3)
c	MM3(j,7)=MM3(j,7)/1.5+MM1(j,1,9)
c	ME3(j,7)=ME3(j,7)/1.5+ME1(j,1,9)
c	MS3(j,7)=MS3(j,7)/1.5+MS1(j,1,9)
c	MM3(j,6)=MM3(j,6)/1.5+MM1(j,1,6)
c	ME3(j,6)=ME3(j,6)/1.5+ME1(j,1,6)
c	MS3(j,6)=MS3(j,6)/1.5+MS1(j,1,6)
c	MM3(j,3)=MM3(j,3)/1.5+MM1(j,1,8)
c	ME3(j,3)=ME3(j,3)/1.5+ME1(j,1,8)
c	MS3(j,3)=MS3(j,3)/1.5+MS1(j,1,8)
c      enddo   

      if(fitopt.eq.'GP'.or.fitopt.eq.'gp')then
c         write(61,332)sqrt(s),ME3(1,5),ME3(2,5)
c         write(62,332)sqrt(s),MS3(1,5),MS3(2,5)
c         write(63,332)sqrt(s),MM3(1,8),MM3(2,8)
c         write(64,332)sqrt(s),MS3(1,8),MS3(2,8)
c         write(65,332)sqrt(s),MM3(1,1),MM3(2,1)
c         write(66,332)sqrt(s),ME3(1,1),ME3(2,1)
c         write(67,332)sqrt(s),MS3(1,1),MS3(2,1)
c         write(68,332)sqrt(s),MM3(1,2),MM3(2,2)
c         write(69,332)sqrt(s),ME3(1,2),ME3(2,2)
c         write(70,332)sqrt(s),MS3(1,2),MS3(2,2)
c         write(71,332)sqrt(s),MM3(1,7),MM3(2,7)
c         write(72,332)sqrt(s),ME3(1,7),ME3(2,7)
c         write(73,332)sqrt(s),MS3(1,7),MS3(2,7)
c         write(74,332)sqrt(s),MM3(1,6),MM3(2,6)
c         write(75,332)sqrt(s),ME3(1,6),ME3(2,6)
c         write(76,332)sqrt(s),MS3(1,6),MS3(2,6)
c         write(77,332)sqrt(s),MM3(1,3),MM3(2,3)
c         write(78,332)sqrt(s),ME3(1,3),ME3(2,3)
c         write(79,332)sqrt(s),MS3(1,3),MS3(2,3)

	w1=ww
         write(61,332)w1,elab,ME3(1,5),ME1(1,1,2),ME1(1,2,2)
         write(62,332)w1,elab,MS3(1,5),MS1(1,1,2),MS1(1,2,2)
         write(63,332)w1,elab,MM3(1,8),MM1(1,1,1),MM1(1,2,1)
         write(64,332)w1,elab,MS3(1,8),MS1(1,1,1),MS1(1,2,1)
         write(65,332)w1,elab,MM3(1,1),MM1(1,1,7),MM1(1,2,7)
         write(66,332)w1,elab,ME3(1,1),ME1(1,1,7),ME1(1,2,7)
         write(67,332)w1,elab,MS3(1,1),MS1(1,1,7),MS1(1,2,7)
         write(68,332)w1,elab,MM3(1,2),MM1(1,1,3),MM1(1,2,3)
         write(69,332)w1,elab,ME3(1,2),ME1(1,1,3),ME1(1,2,3)
         write(70,332)w1,elab,MS3(1,2),MS1(1,1,3),MS1(1,2,3)
         write(71,332)w1,elab,MM3(1,7),MM1(1,1,9),MM1(1,2,9)
         write(72,332)w1,elab,ME3(1,7),ME1(1,1,9),ME1(1,2,9)
         write(73,332)w1,elab,MS3(1,7),MS1(1,1,9),MS1(1,2,9)
         write(74,332)w1,elab,MM3(1,6),MM1(1,1,6),MM1(1,2,6)
         write(75,332)w1,elab,ME3(1,6),ME1(1,1,6),ME1(1,2,6)
         write(76,332)w1,elab,MS3(1,6),MS1(1,1,6),MS1(1,2,6)
         write(77,332)w1,elab,MM3(1,3),MM1(1,1,8),MM1(1,2,8)
         write(78,332)w1,elab,ME3(1,3),ME1(1,1,8),ME1(1,2,8)
         write(79,332)w1,elab,MS3(1,3),MS1(1,1,8),MS1(1,2,8)
         write(80,332)w1,elab,ME3(2,5),ME1(2,1,2),ME1(2,2,2)
         write(81,332)w1,elab,MS3(2,5),MS1(2,1,2),MS1(2,2,2)
         write(82,332)w1,elab,MM3(2,8),MM1(2,1,1),MM1(2,2,1)
         write(83,332)w1,elab,MS3(2,8),MS1(2,1,1),MS1(2,2,1)
         write(84,332)w1,elab,MM3(2,1),MM1(2,1,7),MM1(2,2,7)
         write(85,332)w1,elab,ME3(2,1),ME1(2,1,7),ME1(2,2,7)
         write(86,332)w1,elab,MS3(2,1),MS1(2,1,7),MS1(2,2,7)
         write(87,332)w1,elab,MM3(2,2),MM1(2,1,3),MM1(2,2,3)
         write(88,332)w1,elab,ME3(2,2),ME1(2,1,3),ME1(2,2,3)
         write(89,332)w1,elab,MS3(2,2),MS1(2,1,3),MS1(2,2,3)
         write(90,332)w1,elab,MM3(2,7),MM1(2,1,9),MM1(2,2,9)
         write(91,332)w1,elab,ME3(2,7),ME1(2,1,9),ME1(2,2,9)
         write(92,332)w1,elab,MS3(2,7),MS1(2,1,9),MS1(2,2,9)
         write(93,332)w1,elab,MM3(2,6),MM1(2,1,6),MM1(2,2,6)
         write(94,332)w1,elab,ME3(2,6),ME1(2,1,6),ME1(2,2,6)
         write(95,332)w1,elab,MS3(2,6),MS1(2,1,6),MS1(2,2,6)
         write(96,332)w1,elab,MM3(2,3),MM1(2,1,8),MM1(2,2,8)
         write(97,332)w1,elab,ME3(2,3),ME1(2,1,8),ME1(2,2,8)
         write(98,332)w1,elab,MS3(2,3),MS1(2,1,8),MS1(2,2,8)
      endif
 332  format(2x,2f10.3,3f12.4)
c 332  format(2x,2f10.3,2f12.4)
c 332  format(2x,f10.3,2f12.4)
c     
      RETURN
      END
