      SUBROUTINE MULT_AMPL (S,Q2,kgamma,qpion)

      IMPLICIT NONE
      SAVE

#include "inpar.inc"
#include "gpar.inc"
#include "back_ampl.inc"
#include "multipoles.inc"
#include "amplitudes.inc"
#include "fit_par.inc"
#include "highmult.inc"

      REAL S, Q2,kgamma,qpion

      INTEGER NWRITE
      REAL WWLAST
      DATA NWRITE/0/wWlast/0/
      REAL Br1(2,9),Br3(2,8)
      REAL BM33,BE23,B3M,B3Mn,B6E,B6M,B6Mn
      REAL BM11,w1
      REAL BrM1(2,2,9),BrE1(2,2,9),BrS1(2,2,9)
      REAL BrM3(2,8),BrE3(2,8),BrS3(2,8)
      REAL RMP33,IMP33,coeff
      REAL YOUTmp33r,YOUTmp33i
      REAL XIN,XIN1,elab,WW,Suppress
      REAL YOUTSimr,YOUTSimi,YOUTPimr,YOUTPimi
      REAL YOUTp11r,YOUTs11r,YOUTd13r,YOUTf15r,YOUTp13r,YOUTd15r
      REAL YOUTp11i,YOUTs11i,YOUTd13i,YOUTf15i,YOUTp13i,YOUTd15i
      REAL YOUTp33r,YOUTd33r,YOUTf37r,YOUTs31r,YOUTp31r
      REAL YOUTp33i,YOUTd33i,YOUTf37i,YOUTs31i,YOUTp31i
      REAL mpion(2),qpi(2),ratio
      INTEGER k,i,j,lun
    
      CALL breit_ampl(S,Q2,kgamma,qpion)
      CALL background(S,Q2)
      CALL back_high(S,Q2)
      
      WW = sqrt(s)
      mpion(1) = 0.13498
      mpion(2) = 0.13957  
      
      do i=1,2
        qpi(i)=sqrt(((s-mn22+mpion(i)**2)/2./WW)**2-mpion(i)**2)
      enddo
      
      ratio = qpi(2)/qpi(1)
      
      XIN   = WW
      elab  = (s-mn**2)/2./mn
      
      if (WW.le.1.34) go to 366
      
      go to 367
      
366   CALL interp1(WSa,simr,13,3,XIN,YOUTSimr)
      CALL interp1(WSa,simi,13,3,XIN,YOUTSimi)
      
      do j=1,2
         BreitE1(2,j,2)=BreitE1(1,j,2)*YOUTSimr/(1.-YOUTSimi)
         BreitS1(2,j,2)=BreitS1(1,j,2)*YOUTSimr/(1.-YOUTSimi)
      enddo
      
367   continue

c Breit-Wigner amplitudes
c I=1/2
      
      do i=1,2 		! real(1) imaginary(2)
         do j=1,2	! pA(1/2) (j=1) nA(1/2) (j=2)
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
      
c 1+ P33 M1+
      XIN1=XIN
      CALL interp(WSa,mp33r,51,4,XIN1,YOUTmp33r)
      RMP33=YOUTmp33r
      CALL interp(WSa,mp33i,51,4,XIN1,YOUTmp33i)
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
      Br3(1,1)=YOUTp33r-0.05
      CALL interp(WSa,p33i,51,3,XIN,YOUTp33i)
      Br3(2,1)=YOUTp33i*1.017
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
         do i=1,2	! Real (1) Imag (2)
            MM1(1,i,j)=BM1(i,j)*(1.-Br1(2,j))+BrM1(1,i,j)
            ME1(1,i,j)=BE1(i,j)*(1.-Br1(2,j))+BrE1(1,i,j)
            MS1(1,i,j)=BS1(i,j)*(1.-Br1(2,j))+BrS1(1,i,j)
            MM1(2,i,j)=BM1(i,j)*Br1(1,j)+BrM1(2,i,j)
            ME1(2,i,j)=BE1(i,j)*Br1(1,j)+BrE1(2,i,j)
            MS1(2,i,j)=BS1(i,j)*Br1(1,j)+BrS1(2,i,j)
         enddo
      enddo
      
c I=3/2

      do i=1,2
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
           
c     - construct Re(1) Im(2) parts of (1+) multipole amplitudes A(3/2)
      MM3(1,1)=BM3(1)*(1.-Br3(2,1))+BrM3(1,1) 
      ME3(1,1)=BE3(1)*(1.-Br3(2,1))+BrE3(1,1)
      MS3(1,1)=BS3(1)*(1.-Br3(2,1))+BrS3(1,1)
      MM3(2,1)=BM3(1)*Br3(1,1)+BrM3(2,1)  
      ME3(2,1)=BE3(1)*Br3(1,1)+BrE3(2,1)
      MS3(2,1)=BS3(1)*Br3(1,1)+BrS3(2,1)
      
      MM3(1,1)=resM*RMP33*(WW/1.2283)**(0.3857+0.229-0.095)
      MM3(2,1)=resM*IMP33*(WW/1.2283)**(0.3857+0.229-0.095)
      
      if (WW.ge.1.2283) then
        MM3(1,1)=resM*RMP33*(WW/1.2283)**(-0.096546+0.18467+0.0932)
	MM3(2,1)=resM*IMP33*(WW/1.2283)**(-1.6189+0.39945+0.192)
    	ME3(1,1)= ME3(1,1)*(WW/1.2283)**(-1.465+1.55-0.54)
    	ME3(2,1)= ME3(2,1)*(WW/1.2283)**(5.505-0.66997+0.29)
    	MS3(1,1)= MS3(1,1)*(WW/1.2283)**(2.3835)
    	MS3(2,1)= MS3(2,1)*(WW/1.2283)**(-16.621-2.495)
      endif

c     - construct Re(1) Im(2) parts of other multipole amplitudes A(3/2) 
      do j=2,8
         MM3(1,j)=BM3(j)*(1.-Br3(2,j))+BrM3(1,j)
         ME3(1,j)=BE3(j)*(1.-Br3(2,j))+BrE3(1,j)
         MS3(1,j)=BS3(j)*(1.-Br3(2,j))+BrS3(1,j)
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
      
C	ReacM1(i,j)=0+	i=1,2 for Re,Im;	j=1,2 for pi0,pi+	 
C	ReacM2(i,j)=1-
C	ReacM3(i,j)=1+
C	ReacM4(i,j)=2-
C	ReacM5(i,j)=2+
C	ReacM6(i,j)=3-
C	ReacM7(i,j)=3+

      do i=1,2
        ReacE1(i,1)= ME1(i,1,2)+ME3(i,5)/1.5
        ReacE1(i,2)= sqrt(2.)*(ME1(i,1,2)-ME3(i,5)/3.)
        ReacS1(i,1)= MS1(i,1,2)+MS3(i,5)/1.5
        ReacS1(i,2)= sqrt(2.)*(MS1(i,1,2)-MS3(i,5)/3.)
        ReacM2(i,1)= MM1(i,1,1)+MM3(i,8)/1.5
        ReacM2(i,2)= sqrt(2.)*(MM1(i,1,1)-MM3(i,8)/3.)*ratio
        ReacS2(i,1)= MS1(i,1,1)+MS3(i,8)/1.5
        ReacS2(i,2)= sqrt(2.)*(MS1(i,1,1)-MS3(i,8)/3.)*ratio
        ReacM3(i,1)= MM1(i,1,7)+MM3(i,1)/1.5
        ReacM3(i,2)= sqrt(2.)*(MM1(i,1,7)-MM3(i,1)/3.)*ratio
        ReacE3(i,1)= ME1(i,1,7)+ME3(i,1)/1.5
        ReacE3(i,2)= sqrt(2.)*(ME1(i,1,7)-ME3(i,1)/3.)*ratio
        ReacS3(i,1)= MS1(i,1,7)+MS3(i,1)/1.5
        ReacS3(i,2)= sqrt(2.)*(MS1(i,1,7)-MS3(i,1)/3.)*ratio
        ReacM4(i,1)= MM1(i,1,3)+MM3(i,2)/1.5
        ReacM4(i,2)= sqrt(2.)*(MM1(i,1,3)-MM3(i,2)/3.)*ratio**2
        ReacE4(i,1)= ME1(i,1,3)+ME3(i,2)/1.5
        ReacE4(i,2)= sqrt(2.)*(ME1(i,1,3)-ME3(i,2)/3.)*ratio**2
        ReacS4(i,1)= MS1(i,1,3)+MS3(i,2)/1.5
        ReacS4(i,2)= sqrt(2.)*(MS1(i,1,3)-MS3(i,2)/3.)*ratio**2
        ReacM5(i,1)= MM1(i,1,9)+MM3(i,7)/1.5
        ReacM5(i,2)= sqrt(2.)*(MM1(i,1,9)-MM3(i,7)/3.)*ratio**2
        ReacE5(i,1)= ME1(i,1,9)+ME3(i,7)/1.5
        ReacE5(i,2)= sqrt(2.)*(ME1(i,1,9)-ME3(i,7)/3.)*ratio**2
        ReacS5(i,1)= MS1(i,1,9)+MS3(i,7)/1.5
        ReacS5(i,2)= sqrt(2.)*(MS1(i,1,9)-MS3(i,7)/3.)*ratio**2
        ReacM6(i,1)= MM1(i,1,6)+MM3(i,6)/1.5
        ReacM6(i,2)= sqrt(2.)*(MM1(i,1,6)-MM3(i,6)/3.)*ratio**3
        ReacE6(i,1)= ME1(i,1,6)+ME3(i,6)/1.5
        ReacE6(i,2)= sqrt(2.)*(ME1(i,1,6)-ME3(i,6)/3.)*ratio**3
        ReacS6(i,1)= MS1(i,1,6)+MS3(i,6)/1.5
        ReacS6(i,2)= sqrt(2.)*(MS1(i,1,6)-MS3(i,6)/3.)*ratio**3
        ReacM7(i,1)= MM1(i,1,8)+MM3(i,3)/1.5
        ReacM7(i,2)= sqrt(2.)*(MM1(i,1,8)-MM3(i,3)/3.)*ratio**3
        ReacE7(i,1)= ME1(i,1,8)+ME3(i,3)/1.5
        ReacE7(i,2)= sqrt(2.)*(ME1(i,1,8)-ME3(i,3)/3.)*ratio**3
        ReacS7(i,1)= MS1(i,1,8)+MS3(i,3)/1.5
        ReacS7(i,2)= sqrt(2.)*(MS1(i,1,8)-MS3(i,3)/3.)*ratio**3
      enddo
            
c      if (ww.gt.wwlast) call write_mpoles(s,ratio)
      wwlast = ww
      
      RETURN
      END
