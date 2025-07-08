      SUBROUTINE BREIT_AMPL (S,Q2,kgamma, qpion)
c     
c     _begin_doc
c     RCS ID string
c     $Id$
c     
c     Documentation for subroutine BREIT_AMPL
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
c     Author:         Created:  Wed Dec 11 17:43:38 EST 1996
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
#include "gpar.inc"
#include "resonances.inc"
#include "fit_par.inc"
#include "amplitudes.inc"
#include "delta_dat.inc"
c     _end_inc
c     
c     _begin_var
c     input/output variables:
c     -----------------------
      REAL S, Q2, kgamma, qpion
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
      CHARACTER*100  CRMESS
      PARAMETER (CRNAME='BREIT_AMPL')
      PARAMETER (CRAUTH='')
c     
c     Local User defined variables:
c     -----------------------------
c     INTEGER NWRITE
      INTEGER i, ll, j,k
      REAL M22,M2,En1,En2,kr22,qr22,kr,qr, X2, q22,k22,den
      REAL E2pion,qr2pion,r2pion,qxl,qx2l,Gpion,Ginel,Ggamma
      REAL rvirtual,rpion,WW,Wthr,Eeta,qreta,reta,Breit(2)
      REAL BreitW,Gtotal,Gtotal1,REden,IMden,qxgamma 
      REAL kreal,kRreal,W2pion
      REAL XIN,YOUTOM,YOUTOE,YOUTOS
c     _end_var
c     
c     executable code for routine BREIT_AMPL:
c---- 6----------------------------------------------------------------72
c     -  Calculate Real and Imaginary parts of Breit-Wigner electroproduction amplitudes
c     -  Breit(1),Breit(2) 
c     -  Real and Imaginary parts of Breit-Wigner electroproduction multipoles
c     -  BreitM3(1,j),BreitM3(2,j) etc for I=3/2 
c     -  BreitM1(1,i,j),BreitM1(2,i,j) etc for I=1/2, 
c     -  i=1 for pM(1/2), i=2 for nM(1/2) 
        WW=sqrt(s)
      q22=qpion**2
      k22=kgamma**2
	W2pion=mn+2.*mp
      if (WW.le.W2pion) go to 166
      go to 167
166   q2pion=0.
167   continue
c     -  First do the loop for the P33(1232) resonance      
	DO i=1,1
         M22=MRES(i)**2
         M2=2.*MRES(i)
         En1=(M22+mn22+Q2)/M2  
         En2=(M22+mn22-mp22)/M2
         kr22=En1**2-mn22
         qr22=En2**2-mn22
         kr=sqrt(kr22)
         qr=sqrt(qr22)
         kreal=(s-mn22)/2./WW  
         kRreal=(M22-mn22)/M2
	ratio=kgamma*kRreal/(kr*kreal)  
         X2=XRES(i)**2
         rvirtual=kgamma/kr
         rpion=qpion/qr
         E2pion=(M22-mn22+m2pion**2)/M2
         qr2pion=E2pion**2-m2pion**2
         r2pion=(q2pion/qr2pion)**(Langul(i)+2)
         qxl=((X2+qr22)/(q22+X2))**Langul(i)
         qx2l=((X2+qr2pion)/(q2pion+X2))**(Langul(i)+2)   
         Gpion=GRES(i)*rpion**(2.*Langul(i)+1.)*qxl
         Gtotal=Gpion/sqrt(s/M22)
         qxgamma=((X2+kr22)/(k22+X2))**Lprime(i)
         Ggamma=GRES(i)*rvirtual**(2.*Lprime(i)+1.)*qxgamma
         den=(M22-s)**2+(MRES(i)*Gtotal)**2
         BreitW=MRES(i)*sqrt(Gpion*Ggamma/(rvirtual*rpion))/den
         Breit(1)=BreitW*(M22-s)
         Breit(2)=BreitW*MRES(i)*Gtotal
c     -  Do the loop for real (j=1) and imaginary (j=2) parts
         do j=1,2
            BreitE3(j,i)=Breit(j)*AE3(i)
            BreitM3(j,i)=Breit(j)*AM3(i)
            BreitS3(j,i)=Breit(j)*rvirtual*AS3(i)
         enddo
      ENDDO
c     -  Do the loop for other resonances with isospin 3/2      
      DO i=2,8
         M22=MRES(i)**2	
         M2=2.*MRES(i)
         En1=(M22+mn22+Q2)/M2
         En2=(M22+mn22-mp22)/M2
         kr22=En1**2-mn22
         qr22=En2**2-mn22
         kr=sqrt(kr22)
         qr=sqrt(qr22)
         X2=XRES(i)**2
         rvirtual=kgamma/kr
         rpion=qpion/qr
         E2pion=(M22-mn22+m2pion**2)/M2
         qr2pion=E2pion**2-m2pion**2
         r2pion=(q2pion/qr2pion)**(Langul(i)+2)
         qxl=((X2+qr22)/(q22+X2))**Langul(i)
         qx2l=((X2+qr2pion)/(q2pion+X2))**(Langul(i)+2)
         Gpion=eta(i)*GRES(i)*rpion**(2.*Langul(i)+1.)*qxl         
         Ginel=(1.-eta(i))*GRES(i)*r2pion*qx2l
         Gtotal=Gpion+Ginel 
         qxgamma=((X2+kr22)/(k22+X2))**Lprime(i)
         Ggamma=GRES(i)*rvirtual**(2.*Lprime(i)+1.)*qxgamma
         den=(M22-s)**2+(MRES(i)*Gtotal)**2
         BreitW=MRES(i)*sqrt(Gpion*Ggamma/(eta(i)*rvirtual*rpion))/den
c	 BreitW=BreitW*MRES(i)/WW
         Breit(1)=BreitW*(M22-s)
         Breit(2)=BreitW*MRES(i)*Gtotal
c     -  Do the loop for real (j=1) and imaginary (j=2) parts		 
         do j=1,2
            BreitE3(j,i)=Breit(j)*AE3(i)
            BreitM3(j,i)=Breit(j)*AM3(i)
            BreitS3(j,i)=Breit(j)*rvirtual*AS3(i)
         enddo
      ENDDO
c     
      DO k=9,13
         DO i=1,2
            BreitE3(i,k)=0.
            BreitM3(i,k)=0.
            BreitS3(i,k)=0.
         ENDDO	
      ENDDO	
c     -  Now calculate the same for resonances with isospin 1/2.	
      DO ll=1,7
         i=ll+20
         M22=MRES(i)**2	
         M2=2.*MRES(i)
         En1=(M22+mn22+Q2)/M2
         En2=(M22+mn22-mp22)/M2
         kr22=En1**2-mn22
         qr22=En2**2-mn22
         kr=sqrt(kr22)
         qr=sqrt(qr22)
         X2=XRES(i)**2
         rvirtual=kgamma/kr
         rpion=qpion/qr
         E2pion=(M22-mn22+m2pion**2)/M2
         qr2pion=E2pion**2-m2pion**2
         r2pion=(q2pion/qr2pion)**(Langul(i)+2)
         qxl=((X2+qr22)/(q22+X2))**Langul(i)
         qx2l=((X2+qr2pion)/(q2pion+X2))**(Langul(i)+2)
         Gpion=eta(i)*GRES(i)*rpion**(2.*Langul(i)+1.)*qxl         
         Ginel=(1.-eta(i))*GRES(i)*r2pion*qx2l
         Gtotal=Gpion+Ginel
         qxgamma=((X2+kr22)/(k22+X2))**Lprime(i)
         Ggamma=GRES(i)*rvirtual**(2.*Lprime(i)+1.)*qxgamma
         den=(M22-s)**2+(MRES(i)*Gtotal)**2
         BreitW=MRES(i)*sqrt(Gpion*Ggamma/(eta(i)*rvirtual*rpion))/den
c         BreitW=BreitW*MRES(i)/WW
	 Breit(1)=BreitW*(M22-s)
         Breit(2)=BreitW*MRES(i)*Gtotal
c     -  Do the loop for real (j=1) and imaginary (j=2) parts	
c     -  Do the loop for  pA(1/2) (k=1) and nA(1/2)(k=2) amplitudes
         do j=1,2
            do k=1,2
               BreitE1(j,k,ll)=Breit(j)*AE1(k,ll)
               BreitM1(j,k,ll)=Breit(j)*AM1(k,ll)
               BreitS1(j,k,ll)=Breit(j)*rvirtual*AS1(k,ll)
            enddo
         enddo
      ENDDO
c     
      DO ll=8,13
         do k=1,2
            do j=1,2
               BreitE1(j,k,ll)=0.
               BreitM1(j,k,ll)=0.
               BreitS1(j,k,ll)=0.
            enddo
         enddo
      ENDDO
c     -  Now calculate the same for S11(1535) resonance		
      Wthr=mn+meta
      M22=MRES(22)**2
      M2=2.*MRES(22)
      En1=(M22+mn22+Q2)/M2
      En2=(M22+mn22-mp22)/M2
      kr22=En1**2-mn22
      qr22=En2**2-mn22
      kr=sqrt(kr22)
      qr=sqrt(qr22)
      X2=XRES(22)**2
      rvirtual=kgamma/kr 
      rpion=qpion/qr
      E2pion=(M22-mn22+m2pion**2)/M2
      qr2pion=E2pion**2-m2pion**2
      Eeta=(M22-mn22+meta**2)/M2
      qreta=Eeta**2-meta**2
      reta=abs(qeta/qreta)
      reta=sqrt(reta)
      r2pion=(q2pion/qr2pion)**2
      qx2l=((X2+qr2pion)/(q2pion+X2))**2
      Gtotal1=GRES(22)*rpion
      BreitW=MRES(22)*GRES(22)
      Gpion=0.7*GRES(22)*rpion
      Gtotal=Gpion
c      REden=M22-s+0.3*GRES(22)*MRES(22)*reta
      REden=M22-s
      IMden=MRES(22)*Gtotal
      den=REden**2+IMden**2
      Breit(1)=BreitW*REden/den
      Breit(2)=BreitW*IMden/den
      do j=1,2
         do k=1,2
            BreitE1(j,k,2)=Breit(j)*AE1(k,2)
            BreitS1(j,k,2)=Breit(j)*rvirtual*AS1(k,2)
         enddo
      enddo
      if (WW.ge.Wthr) go to 566
      go to 567
566   REden=M22-s
      IMden=MRES(22)*(Gtotal+0.3*GRES(22)*reta)
      den=REden**2+IMden**2
      Breit(1)=BreitW*REden/den
      Breit(2)=BreitW*IMden/den
      do j=1,2
         do k=1,2
            BreitE1(j,k,2)=Breit(j)*AE1(k,2)
            BreitS1(j,k,2)=Breit(j)*rvirtual*AS1(k,2)
         enddo
      enddo 
567   continue
      RETURN
      END
c
c------------------------------------------------------------------------------




