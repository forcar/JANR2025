      SUBROUTINE MIXING_AMPL (S,Q2)
     
      IMPLICIT NONE
      SAVE
    
#include "gpar.inc"
#include "fit_par.inc"
#include "back_ampl.inc"

      REAL S, Q2, kgamma, qpion
      INTEGER i, ll, j,k
      REAL W,tau,q1,F_dipole 
      REAL on1,on2,op1,om1,op2,om2 
      REAL z,Gmp,F1p,F2p,Gmn,Gen,F2n,B1(2),Lpv,Lmix 
      
c     -  here we calculate multipole amplitudes for PV-PS mixing 

      W   = sqrt(s)
      on1 = W-Egamma
      on2 = W-Epion
      op1 = on1+mn 
      om1 = on1-mn
      op2 = on2+mn
      om2 = on2-mn  
      
      tau   = Q2/4./mn22
      q1    = sqrt(Q2)
      F_dipole = 1./(1.+Q2/0.71)**2
      z     = 1.793/mn2/(1.+1.2*Q2/(1.+1.1*q1)+0.015*Q2+0.001*Q2*4)
      Gmp   = 1.+0.35*q1+2.44*q2+0.5*q1**3+1.04*q1**4+0.34*q1**5
      Gmp   = 2.793/Gmp
      F1p   = Gmp/(1.+mn2*z)
      F2p   = z*F1p
      Gen   = 0.5*Q2/(1.+25.*Q2**2)
      Gmn   = -1.913*F_dipole
      F2n   = (Gmn-Gen)/(1.+tau)/mn2
      Lmix  = 0.445*(1.03231466)
      Lpv   = Lmix**2/(Lmix**2+Epion**2-mp22)
      B1(1) = Lpv*3.212*(F2p-F2n)/mn/W/2.
      B1(2) = Lpv*3.212*(F2p+F2n)/mn/W/2.
      
C -this is loop for (+) and (0) amplitudes. (-)=0.

      do i=1,2
        E0(i) =  sqrt(op1*op2)*(W-mn)*B1(i)
        M1(i) = -sqrt(om1*om2)*(W+mn)*B1(i)
        S0(i) =  sqrt(om1*op2)*op1*B1(i)
        S1(i) = -sqrt(op1*om2)*om1*B1(i)
      enddo

      RETURN
      END
