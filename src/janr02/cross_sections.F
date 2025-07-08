      SUBROUTINE CROSS_SECTIONS (s,azimut,epsilon)
      IMPLICIT NONE
      SAVE

#include "helicity.inc"
#include "cross_sec.inc"
#include "fit_par.inc"

      REAL s,azimut,qkratio,epsilon,sinth,eps,az
      REAL mpion(2),mnuc,qpi
      INTEGER NWRITE
      INTEGER i,j

      mnuc=0.93827
      mpion(1)=0.13498
      mpion(2)=0.13957
            
      DO i=1,2
      
        qpi = sqrt(((s-mnuc**2+mpion(i)**2)/2.)**2/s-mpion(i)**2)
        qkratio = qpi/((s-mnuc**2)/2./sqrt(s))   
        
        sigmaT(i)=0.
        sigmaL(i)=0.
        sigmaTT(i)=0.
        sigmaTL(i)=0.
        Do j=1,2
          sigmat(i) = sigmat(i) + 
     +                H1(i,j)**2.+ H2(i,j)**2.+ H3(i,j)**2.+ H4(i,j)**2.
          sigmaL(i) = sigmaL(i) + H5(i,j)**2. + H6(i,j)**2.
          sigmaTT(i) = sigmaTT(i)-H4(i,j)*H1(i,j)+H3(i,j)*H2(i,j)
          sigmaTL(i) = sigmaTL(i)+H5(i,j)*(H1(i,j)-H4(i,j))+ 
     -                 H6(i,j) * ( H2(i,j) + H3(i,j) )
        EndDo
      sigmaT(i)=sigmaT(i) * qkratio / 2.
      sigmaL(i)=sigmaL(i) * qkratio
      sigmaTT(i)=sigmaTT(i) * qkratio
      sigmaTL(i)=sigmaTL(i) * qkratio / sqrt(2.)
      
      sigma(i) = sigmaT(i) + epsilon * sigmaL(i) + 
     + epsilon * sigmaTT(i) * cos(2. * azimut) + sqrt( 2. * epsilon *
     *(1. + epsilon)) * sigmaTL(i) * cos( azimut )

      ENDDO
      RETURN
      END



