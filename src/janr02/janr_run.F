      SUBROUTINE JANR_RUN(Wi,Q2i,epsilon,costhi,phidi,ierr)

      IMPLICIT NONE
      SAVE

#include "inpar.inc"
#include "gpar.inc"
#include "cross_sec.inc"

      REAL Wi, Q2i, epsilon, costhi, phidi
      INTEGER ierr

      INTEGER i,j
      REAL sinth,theta,im_phi_save(6,3),qkratio,phi
      REAL kgamma,qpion,k22,q22,E2pion,Eeta
      REAL W, costh, phid
      REAL Q2, S

      w		= wi
      q2	= q2i
      costh	= costhi
      phid	= phidi

      call janr_ini_point(Q2)
      
      theta	= acos(costh)
      theta	= theta*180./pi
      sinth	= sqrt(1.-costh**2)
      phi	= phidi*pi/180.

      s		= w**2
      w2	= 2.*w
      E	 	= s-mn22
      Egamma	= (E-Q2)/w2
      Epion	= (E+mp22)/w2
      q22	= Epion**2-mp22
      k22	= (w-Egamma)**2-mn22
      qpion	= sqrt(q22)
      qkratio	= qpion*w2/E
      kgamma	= sqrt(k22)
      qk	= qpion*kgamma*costh-Epion*Egamma
      t		= mp22-Q2+2.*qk
      u		= 2.*mn22+mp22-Q2-s-t
      E2pion	= (E+m2pion**2)/w2
      q2pion	= E2pion**2-m2pion**2
      Eeta	= (E+meta**2)/w2
      qeta	= Eeta**2-meta**2
      
      call phi_ampl(s,costh,Q2,kgamma,qpion)

      call hel_ampl(s,Q2,theta,kgamma)

      call response_functions(s,phi,epsilon)

      call cross_sections(s,phi,epsilon)
      
      ierr = 0
      
      RETURN
      END









