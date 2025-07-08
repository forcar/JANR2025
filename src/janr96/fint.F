      Function Fint(x)
c
      IMPLICIT NONE
      SAVE
#include "rephi.inc"
      common/ifint/lm,is,i,j,w,q2
c
      Real Fint,x
      Real W,Q2
      Integer i,j,is,lm,ii
      real p,p1,p2,rim,PHis(6,3,2)
      integer ierr
c
      call janr_rephi(W,Q2,x,ierr)
c
      if(is.eq.3)then
         do ii=1,6
            PHis(ii,is,j)=Phisot(j,ii,3)
         enddo
      elseif(is.eq.2)then
         do ii=1,6
            PHis(ii,is,j)=Phisot(j,ii,1)+2.*Phisot(j,ii,2)
         enddo
      elseif(is.eq.1)then
         do ii=1,6
            PHis(ii,is,j)=Phisot(j,ii,1)-Phisot(j,ii,2)
         enddo
      endif
c
      goto(21,22,23,24)lm+1
 21   continue
      p=1.
      p1=0.
      p2=0.
      go to 111
c
 22   continue
      p=x
      p1=-sqrt(1.-x**2)
      p2=0.
      go to 111
c
 23   continue
      p=3.*x**2/2.-0.5
      p1=-3.*x*sqrt(1-x**2)
      p2=3.*(1.-x**2)
      go to 111
c
 24   continue
      p=(5.*x**2-3.)*x/2.
      p1=-sqrt(1.-x**2)*(15.*x**2-3.)/2.
      p2=15.*x*(1-x**2)
c
 111  continue
      fint=0.0
      goto(11,12,13,14,15,16)i
c
 11   continue
c - Integeral 1 (I1)
      if(lm.lt.1)return
      rim=p1*sqrt(1.-x**2)*(PHis(3,is,j)+x*PHis(4,is,j))
      fint=rim*(2.*lm+1.)/(2.*lm*(lm+1.))
      go to 113
c
 12   continue 
c - Integeral 2 (I2)
      rim=p*(2.*PHis(1,is,j)-2.*x*PHis(2,is,j)+(1.-x**2)*PHis(4,is,j))
      fint=rim*(2.*lm+1.)/2.
      go to 113
c
 13   continue
c - Integeral 3 (I3)
      if(lm.lt.1)return
      rim=p1*sqrt(1.-x**2)*PHis(2,is,j)
      fint=-rim*(2.*lm+1.)/(2.*lm*(lm+1.))
      go to 113
c
 14   continue
c - Integeral 4 (I4)
      if(lm.lt.2)return
      rim=p2*(1.-x**2)*PHis(4,is,j)
      fint=rim*(2.*lm+1.)/(2.*lm*(lm+1.)*(lm-1)*(lm+2))
      go to 113
c
 15   continue 
c - Integeral 5 (I5)
      rim=p*(PHis(5,is,j)+x*PHis(6,is,j))
      fint=rim*(2.*lm+1.)/2.
      go to 113
c
 16   continue
c - Integeral 6 (I6)
      if(lm.lt.1)return
      rim=p1*sqrt(1.-x**2)*PHis(6,is,j)
      fint=-rim*(2.*lm+1.)/(2.*lm*(lm+1.))
c
113   continue
      return
      end
