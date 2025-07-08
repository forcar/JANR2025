      subroutine janr_spp(q2,w,eps,csm,phi,chan,crs,sltp)
      
      real q2,w,eps,csm,phi,crs,sltp
      integer chan,ierr

#include "cross_sec.inc"
#include "fit_par.inc"
#include "res_func.inc"

      pi_type = chan
      
      cm1 = 4.7276

c Fit 1
c      cm1=4.81
c      ce1=-0.138
c      cs1=-0.252
c      cm2=5.01
      
c Fit 2      
c      cm1=4.839
c      ce1=-0.092
c      cs1=-0.243
c      cm2=5.85
      
      ierr = 0
      call janr_run(w,q2,eps,csm,phi,ierr)
      crs = sigma(chan)
      sltp = sigmaTLP(chan)
      end
      
      