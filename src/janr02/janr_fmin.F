      SUBROUTINE JANR_FMIN(npar,grad,fval,xval,iflag,futil)
      IMPLICIT NONE
      SAVE

#include "cross_sec.inc"
#include "res_func.inc"
#include "fit_par.inc"
#include "janr_data.inc"
#include "janrtab.inc"

      DOUBLE PRECISION grad,xval,fval
      dimension grad(35),xval(35)
      REAL  futil
      INTEGER iflag, npar
      INTEGER i,itr_number,Npoints_max,i_point,ierr
      integer itr_count,nhic,itr_bad,ip,N_zero,N_bad,jitr
      integer  obsv
      REAL diff,cross_calc,fitted_points
      REAL W,costh,Q2,epsilon,phi,sec,errors,fval_last
      REAL timdd, timdd1, timdd2, time_add, ac, rad,nb
      integer fitop
      logical bad_xsec,bad_rad,bad_acc,bad_stat
      real amod
      real var_calc, var_meas , sigma_total
      real c2weight
      integer arg0, arg1, argt
      double precision norm_scale(10) 
      data norm_scale / 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. /
      
      IF ( iflag .EQ. 1 ) THEN
        OPEN (UNIT=f_out,FILE=inffile,status='UNKNOWN')
        CALL TIMEX(timdd2)
        time_add   = 0.
        fitop      = 1
        itr_number = 0
        itr_count  = 0
        Npoints_max= 0
        N_zero     = 0
        N_bad      = 0
        nhic       = 0
        
        IF( Datafile.NE.'JUNK' .OR. Datafile.NE.'junk' ) THEN
          OPEN (UNIT=44,FILE=Datafile,status='OLD')
          do i=1,MAXpoints
          
            if(datform.eq.'o'.or.datform.eq.'O')then
            
              READ(44,*,err=999,end=998)W,Q2,epsilon,costh,phi,sec,errors
              obsv        = 1
              pi_type     = pi_type2fit 
              Npoints_max = Npoints_max+1
              bad_xsec    = sec.le.0..and.errors.le.0
              if(bad_xsec)then
                point_stat(Npoints_max) = 1
                N_zero = N_zero+1
              else
                point_stat(Npoints_max) = 0
              endif
              
            elseif(datform.eq.'n'.or.datform.eq.'N')then
            
              READ(44,*,err=999,end=998)W,Q2,epsilon,costh,phi,sec,errors,
     +        nb,ac,rad,pi_type,obsv
              Npoints_max  = Npoints_max+1
              bad_xsec     = sec.le.0..and.errors.le.0
              bad_stat     = nb.lt.5
              bad_acc      = ac.lt.0.05.or.(ac.lt.0.1.and.nb.lt.10)
              bad_rad      = rad.lt.1.
              if(bad_xsec)then
                point_stat(Npoints_max) = 1
                N_zero = N_zero+1
              else
                point_stat(Npoints_max) = 0
              endif
              arg0 = 10**(obsv) 
              arg1 = 10**(obsv-1)
              pi_type2fit =  amod( float(obsv2fit), float(arg0) ) / arg1
              
            endif

            if(w.lt.wmin.or.w.gt.wmax)             point_stat(Npoints_max)=1
            if(q2.lt.q2min.or.q2.gt.q2max)         point_stat(Npoints_max)=1
            if(costh.lt.cosmin.or.costh.gt.cosmax) point_stat(Npoints_max)=1
            if(phi.lt.phimin.or.phi.gt.phimax)     point_stat(Npoints_max)=1

            if (((pi_type.ne.pi_type2fit).AND.(pi_type2fit.ne.9)).OR. 
     +           (pi_type2fit.eq.0))               point_stat(Npoints_max)=1 

            Wd(Npoints_max)       = W
            Q2d(Npoints_max)      = Q2
            epsilond(Npoints_max) = epsilon
            costhd(Npoints_max)   = costh
            phi_str(Npoints_max)  = phi
            cross_sec(Npoints_max)= sec
            st_error(Npoints_max) = errors
            pi_typed(Npoints_max) = pi_type
            obsvd(Npoints_max)    = obsv
            
          enddo   
           
 998      close(44)
          
          OPEN (UNIT=53,FILE='start_tab',status='UNKNOWN')
          Write(lmino,*)'Number of data points =  ',Npoints_max
          Write(lmino,*)'Number of NON zero data points =  ',Npoints_max-N_zero
          write(f_out,*)'Kinematics:   W  =',wmin,' - ',wmax
          write(f_out,*)'              Q2 =',q2min,' - ',q2max
          Write(f_out,*)'Number of data points =  ',Npoints_max
          Write(f_out,*)'Number of NON zero data points =  ',Npoints_max-N_zero
          
c     Write(6,*)
c     +'Order is -> W, Q2, epsilon, cos(*), phi*, sigma, err_st'
c     do i=1, Npoints_max       
c     write(6,*)Wd(i),Q2d(i),epsilond(i),costhd(i),
c     +azim_gradd(i),cross_sec(i),st_error(i)
c     enddo

          DO I_point=1,Npoints_max
            do i=1,Nbad
              if(i_point.gt.Ibad(1,i).and.I_point.lt.Ibad(2,i)) then
                point_stat(I_point) = 1
                N_bad = N_bad+1
              endif
            enddo  
          ENDDO
          Write(lmino,*)'Number of data points for first iteration =  ',Npoints_max-N_zero-N_bad
        ENDIF
        
        write(6,111)(xval(i),i=1,parms)
111     format(20x,'Value of parameters at the begining => 1-30 ',/,
     +'   cM1 =',f9.5,'   cE1 =',f9.5,'   cS1 =',f9.5,/,
     +'   cM2 =',f9.5,'   cE2 =',f9.5,'   cS2 =',f9.5,/,
     +'   cM3 =',f9.5,'   cE3 =',f9.5,'   cS3 =',f9.5,/,
     +'   cMP1 =',f9.5,'   cSP1 =',f9.5,'   cEP2 =',f9.5,'   cSP2 =',f9.5,/,
     +'   cMP3 =',f9.5,'   cEP3 =',f9.5,'   cSP3 =',f9.5,/,
     +'   cMP4 =',f9.5,'   cEP4 =',f9.5,'   cSP4 =',f9.5,/,
     +'   cEP5 =',f9.5,'   cSP5 =',f9.5,/,
     +'   cMP6 =',f9.5,'   cEP6 =',f9.5,'   cSP6 =',f9.5,/,
     +'   cMP7 =',f9.5,'   cEP7 =',f9.5,'   cSP7 =',f9.5,/,
     +'   cMP12 =',f9.5,'   cEP12 =',f9.5,'   cSP12 =',f9.5,/)
     
      ENDIF
      
      IF ( iflag .EQ. 2 ) THEN
         print *,'Changing option'
         fitop=0
      ELSEIF ( iflag .EQ. 3 ) THEN
         Write(f_out,*)'Number of data points to fit =  ',ndf
         write(f_out,*)'Number of iterations = ',itr_number,'  Chi2/point',fval_last/float(ndf)
         call prterr
         itr_number=0
         OPEN (UNIT=53,FILE=tabfile,status='UNKNOWN')
      ENDIF

      IF(iflag.eq.1.or.iflag.ge.3)THEN
        timdd1   = timdd2
        CALL TIMEX(timdd2)
        timdd    = timdd2-timdd1
        time_add = time_add+timdd
        write(6,*)
        write (*,*)'time for iteration = ',timdd,'  Total time = ',time_add
        fitted_points=0
        fval=0.
        cM1=sngl(xval(1))
        cE1=sngl(xval(2))
        cS1=sngl(xval(3))
        cM2=sngl(xval(4))
        cE2=sngl(xval(5))
        cS2=sngl(xval(6))
        cM3=sngl(xval(7))
        cE3=sngl(xval(8))
        cS3=sngl(xval(9))
        cMP1=sngl(xval(10))
        cSP1=sngl(xval(11))
        cEP2=sngl(xval(12))
        cSP2=sngl(xval(13))
        cMP3=sngl(xval(14))
        cEP3=sngl(xval(15))
        cSP3=sngl(xval(16))
        cMP4=sngl(xval(17))
        cEP4=sngl(xval(18))
        cSP4=sngl(xval(19))
        cEP5=sngl(xval(20))
        cSP5=sngl(xval(21))
        cMP6=sngl(xval(22))
        cEP6=sngl(xval(23))
        cSP6=sngl(xval(24))
        cMP7=sngl(xval(25))
        cEP7=sngl(xval(26))
        cSP7=sngl(xval(27))
        cMP12=sngl(xval(28))
        cEP12=sngl(xval(29))
        cSP12=sngl(xval(30))
        norm_scale(2)=(xval(31))

        DO I_point=1,Npoints_max
          if(point_stat(I_point).eq.0)then
            W       = Wd(i_point)
            Q2      = Q2d(i_point)
            epsilon = epsilond(i_point)
            costh   = costhd(i_point)
            phi     = phi_str(i_point)
            pi_type = pi_typed(i_point)
            obsv    = obsvd(i_point)
               
            CALL JANR_RUN(W,Q2,epsilon,costh,phi,ierr)
               
            IF(ierr.eq.0) THEN
               
              if (obsv.eq.1) then	! Cross section
                var_calc = sigma(pi_type)
                var_meas = norm_scale(pi_type) * cross_sec(I_point)
                if(pi_type.eq.2) then 
                  c2weight = 1.
                else
                  c2weight = 1.
                endif
                   
                diff=c2weight*((var_calc-var_meas)/st_error(i_point))**2
                  
              else if (obsv.eq.2) then	!Sig_ltp
                sigma_total  = sigma(pi_type)
                var_calc = sigmaTLP(pi_type)*0.67*sqrt(2.*epsilon*(1.-epsilon))
     +                    *sin(phi*acos(0.)/90.)/sigma_total 
                var_meas = cross_sec(I_point)
                c2weight = 1.
c                 c2weight = 2.

                diff=c2weight*((var_calc-var_meas )/st_error(i_point))**2
     
              endif 

c                  cross_calc = sigma(pi_type)
c                  diff=((cross_calc-cross_sec(I_point))/st_error(i_point))**2

              if(diff.gt.chimax.and.itr_number.eq.0.and.fitop.eq.0)then
                nhic = nhic+1
                point_stat(I_point) = 1
                write(55,134)nhic,I_point,W,Q2,costh,phi,cross_sec(I_point),var_calc,diff
                itr_bad = nhic
              else
                if(itr_number.eq.0.and.(iflag.eq.1.or.iflag.eq.3))then
                write(53,133)I_point,W,Q2,epsilon,costh,phi,cross_sec(I_point),
     +st_error(i_point),var_calc,diff
                endif
                
                fitted_points = fitted_points + 1.
                fval          = fval + diff
                fval_last     = fval
              endif
            Else
              write(6,*)'Bad calculation',w,q2,costh,phi,i_point
            ENDIF
          endif
        ENDDO
        
        print *,Npoints_max,fitted_points,nhic
        if(itr_number.eq.0)then
          Write(lmino,*)'The number of data points for fit =  ',Npoints_max-N_zero-N_bad-nhic
          if(iflag.eq.1)close(53)
          ndf=fitted_points
        endif
         
        if(iflag.ne.1)itr_number=itr_number+1
        if(iflag.ne.1)itr_count=itr_count+1
        write(6,*)'Iteration = ',itr_number,'    Value to minimize',fval,'  Chi2/point',fval/float(ndf)
        write(6,*)'Number fit parameters',int(parms-fix_param),' Number of parameters',parms,fix_param,
     +' Number of points for fit',int(fitted_points)
        write(26,*)'Iteration = ',itr_number,'    Value to minimize',fval,'  Chi2/point',fval/float(ndf)
        do ip=1,parms
          if((par_stat(ip).le.0.and.fitop.eq.0).or.(abs(par_stat(ip)).eq.1.and.fitop.eq.1))then
            write(6,*)pname(ip),'  ',sngl(xval(ip))
            write(26,*)pname(ip),'  ',xval(ip)
          endif
        enddo
      
      ENDIF
 133  format(i7,8f9.4,f10.5)
 134  format(2i7,7f8.3)
      RETURN
      
 999  write(6,*)'Error while reading the data file'
 
      RETURN
      END


