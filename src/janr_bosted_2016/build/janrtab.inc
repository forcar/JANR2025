      INTEGER maxevents, maxiter, randomseed
      REAL wmin,wmax,wstep
      REAL q2min,q2max,q2step
      REAL cosmin,cosmax,cosstep
      REAL phimin,phimax,phistep
      REAL ebeam, intAccur
      INTEGER maxq2grid, maxwgrid, maxcosgrid, maxphigrid
      PARAMETER (maxq2grid = 20)
      PARAMETER (maxwgrid = 100)
      PARAMETER (maxcosgrid = 20)
      PARAMETER (maxphigrid = 30)
      REAL cs_tab(maxq2grid,maxwgrid,maxcosgrid,maxphigrid)
      CHARACTER*60 tabfile,inffile

      COMMON/janrtab/ebeam,wmin,wmax,wstep,q2min,q2max,q2step,intaccur,
     +cosmin,cosmax,cosstep,phimin,phimax,phistep,cs_tab,tabfile,inffile




