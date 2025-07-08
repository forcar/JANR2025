c      include file FIT_PAR
c
c_begin_doc
c  RCS ID string
c  $Id: fit_par.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $
c
c  Documentation for include file FIT_PAR
c
c  Purpose:
c  --------
c
c  Used in routines:
c  ----------------
c
c  Notes:
c  ------
c
c  Author:         Created:  Wed Dec 11 17:25:11 EST 1996
c  -------
c
c  Major revisions:
c  ----------------
c
c_end_doc
c
c
c  Remember NO NESTED INCLUDES
c
c  COMMON BLOCK NAME should be the same as include file name!
c
c_begin_var
c  variables for include file FIT_PAR:
c  -----------------------------------
c  RCS information: 
      CHARACTER*(132)  CIFIT_PAR
      DATA CIFIT_PAR/
     1'$Id: fit_par.inc,v 1.1 2003/11/20 05:16:08 lcsmith Exp $'
     2/
c----6----------------------------------------------------------------72
c  
c  Variables contained in FIT_PAR common block
c  --------------------------------------------
c
c_end_var
      INTEGER MAXpar
      PARAMETER (MAXpar = 99)
      CHARACTER*10 pname(MAXpar)
      REAL start_value(MAXpar),step_size(MAXpar),low_lim(MAXpar),up_lim(MAXpar)
      INTEGER pi_type, pi_type2fit, obsv2fit, ndf
      INTEGER par_stat(MAXpar),plist(MAXpar),parms,lmini,lmino,lminf,f_out
      PARAMETER (f_out=49)
      DOUBLE PRECISION cM1,cE1,cS1								
      DOUBLE PRECISION cM2,cE2,cS2								
      DOUBLE PRECISION cM3,cE3,cS3								
      DOUBLE PRECISION cMP1,cSP1,cEP2,cSP2,cEP5,cSP5
      DOUBLE PRECISION cMP3,cEP3,cSP3
      DOUBLE PRECISION cMP4,cEP4,cSP4
      DOUBLE PRECISION cMP6,cEP6,cSP6
      DOUBLE PRECISION cMP7,cEP7,cSP7
      DOUBLE PRECISION cMP12,cEP12,cSP12,a_pip
      REAL ChiMax,fix_param,tolorance, maxcalls,errorup
c
      COMMON/FIT_PAR/cM1,cE1,cS1,cM2,cE2,cS2,cM3,cE3,cS3,
     +cMP1,cSP1,cEP2,cSP2,cEP5,cSP5,cMP3,cEP3,cSP3,
     +cMP4,cEP4,cSP4,cMP6,cEP6,cSP6,cMP7,cEP7,cSP7,cMP12,
     +cEP12,cSP12,a_pip
	common/fit_parc/
     +ChiMax,fix_param,tolorance, maxcalls,errorup,pi_type2fit,obsv2fit,
     +pi_type
     +,start_value,step_size,low_lim,up_lim,par_stat,plist,
     +parms,lmini,lmino,lminf,ndf,pname
c----6----------------------------------------------------------------72




