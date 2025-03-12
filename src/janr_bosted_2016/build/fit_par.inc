      INTEGER MAXpar
      PARAMETER (MAXpar = 99)
      CHARACTER*10 pname(MAXpar)
      REAL start_value(MAXpar),step_size(MAXpar)
      REAL low_lim(MAXpar),up_lim(MAXpar)
      INTEGER pi_type, pi_type2fit, obsv2fit, ndf
      INTEGER par_stat(MAXpar),plist(MAXpar),parms
      INTEGER lmini,lmino,lminf,f_out
      PARAMETER (f_out=50)
      DOUBLE PRECISION cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10				
      DOUBLE PRECISION cc11,cc12,cc13,cc14,cc15,cc16							
      DOUBLE PRECISION cc17,cc18,cc19,cc20								
      DOUBLE PRECISION cc21,cc22,cc23,cc24,cc25								
      DOUBLE PRECISION cc26,cc27,cc28,cc29,cc30								
      DOUBLE PRECISION cc31,cc32,cc33,cc34,cc35								
      DOUBLE PRECISION cc36,cc37,cc38,cc39,cc40								
      DOUBLE PRECISION cc41,cc42,cc43,cc44,cc45,cc46,cc47,cc48,cc49
      REAL ChiMax,fix_param,tolorance, maxcalls,errorup

      COMMON/FIT_PAR/cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10,
     +cc11,cc12,cc13,cc14,cc15,cc16,cc17,cc18,cc19,cc20,
     +cc21,cc22,cc23,cc24,cc25,cc26,cc27,cc28,cc29,cc30,
     +cc31,cc32,cc33,cc34,cc35,cc36,cc37,cc38,cc39,cc40,
     +cc41,cc42,cc43,cc44,cc45,cc46,cc47,cc48,cc49
	common/fit_parc/
     +ChiMax,fix_param,tolorance, maxcalls,errorup,pi_type2fit,
     +obsv2fit,pi_type
     +,start_value,step_size,low_lim,up_lim,par_stat,plist,
     +parms,lmini,lmino,lminf,ndf,pname




