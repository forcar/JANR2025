      INTEGER MAXpoints, Ibad(2,100), Nbad
      PARAMETER ( MAXpoints = 100000 )
      CHARACTER*60 Datafile
      CHARACTER*12 datform
      INTEGER point_stat(MAXpoints)
      REAL Wd(MAXpoints), costhd(MAXpoints), Q2d(MAXpoints)
      REAL phi_str(MAXpoints),epsilond(MAXpoints)
      REAL cross_sec(MAXpoints),st_error(MAXpoints)
      integer pi_typed(MAXpoints), obsvd(MAXpoints)
      Common/janr_data/Wd,costhd,Q2d,phi_str,epsilond,cross_sec,
     +st_error,Datafile,Ibad,Nbad,point_stat,datform, pi_typed, obsvd
 
