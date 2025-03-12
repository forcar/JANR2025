# JANR2025
## Contents
**JANE.tar**: JANR release from Inna Aznauryan dated Sept 25 2014.  Used by Nickolay Markov and possibly Kijun Park.  

Major code changes include: 
- replace **getbreit.F** with **breit_ampl.F**.
- changes to eta channel for S11(1535)
- extensive changes in multipole amplitude Q2 parameterization and expansion of number of fitting coefficients to accomodate background fits.
- removal of Regge amplitudes which are deemed by Inna as unnecessary for Q2 > 0.4 GeV.  However this means code is not usable for lower Q2.
    
**src/janr_bosted_2016**: 

Peter Bosted used Inna's 2014 release to generate predictions of beam-target and target asymmetry
observables for CLAS eg1-dvcs, eg1b, eg1c, eg4, e16 experiments.  The Bosted version of Inna's code removed MINUIT fitting subroutines and
condensed the remaining separate subroutine files into one package.  Otherwise this version is identical to Inna's 2014 release.  

Note that this code was initially used incorrectly by Bosted as I summarized here:  https://clasweb.jlab.org/wiki/index.php/Discussion_of_Zheng/Bosted_Papers

**src/janr03**: 

Archive of my attempt in 2003 to clean up JANR I/O and code organization and using the CERNLIB PAW package for interactive fitting visualization.  The description is given here: https://userweb.jlab.org/~lcsmith/galileo4/codes/janr2003/

**NOTES:**
