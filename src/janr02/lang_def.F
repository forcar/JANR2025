      SUBROUTINE LANG_DEF

      IMPLICIT NONE
      SAVE

#include "resonances.inc"
#include "fit_par.inc"

      INTEGER i,k

c- Langul(i) is the degree of photon momenta in Breit-Wigner formula for 
c- M(l+),M(l-): Langul(i)=l.
      Langul(1)=1
      Langul(2)=2
      Langul(3)=3
      Langul(4)=1  
      Langul(5)=0
      Langul(6)=3
      Langul(7)=2
      Langul(8)=1    
      Langul(9)=0
      Langul(10)=1
      Langul(11)=1
      Langul(12)=2
      Langul(13)=3

      Langul(21)=1
      Langul(22)=0
      Langul(23)=2  
      Langul(24)=2
      Langul(25)=0
      Langul(26)=3    
      Langul(27)=1
      Langul(28)=3
      Langul(29)=2
      Langul(30)=1
      Langul(31)=1
      Langul(32)=2  
      Langul(33)=3
	
c  Lprime(i) is the degree of photon momenta in Breit-Wigner formula 
c  for E(l+), E(l-):
c  Lprime(i)=l for E(l+), 
c  Lprime(i)=l-2 for E(l-) at l>1,
c  Lprime(i)=l for E(l-) at l=1.    	
      Lprime(1)=1    	
      Lprime(2)=0
      Lprime(3)=3
      Lprime(4)=1
      Lprime(5)=0
      Lprime(6)=1
      Lprime(7)=2
      Lprime(8)=1
      Lprime(9)=0
      Lprime(10)=1
      Lprime(11)=1
      Lprime(12)=0    	
      Lprime(13)=1
	
      Lprime(21)=1
      Lprime(22)=0
      Lprime(23)=0
      Lprime(24)=0
      Lprime(25)=0
      Lprime(26)=1
      Lprime(27)=1
      Lprime(28)=3
      Lprime(29)=2
      Lprime(30)=1
      Lprime(31)=1
      Lprime(32)=0    	
      Lprime(33)=1
      
c     eta are branching ratios to pion-nucleon channel
      eta(1)=1.
      eta(2)=0.16
      eta(3)=0.42
      eta(4)=0.1
      eta(5)=0.29
      eta(6)=0.1
      eta(7)=0.12
      eta(8)=0.01
      eta(9)=0.13
      eta(10)=0.16
      eta(11)=0.18
      eta(12)=0.01
      eta(13)=0.01


      eta(21)=0.6
      eta(22)=0.4
      eta(23)=0.54
      eta(23)=1.0
      eta(24)=0.1
      eta(25)=0.6
      eta(26)=0.62
      eta(27)=0.15
      eta(28)=0.06
      eta(29)=0.38
      eta(30)=0.18
      eta(31)=0.2
      eta(32)=0.1
      eta(33)=0.01

      RETURN
      END
