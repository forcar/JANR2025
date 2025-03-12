        program convert
        real a1,a2,a3,a4,a5,a6,a7
        open(unit=2,file='error.dat',status='old')
           open (30,FILE='error_348.dat')
110     read(2,*,end=114)a1,a2,a3,a4,a5
	a7=sqrt(a2**2+a3**2+a4**2+a5**2)
           write(30,133) a2,a3,a4,a5,a7
        go to 110
114     continue
133    FORMAT(7f9.3)
        stop
        end


