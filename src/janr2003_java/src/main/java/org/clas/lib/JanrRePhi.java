package org.clas.lib;

import static java.lang.Math.*;

public class JanrRePhi extends Constants {
	
	Etta  e = new Etta();
	Born  b = new Born(); 
	RePhi r = new RePhi();
	
    double s, qkratio, kgamma, qpion, k22, q22;
   
    public JanrRePhi() {
    	
    }
    
    public void janrRePhi(double w, double q2, double costh, int ierr) {
    	
        s  = pow(w, 2);
        w2 = 2.0 * w;
        E  = s - mn22;
        
        Egamma = (E - q2) / w2;
        Epion  = (E + mp22) / w2;
        q22    = pow(Epion, 2) - mp22;
        k22    = pow(w - Egamma, 2) - mn22;
        qpion  = sqrt(q22);
        
        qkratio = qpion * w2 / E;
        kgamma  = sqrt(k22);
        
        qk = qpion * kgamma * costh - Epion * Egamma;
        t  = mp22 - q2 + 2.0 * qk;
        u  = 2.0 * mn22 + mp22 - q2 - s - t;
        
        e.etta();
        b.born(s,q2);
        r.rePhi(s,q2);
        
        ierr = 0;
    }

}
