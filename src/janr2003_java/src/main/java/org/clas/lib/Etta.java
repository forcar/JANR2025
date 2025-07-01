package org.clas.lib;

public class Etta extends Constants {
	
	public Etta() {
		
	}

    public void etta() {
        int[] ETampl = new int[6];
        int[] ETisot = new int[3];
        
        // Definition of coefficients that denote the symmetry of Invariant 
        // Amplitudes under the interchange "s<==>u"
        ETampl[0] = 1;  // Note: Java arrays are 0-based
        ETampl[1] = 1;
        ETampl[5] = 1;
        
        for (int i = 2; i < 5; i++) {
            ETampl[i] = -1;
        }

        // Coefficients which determine the symmetry of Isospin Amplitudes 
        // (+),(-),(0) under the interchange "s<==>u"
        ETisot[0] = 1;
        ETisot[1] = -1;
        ETisot[2] = 1;

        // Products of two symmetry coefficients
        for (int k = 0; k < 6; k++) {
            for (int m = 0; m < 3; m++) {
                ett[k][m] = ETampl[k] * ETisot[m];
            }
        }
    }

}
