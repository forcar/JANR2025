package org.clas.lib;

import static java.lang.Math.*;

public class Born extends Constants{
	
	GAMMF g = new GAMMF();

    public Born() {
    	
    }
    
    public void born(double s, double q2) {
    	
    	// Calculate Born term contribution (nucleon and pion poles) into invariant amplitudes.
    	// zero-indexing arrays only
    	
        double d = 1.0 / (s - mn22);
        double d1 = u - mn22;
        double tau = q2 / 4.0 / mn22;
        double q1 = sqrt(q2);
        
        // Pion Formfactor
        double f_pion = 1.0 / (1.0 + q2 / 0.54);
        
        // Dipole formula
        double f_dipole = 1.0 / pow(1.0 + q2 / 0.71, 2);
        
        // z=F2p/F1p from the data on Gep/Gmp: PRL 84(2000)1398
        double z = 1.793 / mn2 / (1.0 + 1.2 * q2 / (1.0 + 1.1 * q1) + 0.015 * q2 + 0.001 * q2 * 4);
        
        // Bosted parameterization of Gmp Phys. Rev. C 51, 409 (1995)
        double gmp = 1.0 + 0.35 * q1 + 2.44 * q2 + 0.5 * pow(q1, 3) + 1.04 * pow(q1, 4) + 0.34 * pow(q1, 5);
        gmp = 2.793 / gmp;
        
        // Proton Pauli Formfactors
        double f1p = gmp / (1.0 + mn2 * z);
        double f2p = z * f1p;
        
        // Gen and Gmn from Kees de Jager JLAB-PHY-00-01
        double gen = 1.6 * 0.5 * q2 / (1.0 + 25.0 * q2 * q2);
        double gmn = -1.913 * f_dipole;
        
        // Neutron Pauli Formfactors
        double f1n = (gen + tau * gmn) / (1.0 + tau);
        double f2n = (gmn - gen) / (1.0 + tau) / mn2;
        
        // Isovector and Isoscalar Pauli Formfactors (1 is V, 2 is S)
        double[] f_1 = new double[2];
        double[] f_2 = new double[2];
        f_1[0] = f1p - f1n;
        f_1[1] = f1p + f1n;
        f_2[0] = f2p - f2n;
        f_2[1] = f2p + f2n;
        
        // This loop counts for isotopic amplitudes (+) and (-),
        // which are determined by isovector form factors
        double[] ef1 = new double[3];
        double[] ef2 = new double[3];
        for (int k = 0; k < 2; k++) {
            ef1[k] = f_1[0];
            ef2[k] = f_2[0];
        }
        // The following are isotopic form factors (0) 
        // which are equal to isoscalar form factors
        ef1[2] = f_1[1];
        ef2[2] = f_2[1];
        
        // This loop counts for isotopic states (+),(-),(0)	
        double[][] re = new double[6][3];
        double[][] rm = new double[6][3];
        
        for (int k = 0; k < 3; k++) {
            re[0][k] =  ef1[k];
            re[1][k] = -ef1[k];
            re[2][k] = -ef1[k] / 2.0;
            re[3][k] = 0.0;
            re[4][k] = 0.0;
            re[5][k] = 0.0;
            
            rm[0][k] = mn2 * ef2[k];
            rm[1][k] = 0.0;
            rm[2][k] = 0.0;
            rm[3][k] = ef2[k];
            rm[4][k] = -qk * ef2[k];
            rm[5][k] = 2.0 * ef2[k];
        }
        
        // The coefficient 3.134 (microbarn)**(1/2) comes from
        // g*e/(8piGeV), where g**2/4pi=13.8, e**2/4pi=1/137.
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 3; j++) {
                re[i][j] = 3.134 * re[i][j] * (d + ett[i][j] / d1);
                rm[i][j] = 3.134 * rm[i][j] * (d + ett[i][j] / d1);
            }
        }
        
        // Pion pole term (previously used f_pion = f_1[0])        
        re[2][1] = re[2][1] - 6.268 * f_pion / (t - mp22);
                
        //Rho, Omega exchange
        double   gvec = f_dipole;
        double lamdom = 1.2;
        double lamdro = 1.5;
        
        double gvecom = gvec * pow(lamdom, 2) / (pow(lamdom, 2) - t) / (t - 0.612);
        double gvecro = gvec * pow(lamdro, 2) / (pow(lamdro, 2) - t) / (t - 0.591);
        
        double gvec1 = gvec / (t - 0.591);
        double   gro =  0.75 * gvecro;
        double  gro2 =   2.2 * gvecro;
        double   gom =  22.6 * gvecom;
        double  gom2 = -6.44 * gvecom;
        double  gro1 =   1.3 * gvec1;
        double gro21 =   2.2 * gvec1;
        double  gom1 =  15.0 * gvec1;
        double gom21 = 0.0;
        
        // Regge amplitudes
        double[][]  rro = new double[6][3];
        double[][]  rom = new double[6][3];
        double[][] rro1 = new double[6][3];
        double[][] rom1 = new double[6][3];
        
        rro[0][2] = mn2 * gro + t * gro2;
        rro[1][2] = 0.5 * (q2 + mp22 - t) * gro2;
        rro[2][2] = (u - s) * gro2 / 4.0;
        rro[3][2] = 0.0;
        rro[4][2] = 0.0;
        rro[5][2] = 2.0 * gro;
        
        rom[0][0] = mn2 * gom + t * gom2;
        rom[1][0] = 0.5 * (q2 + mp22 - t) * gom2;
        rom[2][0] = (u - s) * gom2 / 4.0;
        rom[3][0] = 0.0;
        rom[4][0] = 0.0;
        rom[5][0] = 2.0 * gom;
        
        rro1[0][2] = mn2 * gro1 + t * gro21;
        rro1[1][2] = 0.5 * (q2 + mp22 - t) * gro21;
        rro1[2][2] = (u - s) * gro21 / 4.0;
        rro1[3][2] = 0.0;
        rro1[4][2] = 0.0;
        rro1[5][2] = 2.0 * gro1;
        
        rom1[0][0] = mn2 * gom1 + t * gom21;
        rom1[1][0] = 0.5 * (q2 + mp22 - t) * gom21;
        rom1[2][0] = (u - s) * gom21 / 4.0;
        rom1[3][0] = 0.0;
        rom1[4][0] = 0.0;
        rom1[5][0] = 2.0 * gom1;
        
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 2; j++) {
                int k = j + 1;
                rro[i][j]  = 0.0;
                rom[i][k]  = 0.0;
                rro1[i][j] = 0.0;
                rom1[i][k] = 0.0;
            }
        }
        
        double   sregge = (s - u) / 2.0;
        
        double alpha_pi = 0.7 * (t - mp22);
        double alpha_ro = 0.55 + 0.8 * t;
        
        double reggepi = 0.7 * PI * (t - mp22) * pow(sregge, alpha_pi);
        reggepi = reggepi / (2.0 * sin(PI * alpha_pi) * g.gammaf(1.0 + alpha_pi));
        double reggero = 0.8 * PI * (t - 0.591)* pow(sregge, alpha_ro - 1.0);
        reggero = reggero / (2.0 * sin(PI * alpha_ro) * g.gammaf(alpha_ro));
        
        double reggepi_re =  (1.0 + cos(PI * alpha_pi)) * reggepi;
        double reggeb_re  = (-1.0 + cos(PI * alpha_pi)) * reggepi;
        double reggero_re = (-1.0 + cos(PI * alpha_ro)) * reggero;
        double reggea_re  =  (1.0 + cos(PI * alpha_ro)) * reggero;
        
        double[][] rpi_re = new double[6][3];
        double[][] rro_re = new double[6][3];
        double[][] rom_re = new double[6][3];
        
        for (int i = 0; i < 6; i++) {
            rpi_re[i][0] = 0.0;
            rpi_re[i][1] = re[i][1] * reggepi_re;
            rpi_re[i][2] = re[i][1] * reggeb_re;
            
            rro_re[i][0] = 0.0;
            rro_re[i][2] = rro1[i][2] * reggero_re;
            rro_re[i][1] = rro1[i][2] * reggea_re;
            
            rom_re[i][0] = rom1[i][0] * reggero_re;
            rom_re[i][1] = 0.0;
            rom_re[i][2] = 0.0;
        }
        
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 3; j++) {
                br[0][i][j] =    rom[i][j] +    rro[i][j] +     re[i][j] + rm[i][j];
                br[1][i][j] = rpi_re[i][j] + rom_re[i][j] + rro_re[i][j];
            }
        }
    }

}
