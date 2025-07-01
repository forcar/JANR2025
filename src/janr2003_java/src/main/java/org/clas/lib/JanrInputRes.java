package org.clas.lib;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class JanrInputRes extends Constants {

    private double[][] resdat = new double[10][MAXres];

    String resifile;

    public JanrInputRes(String resifile) {
        this.resifile = resifile;
        initializeConstants();
    }

    private void initializeConstants() {
        pi = 2.0 * Math.acos(0.0);
        a2 = Math.sqrt(2.0);
        mn = 0.93827;
        mp = 0.13498;

        meta = 0.547;
        m2pion = 2.0 * mp;
        mn22 = mn * mn;
        mp22 = mp * mp;
        mn2 = 2.0 * mn;
    }

    public void loadResonanceData() throws IOException {

        try (BufferedReader br = new BufferedReader(new FileReader(resifile))) {
            String line;
            int k = 0;
            while ((line = br.readLine()) != null && k < MAXres) {
                // Split line into tokens (assuming whitespace or comma separated)
                String[] tokens = line.trim().split("[,\\s]+");
                if (tokens.length < 10) continue; // skip incomplete lines

                for (int i = 0; i < 10; i++) {
                    resdat[i][k] = Double.parseDouble(tokens[i]);
                }
                k++;
            }
        }

        for (int i = 0; i < MAXres; i++) {
            int ii = (int) resdat[0][i];
            if (ii < 0 || ii >= MAXres) continue; // prevent index out of bounds

            ires[ii]   = (int) resdat[1][i];
            jres[ii]   = (int) resdat[2][i];
            lres[ii]   = (int) resdat[3][i];
            mres[ii]   = resdat[4][i] / 1000.0;
            gres[ii]   = resdat[5][i] / 1000.0;
            eta[ii]    = resdat[6][i];
            xres[ii]   = resdat[7][i];
            langul[ii] = (int) resdat[8][i];
            lprime[ii] = (int) resdat[9][i];
        }
    }

    public void printResonances() {
        for (int i = 0; i < MAXres; i++) {
            System.out.printf("ires[%d]=%d, jres[%d]=%d, lres[%d]=%d, mres[%d]=%.6f, gres[%d]=%.6f, eta[%d]=%.6f, xres[%d]=%.6f, langul[%d]=%d, lprime[%d]=%d%n",
                i, ires[i], i, jres[i], i, lres[i], i, mres[i], i, gres[i], i, eta[i], i, xres[i], i, langul[i], i, lprime[i]);
        }
    }

    public void main(String[] args) {
        String filename = "resifile.dat"; 
        JanrInputRes loader = new JanrInputRes(filename);
        try {
            loader.loadResonanceData();
            loader.printResonances();
        } catch (IOException e) {
            System.err.println("Error reading resonance data: " + e.getMessage());
        }
    }
}

/*
Resonances with (****),(***),(**) are included
-(1)	P33(1232)
-(2)	D33(1700)
-(3)	F37(1950)
-(4)	P33(1600)
-(5)	S31(1620)
-(6)	F35(1905)
-(7)	D35(1930)
-(8)	P31(1925) missing resonance
-(9)	S31(1900)
-(10)	P33(1920)+ missing resonance (1975)
-(11)	P31(1910)
-(12)	D33(1940)
-(13)	F35(2000)
-(21)	P11(1440)
-(22) 	S11(1535)
-(23) 	D13(1520)
-(24) 	D13(1700)
-(25) 	S11(1650)
-(26) 	F15(1680)
-(27) 	P13(1720)
-(28) 	F17(1990)
-(29) 	D15(1675)
-(30)	P11(1710)
-(31)	P13(1900)+3 missing resonances (1870),(1980),(1955)
-(32)	D13(2080)
-(33)	F15(2000)+ missing resonance (1955)

Langul - photon ang. momenta in B-W 
Lprime(i) is the degree of photon momenta in Breit-Wigner formula 
for E(l+), E(l-):
E(l+):     Lprime(i)=l 
E(l-) l>1: Lprime(i)=l-2
E(l-) l=1: Lprime(i)=l for E(l-)
*/      

