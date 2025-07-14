package org.clas.lib;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

public class JanrInputRes extends Constants {

    private double[][] resdat = new double[10][MAXres];

    public JanrInputRes() {
        initializeConstants();
    }

    public void initializeConstants() {
        System.out.println("JanrInputRes.initializeConstants()");
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

    public void loadResonanceData(String file) {
        System.out.println("JanrInputRes.loadResonanceData("+file+")");
        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line; int k = 0;
            while (k < Nres) {
                line = br.readLine();
                if (line.trim().length() ==  0) continue;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 10) break;
                for (int i = 0; i < 10; i++) {
                    resdat[i][k] = Double.parseDouble(tokens[i]);
//                    System.out.println(i+" "+k+" "+resdat[i][k]);
                }
                k++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        for (int i = 0; i < MAXres; i++) { //[0] not used since indices are hard-coded later
            int ii = (int) resdat[0][i];
            if (ii < 0 || ii >= MAXres) continue; 
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
        
        System.out.println("MAXres = "+MAXres);
//        for (int i=0; i<MAXres; i++) System.out.println(i+" "+(int)resdat[1][i]+" "+ires[i]+" "+jres[i]+" "+(int)lres[i]+" "+mres[i]);        
//        for (int i=0; i<MAXres; i++) System.out.println(i+" "+(int)resdat[0][i]+" "+(int)resdat[1][i]+" "+(int)ires[i]+" "+(int)jres[i]+" "+(int)lres[i]+" "+mres[i]);



    }

    public static void main(String[] args) {   	
       JanrInputRes jir = new JanrInputRes();
       jir.loadResonanceData("run/res-q0.4.inp");
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

