package org.clas.lib;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.NoSuchElementException;
import java.util.Scanner;

public class Reader extends Constants {
	
	boolean debug=false;
	
	public Reader() {
		initializeConstants(); 
	}
	
    public void initializeConstants() { 
        System.out.println("Janr03.Reader.initializeConstants()");
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
	
    public void InputFile(String file) {
        System.out.println("Janr03.Reader.InputFile("+file+")");
        try {
            Scanner fileScanner = new Scanner(new File(file));
        	while (fileScanner.hasNextLine()) {
               fitopt = fileScanner.next();  
               grfopt = fileScanner.next();  
                   n1 = fileScanner.nextInt();  
                   n2 = fileScanner.nextInt();  
                nebac = fileScanner.nextInt();  
             obsv2fit = fileScanner.nextInt();
                ebeam = fileScanner.nextFloat();
                 wmin = fileScanner.nextFloat();
                 wmax = fileScanner.nextFloat();
                wstep = fileScanner.nextFloat();
                q2min = fileScanner.nextFloat();
                q2max = fileScanner.nextFloat();
               q2step = fileScanner.nextFloat();
               cosmin = fileScanner.nextFloat();
               cosmax = fileScanner.nextFloat();
              cosstep = fileScanner.nextFloat();
               phimin = fileScanner.nextFloat();
               phimax = fileScanner.nextFloat();
              phistep = fileScanner.nextFloat();
             IntAccur = fileScanner.nextFloat();
              inffile = fileScanner.next();
             maxcalls = fileScanner.nextInt();
            tolorance = fileScanner.nextFloat();
              errorup = fileScanner.nextFloat();
                parms = fileScanner.nextInt();
             resifile = janrPath+fileScanner.next();
             parifile = janrPath+fileScanner.next();
             parofile = fileScanner.next();
              datform = fileScanner.nextInt();
                lmini = fileScanner.nextInt();
                lmino = fileScanner.nextInt();
                lminf = fileScanner.nextInt();
               ChiMax = fileScanner.nextFloat(); 
            for (int i = 0; i < 10; i++) fname.add(fileScanner.next());
        	} 
        	fileScanner.close();
        } catch (NoSuchElementException e) {
        } catch (FileNotFoundException e) {
        	e.printStackTrace();
        }

        nplt = fname.size();

    }
       
    public void ParFile(String file) {
        System.out.println("Janr03.Reader.ParFile("+file+")");
    	try (Scanner parScanner = new Scanner(new File(file))) {
            if (parms > 0 && parms <= MAXpar) {
                for (int i = 0; i < parms; i++) {
                          plist[i] = parScanner.nextInt();
                          pname[i] = parScanner.next();
                    start_value[i] = parScanner.nextFloat(); //this must be read in before calculations!  
                      step_size[i] = parScanner.nextFloat();
                       par_stat[i] = parScanner.nextInt();
                        low_lim[i] = parScanner.nextFloat();
                         up_lim[i] = parScanner.nextFloat();
                       res_name[i] = parScanner.next();
            if(debug) System.out.println(plist[i]+" "+pname[i]+" "+start_value[i]+" "+step_size[i]+" "+par_stat[i]+" "+low_lim[i]+" "+up_lim[i]);
                }
            } else {
                System.out.println("janrParFile: Wrong Number of parameters=" + parms);
                return;
            }
        } catch (IOException e) {
            e.printStackTrace();
        } 
    	
         cm1 = start_value[0];
         ce1 = start_value[1];
         cs1 = start_value[2];
        cmp1 = start_value[3];
        csp1 = start_value[4];
        cep2 = start_value[5];
        csp2 = start_value[6];
        cmp3 = start_value[7];
        cep3 = start_value[8];
        csp3 = start_value[9];
        cep4 = start_value[10];
        csp4 = start_value[11];
        cep5 = start_value[12];
        csp5 = start_value[13];
        cmp6 = start_value[14];
        cep6 = start_value[15];
        csp6 = start_value[16];
        cmp7 = start_value[17];
        cep7 = start_value[18];
        csp7 = start_value[19];
        cmp8 = start_value[20];
        cep8 = start_value[21];
        csp8 = start_value[22];
        cmp9 = start_value[23];
        cep9 = start_value[24];
       cmp10 = start_value[25];
       cep10 = start_value[26];
       cmp11 = start_value[27];
       cep11 = start_value[28];
         cm2 = start_value[29];
         ce2 = start_value[30];
         cs2 = start_value[31];
         cm3 = start_value[32];
         ce3 = start_value[33];
         cs3 = start_value[34];
        cspa = start_value[35];
    }
    
    public void ResonanceData(String file) {
        System.out.println("Janr03.Reader.ResonanceData("+file+")");

        double[][] resdat = new double[10][MAXres];

        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line; int k = 0;
            while (k < Nres) {
                line = br.readLine();
                if (line.trim().length() ==  0) continue;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 10) break;
                for (int i = 0; i < 10; i++) {
                    resdat[i][k] = Double.parseDouble(tokens[i]);
                }
                k++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        
        // Format of input resonance file (e.g.-res-q0.4.inp) 
        //  ID  i j l m      w   b    x     l1 l2
        //  ID: id number used in JANR  
        //  i: isospin of resonance (3=3/2 1=1/2) 
        //  j: total angular momentum of resonance (3=3/2 1=1/2)
        //  l: angular momentum of decay pion (l=0,1,2...)
        //  m: mass of resonance (MeV)
        //  w: width of resonance (MeV)
        //  b: branching ratio for pi-N channel
        //  x: cutoff used in Breit-Wigner form factors
        // l1: angular momentum of photon multipole used in B-W threshold factors
        // l2: power used in B-W pi-N threshold factors
      
        //-Resonances with (****),(***),(**) are included
        //-(1)	P33(1232)
        //-(2)	D33(1700)
        //-(3)	F37(1950)
        //-(4)	P33(1600)
        //-(5)	S31(1620)
        //-(6)	F35(1905)
        //-(7)	D35(1930)
        //-(8)	P31(1925) missing resonance
        //-(9)	S31(1900)
        //-(10)	P33(1920)+ missing resonance (1975)
        //-(11)	P31(1910)
        //-(12)	D33(1940)
        //-(13)	F35(2000)
        
        //-(21)	P11(1440)
        //-(22) S11(1535)
        //-(23) D13(1520)
        //-(24) D13(1700)
        //-(25) S11(1650)
        //-(26) F15(1680)
        //-(27) P13(1720)
        //-(28) F17(1990)
        //-(29) D15(1675)
        //-(30)	P11(1710)
        //-(31)	P13(1900)+3 missing resonances (1870),(1980),(1955)
        //-(32)	D13(2080)
        //-(33)	F15(2000)+ missing resonance (1955)

        // Langul: photon angular momentum in B-W 
        // Lprime(i): degree of photon momenta in Breit-Wigner formula 
        // For E(l+), E(l-):
        // E(l+):     Lprime(i)=l 
        // E(l-) l>1: Lprime(i)=l-2
        // E(l-) l=1: Lprime(i)=l for E(l-)   
        
        for (int i = 0; i < MAXres; i++) { //[0] not used since indices are hard-coded later
             int ii  = (int) resdat[0][i];
             if (ii < 0 || ii >= MAXres) continue; 
            ires[ii] = (int) resdat[1][i];
            jres[ii] = (int) resdat[2][i];
            lres[ii] = (int) resdat[3][i];
            mres[ii] =       resdat[4][i] / 1000.0;
            gres[ii] =       resdat[5][i] / 1000.0;
             eta[ii] =       resdat[6][i];
            xres[ii] =       resdat[7][i];
          langul[ii] = (int) resdat[8][i];
          lprime[ii] = (int) resdat[9][i];
        }
        
     }
    
	public void Multipoles() {
        System.out.println("Janr03.Reader.Multipoles");
		loadMultipoleData("mp33",51,6,mp33r,mp33i);
		loadMultipoleData( "s11",51,9,s11r,s11i);
		loadMultipoleData( "s31",51,9,s31r,s31i);
		loadMultipoleData( "p11",51,9,p11r,p11i);
		loadMultipoleData( "p13",51,9,p13r,p13i);
		loadMultipoleData( "p31",51,9,p31r,p31i);
		loadMultipoleData( "p33",51,9,p33r,p33i);
		loadMultipoleData( "d13",51,9,d13r,d13i);
		loadMultipoleData( "d15",51,9,d15r,d15i);
		loadMultipoleData( "d33",51,9,d33r,d33i);
		loadMultipoleData( "f15",51,9,f15r,f15i);
		loadMultipoleData( "f37",51,9,f37r,f37i);
		loadMultipoleData( "pim",13,3,pimr,pimi);
		loadMultipoleData( "sim",13,3,simr,simi);
	}
    
	public void loadMultipoleData(String mpol, int nmax, int len, double re[], double im[]) {
		boolean mp = mpol=="mp33", pm = mpol=="pim"||mpol=="sim";
        try (BufferedReader br = new BufferedReader(new FileReader(janrPath+"multipols/"+mpol+".dat"))) {
        	String line; int i=0; 
            while(i < nmax) {
                line = br.readLine();
                if (line.trim().length() ==  0) continue;
                String[] tokens = line.trim().split("\\s+");
                if (tokens.length < len) break;
                WSa[i]= Double.parseDouble(tokens[0]) / 1000.0;
                re[i] = Double.parseDouble(tokens[mp||pm?1:5])/(mp?52.437:1.0);
                im[i] = Double.parseDouble(tokens[mp||pm?2:6])/(mp?52.437:1.0); 
                // at2, at3, at4 are not used 
                i++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }		
	}

	public void MultTables() {
        System.out.println("Janr03.Reader.MultTables");
        try {
            for (int L = 0; L < 4; L++) {
                BufferedReader[] readers = new BufferedReader[6];
                try {
                    // Open the six files for given L and skip first two header lines
                    for (int j = 0; j < 6; j++) {
                        int n = L * 6 + j;
                        readers[j] = new BufferedReader(new FileReader(janrPath+"multipols/" + fnamr[n]));
                        readers[j].readLine();
                        readers[j].readLine();
                    }
                    // Read data lines
                    for (int iw = 0; iw < Maxmpoints; iw++) {
                        String line;
                        String[] tok;

                        // rm"L"p.dat -> amp
                        line = readers[0].readLine();
                        if (line == null) break;
                        tok = line.trim().split("\\s+");
                        double w   = Double.parseDouble(tok[0]);
                        double m23 = Double.parseDouble(tok[2]);
                        double m12 = Double.parseDouble(tok[3]);
                        double ms  = Double.parseDouble(tok[4]);
                        amp[iw][L][0] = ms;
                        amp[iw][L][1] = m12;
                        amp[iw][L][2] = m23;

                        // rm"L"m.dat -> amm
                        line = readers[1].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        amm[iw][L][0] = ms;
                        amm[iw][L][1] = m12;
                        amm[iw][L][2] = m23;

                        // re"L"p.dat -> aep
                        line = readers[2].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        aep[iw][L][0] = ms;
                        aep[iw][L][1] = m12;
                        aep[iw][L][2] = m23;

                        // re"L"m.dat -> aem
                        line = readers[3].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        aem[iw][L][0] = ms;
                        aem[iw][L][1] = m12;
                        aem[iw][L][2] = m23;

                        // rs"L"p.dat -> asp
                        line = readers[4].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        asp[iw][L][0] = ms;
                        asp[iw][L][1] = m12;
                        asp[iw][L][2] = m23;

                        // rs"L"m.dat -> asm
                        line = readers[5].readLine();
                        tok  = line.trim().split("\\s+");
                        w   = Double.parseDouble(tok[0]);
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        asm[iw][L][0] = ms;
                        asm[iw][L][1] = m12;
                        asm[iw][L][2] = m23;

                        // record w (should be same across files)
                        wtab[iw] = w;
                    }
                } finally {
                    // Close all readers
                    for (BufferedReader br : readers) {
                        if (br != null) br.close();
                    }
                }
            }
            
        } catch (IOException e) {
            e.printStackTrace();
        }		

	}
	
	public void ReggeTables() {
        System.out.println("Janr03.Reader.ReggeTables");
        try {
            for (int L = 0; L < 4; L++) {
                BufferedReader[] readers = new BufferedReader[6];
                try {
                    // Open the six files for given L and skip first two header lines
                    for (int j = 0; j < 6; j++) {
                        int n = L * 6 + j;
                        readers[j] = new BufferedReader(new FileReader(janrPath+"multipols/" + fnamrr[n]));
                        readers[j].readLine();
                        readers[j].readLine();
                    }
                    // Read data lines
                    for (int iw = 0; iw < Maxmpoints; iw++) {
                        String line;
                        String[] tok;

                        // rm"L"p.dat -> amp
                        line = readers[0].readLine();
                        if (line == null) break;
                        tok = line.trim().split("\\s+");
                        double  w  = Double.parseDouble(tok[0]);
                        double m23 = Double.parseDouble(tok[2]);
                        double m12 = Double.parseDouble(tok[3]);
                        double ms  = Double.parseDouble(tok[4]);
                        ramp[iw][L][0] = ms;
                        ramp[iw][L][1] = m12;
                        ramp[iw][L][2] = m23;

                        // rm"L"m.dat -> amm
                        line = readers[1].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        ramm[iw][L][0] = ms;
                        ramm[iw][L][1] = m12;
                        ramm[iw][L][2] = m23;

                        // re"L"p.dat -> aep
                        line = readers[2].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        raep[iw][L][0] = ms;
                        raep[iw][L][1] = m12;
                        raep[iw][L][2] = m23;

                        // re"L"m.dat -> aem
                        line = readers[3].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        raem[iw][L][0] = ms;
                        raem[iw][L][1] = m12;
                        raem[iw][L][2] = m23;

                        // rs"L"p.dat -> asp
                        line = readers[4].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        rasp[iw][L][0] = ms;
                        rasp[iw][L][1] = m12;
                        rasp[iw][L][2] = m23;

                        // rs"L"m.dat -> asm
                        line = readers[5].readLine();
                        tok  = line.trim().split("\\s+");
                        w   = Double.parseDouble(tok[0]);
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        rasm[iw][L][0] = ms;
                        rasm[iw][L][1] = m12;
                        rasm[iw][L][2] = m23;

                        // record w (should be same across files)
                        wtab[iw] = w;
                    }
                } finally {
                    // Close all readers
                    for (BufferedReader br : readers) {
                        if (br != null) br.close();
                    }
                }
            }
            
        } catch (IOException e) {
            e.printStackTrace();
        }
				
	}	
	
	public void HighMultTables() {
        System.out.println("Janr03.Reader.HighMultTables");
        try {
            for (int L = 4; L < 6; L++) {
                BufferedReader[] readers = new BufferedReader[6];
                try {
                    // Open the six files for given L and skip first two header lines
                    for (int j = 0; j < 6; j++) {
                        int n = (L-4) * 6 + j;
                        readers[j] = new BufferedReader(new FileReader(janrPath+"multipols/" + fnamhigh[n]));
                        readers[j].readLine();
                        readers[j].readLine();
                    }
                    // Read data lines
                    for (int iw = 0; iw < Maxmpoints; iw++) {
                        String line;
                        String[] tok;

                        // rm"L"p.dat -> amp
                        line = readers[0].readLine();
                        if (line == null) break;
                        tok = line.trim().split("\\s+");
                        double  w  = Double.parseDouble(tok[0]);
                        double m23 = Double.parseDouble(tok[2]);
                        double m12 = Double.parseDouble(tok[3]);
                        double ms  = Double.parseDouble(tok[4]);
                        amp[iw][L][0] = ms;
                        amp[iw][L][1] = m12;
                        amp[iw][L][2] = m23;

                        // rm"L"m.dat -> amm
                        line = readers[1].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        amm[iw][L][0] = ms;
                        amm[iw][L][1] = m12;
                        amm[iw][L][2] = m23;

                        // re"L"p.dat -> aep
                        line = readers[2].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        aep[iw][L][0] = ms;
                        aep[iw][L][1] = m12;
                        aep[iw][L][2] = m23;

                        // re"L"m.dat -> aem
                        line = readers[3].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        aem[iw][L][0] = ms;
                        aem[iw][L][1] = m12;
                        aem[iw][L][2] = m23;

                        // rs"L"p.dat -> asp
                        line = readers[4].readLine();
                        tok  = line.trim().split("\\s+");
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        asp[iw][L][0] = ms;
                        asp[iw][L][1] = m12;
                        asp[iw][L][2] = m23;

                        // rs"L"m.dat -> asm
                        line = readers[5].readLine();
                        tok  = line.trim().split("\\s+");
                        w   = Double.parseDouble(tok[0]);
                        m23 = Double.parseDouble(tok[2]);
                        m12 = Double.parseDouble(tok[3]);
                        ms  = Double.parseDouble(tok[4]);
                        asm[iw][L][0] = ms;
                        asm[iw][L][1] = m12;
                        asm[iw][L][2] = m23;

                        // record w (should be same across files)
                        wtab[iw] = w;
                    }
                } finally {
                    // Close all readers
                    for (BufferedReader br : readers) {
                        if (br != null) br.close();
                    }
                }
            }
            
        } catch (IOException e) {
            e.printStackTrace();
        }
				
	
	}
	
    public static void main(String[] args) {
    	Reader r = new Reader(); 
    	r.debug = false;
        r.InputFile(infile);
        r.ResonanceData(resifile);
        r.ParFile(parifile); 
        r.Multipoles();
    	r.MultTables();
    	r.ReggeTables();
    	r.HighMultTables();
    }

}
