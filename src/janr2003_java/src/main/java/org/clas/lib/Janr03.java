package org.clas.lib;

import java.util.List;

public class Janr03 extends Constants {
	
    Reader          r = new Reader();
    Writer          w = new Writer();
    JanrIniPoint  jip = new JanrIniPoint();
    JanrRun        jr = new JanrRun();

    public Janr03() {  
        janr_init();
    }
        
    public void janr_init() {

    	System.out.println("Janr03.janr_init() Starting\n");
    	
        r.InputFile(infile);
        r.ResonanceData(resifile);
        r.ParFile(parifile); 
        r.Multipoles();

        if (grfopt != "N") hjanrInit(fname);                

        double qavg = 0.5f * (q2min + q2max);
        
        jip.janrIniPoint(qavg); 
        
        w.MultAna(qavg);        
        r.MultTables();
        r.ReggeTables();
        r.HighMultTables();  
        
        w.generateTables();
        
        switch (fitopt) {
//        case "F": jf.fit("VE"); break;
        case "T": janrOutput(); break;
        case "M": janrOutputMpole();      
        }
        
    	System.out.println("Janr03.janr_init() Finished\n");
    }
    
    public void loadpar(int n) {
    	
        double[] xload = new double[36]; 
        
        switch (n) {
        case 1: for (int i = 0; i <= 35; i++)  xload[i] = start_value[i]; break;
        case 2: for (int i = 0; i <= 35; i++)  xload[i] = xnew[i];        break;
        case 3: for (int i = 0; i <= 35; i++)   xnew[i] = start_value[i]; break;
        case 4: for (int i = 0; i <= 35; i++) start_value[i] = xnew[i];
        }

        cm1  = xload[0];
        ce1  = xload[1];
        cs1  = xload[2];
        cmp1 = xload[3];
        csp1 = xload[4];
        cep2 = xload[5];
        csp2 = xload[6];
        cmp3 = xload[7];
        cep3 = xload[8];
        csp3 = xload[9];
        cep4 = xload[10];
        csp4 = xload[11];
        cep5 = xload[12];
        csp5 = xload[13];
        cmp6 = xload[14];
        cep6 = xload[15];
        csp6 = xload[16];
        cmp7 = xload[17];
        cep7 = xload[18];
        csp7 = xload[19];
        cmp8 = xload[20];
        cep8 = xload[21];
        csp8 = xload[22];
        cmp9 = xload[23];
        cep9 = xload[24];
        cmp10 = xload[25];
        cep10 = xload[26];
        cmp11 = xload[27];
        cep11 = xload[28];
        cm2 = xload[29];
        ce2 = xload[30];
        cs2 = xload[31];
        cm3 = xload[32];
        ce3 = xload[33];
        cs3 = xload[34];
        cspa = xload[35];
    }

    private static void janrOutput() {
    	
    }
    
    private static void janrOutputMpole() {
    	
    }
    
    private static void hjanrInit(List fname) {
    	
    }
    	
    public static void main(String[] args) {     	
    	
    }
     
}


