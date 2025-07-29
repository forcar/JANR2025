package org.clas.lib;

public class Janr03 extends Constants {
	
    Reader          r = new Reader();
    Writer          w = new Writer();
    JanrIniPoint  jip = new JanrIniPoint();

    public Janr03() {
    	janr_main();
    }
        
    public void janr_main() {

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
        case "F": janrFit(); break;
        case "T": janrOutput(); break;
        case "M": janrOutputMpole();      
        }

    }

    private static void janrFit() {
    	
    }

    private static void janrOutput() {
    	
    }
    
    private static void janrOutputMpole() {
    	
    }
    
    private static void hjanrInit(String[] fname) {
    	
    }
    	
    public static void main(String[] args) { 
    	Janr03 j = new Janr03();
    };
     
}


