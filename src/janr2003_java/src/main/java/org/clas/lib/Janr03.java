package org.clas.lib;

import java.util.List;
import org.jlab.groot.math.Func1D;

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
        case "F": janrFit(); break;
        case "T": janrOutput(); break;
        case "M": janrOutputMpole();      
        }
    	System.out.println("Janr03.janr_init() Finished\n");

    }
    
    public class JanrFunc extends Func1D {
    	
    	public JanrFunc(String name, double min, double max) {
    		super(name, min, max);    		
    	}
    	
    	@Override
    	public double evaluate(double x) {
    		double sum = 0.0;
    		return sum;
    	}
    	
    }

    private static void janrFit() {
    	
    }

    private static void janrOutput() {
    	
    }
    
    private static void janrOutputMpole() {
    	
    }
    
    private static void hjanrInit(List fname) {
    	
    }
    	
    public static void main(String[] args) { 
    	Janr03 j = new Janr03();
    };
     
}


