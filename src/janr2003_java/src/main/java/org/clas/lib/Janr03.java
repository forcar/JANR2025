package org.clas.lib;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;

public class Janr03 extends Constants {
	
	Scanner filescanner;
    String infile   = janrPath+"janr_input_e1b";
    String resfile  = janrPath+"res-q0.4.inp";
    String parifile = null, resifile=null;
    JanrInputRes  jir = new JanrInputRes();

    public Janr03() {
    	janr_main();
    }
        
    public void janr_main() {

        janrInit();
        
        System.out.println("Janr03.janr_main(): Begin reading tables");
        readMultables();
        reggeMultables();
        highMultables();        
        System.out.println("Janr03.janr_main(): End reading tables");
        
        if (fitopt.equals("F")) {
            janrFit();
        }
        if (fitopt.equals("T")) {
            janrOutput();
        }
        if (fitopt.equals("M") || fitopt.equals("T")) {
            janrOutputMpole();
        }	
    }
    

    public void janrInit() {
    	
    	double qavg;
    	 	
        System.out.println("Janr03.janrInit()");
        
        janrInFile(infile);
        jir = new JanrInputRes();                jir.loadResonanceData(resfile);
        janrParFile(parifile); 
        JanrInputSAID jis = new JanrInputSAID(); jis.loadMultipoles();

        if (grfopt != "N") hjanrInit(fname);                

        qavg = 0.5f * (q2min + q2max);                
        JanrIniPoint ini = new JanrIniPoint(qavg);
        
        MultAna       ma = new MultAna(); ma.multAna(qavg);

    }
    
    public void janrInFile(String file) {
        System.out.println("janrInFile("+file+")");
        try (Scanner fileScanner = new Scanner(new File(file))) {
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
            phimin= fileScanner.nextFloat();
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
            for (int i = 0; i < 4; i++) fname[i] = fileScanner.next(); 
        } catch (IOException e) {
            e.printStackTrace();
        }    	    	
    }
    
    public void janrParFile(String file) {
        System.out.println("janrParFile("+file+")");
    	try (Scanner parScanner = new Scanner(new File(file))) {
            if (parms > 0 && parms <= MAXpar) {
                for (int i = 0; i < parms; i++) {
                    plist[i] = parScanner.nextInt();
                    pname[i] = parScanner.next();
                    start_value[i] = parScanner.nextFloat();
                    step_size[i] = parScanner.nextFloat();
                    par_stat[i] = parScanner.nextInt();
                    low_lim[i] = parScanner.nextFloat();
                    up_lim[i] = parScanner.nextFloat();
                }
            } else {
                System.out.println("janrParFile: Wrong Number of parameters=" + parms);
                return;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }    	
    }



    private static void readMultables() {

    }    
    
    private static void reggeMultables() {

    }

    private static void highMultables() {

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


