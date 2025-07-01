package org.clas.lib;

import java.io.File;
import java.io.IOException;
import java.util.Scanner;


public class Janr03 extends Constants {
	
    public static void main(String[] args) { };
    
    
    public void janr_main() {

        janrInit();
        
        System.out.println("JANR_MAIN: Begin reading tables");

        readMultables();
        reggeMultables();
        highMultables();
        
        System.out.println("JANR_MAIN: End reading tables");
        
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
    	String resifile;
    	String parifile;
   	
    	try (Scanner sc = new Scanner(System.in)) {

          String infile = sc.nextLine().trim();

          System.out.println("JANR_INIT: Reading the input file " + infile);

                try (Scanner fileScanner = new Scanner(new File(infile))) {
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
                    resifile = fileScanner.next();
                    parifile = fileScanner.next();
                    parofile = fileScanner.next();
                    datform = fileScanner.nextInt();
                    lmini = fileScanner.nextInt();
                    lmino = fileScanner.nextInt();
                    lminf = fileScanner.nextInt();
                    ChiMax = fileScanner.nextFloat();
                    for (int i = 0; i < 4; i++) fname[i] = fileScanner.next(); 
                }

                JanrInputRes loader1 = new JanrInputRes(resifile);
                loader1.loadResonanceData();

                try (Scanner parScanner = new Scanner(new File(parifile))) {
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
                        System.out.println("JANR_INIT: Wrong Number of parameters=" + parms);
                        return;
                    }
                }

                JanrInputSAID loader2 = new JanrInputSAID();
                loader2.loadMultipoleData();

                if (grfopt != "N") hjanrInit(fname);                

                qavg = 0.5f * (q2min + q2max);
                
                JanrIniPoint ini = new JanrIniPoint(qavg);
                MultAna       ma = new MultAna(qavg); 

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
 
}


