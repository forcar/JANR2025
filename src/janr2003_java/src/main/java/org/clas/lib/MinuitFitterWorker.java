package org.clas.lib;

import javax.swing.*;

import org.freehep.math.minuit.FunctionMinimum;
import org.freehep.math.minuit.MnMigrad;
import org.freehep.math.minuit.MnScan;
import org.freehep.math.minuit.MnUserParameters;

public class MinuitFitterWorker extends SwingWorker<FunctionMinimum, String> {
	
    private final JanrFunc fcn;
    private MnUserParameters upar;
	
    public MinuitFitterWorker(JanrFunc fcn) {
        this.fcn = fcn;
    }
	
    public Boolean FITPRINTOUT = true;
    
    @Override
    protected FunctionMinimum doInBackground() throws Exception {
    	// This runs in a background thread, so it won't freeze the GUI.
	    int npars = 30;
	    publish("Starting fit...");          
	    upar = new MnUserParameters();   
	        
	    for(int loop = 0; loop < npars; loop++){
	    	double val = fcn.j.start_value[loop];
	        upar.add(fcn.j.pname[loop],val,0.5,-0.2,0.2);
	        if(loop>=23) upar.fix(loop);
	    }
	           	        
	    MnScan  scanner = new MnScan(fcn,upar);
	    FunctionMinimum scanmin = scanner.minimize(); 
	        
	    System.err.println("******************");
	    System.err.println("*   SCAN RESULTS  *");
	    System.err.println("******************");
	    System.out.println("minimum : " + scanmin);
	    System.out.println("pars    : " + upar);
	    System.out.println(upar);
	    System.err.println("*******************************************");
	        
	    MnMigrad migrad = new MnMigrad(fcn, upar);  	        
	    FunctionMinimum min = migrad.minimize();
	    publish("Fit complete.");
	    MnUserParameters userpar = min.userParameters();
	        
	    for(int loop = 0; loop < npars; loop++){   	           
	         fcn.j.xnew[loop] = userpar.value(fcn.j.pname[loop]);
	        fcn.j.xnewe[loop] = userpar.error(fcn.j.pname[loop]);
	    }
	        
	    System.out.println(upar);
	    System.err.println("******************");
	    System.err.println("*   FIT RESULTS  *");
	    System.err.println("******************");
	            
	    System.err.println(min);
        
	    System.out.println(fcn.getBenchmarkString());
	        
		return min;
    }
    
    @Override
    protected void process(java.util.List<String> chunks) {
        // This runs on the Event Dispatch Thread (EDT) and updates the GUI.
        for (String chunk : chunks) {
            System.out.println(chunk); // For simple console updates
        }
    }

    @Override
    protected void done() {
        // This runs on the EDT after doInBackground() finishes.
        try {
            FunctionMinimum min = get();
            System.out.println("Final parameters:");
            System.out.println("Minimum FCN value: " + min.fval());
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
