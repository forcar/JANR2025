package org.clas.lib;

import org.freehep.math.minuit.FCNBase;

public class JanrFunc implements FCNBase {
	
    private int       numberOfCalls = 0;
    private long      startTime     = 0L;
    private long      endTime       = 0L;
    Janr03            j             = null;
	
	public JanrFunc(Janr03 j) {
		this.j = j;	
        startTime  = System.currentTimeMillis();
	}
	
	@Override
	public double valueOf(double xval[]) {
		for (int i=0; i<xval.length; i++) j.xnew[i] = xval[i];
		double chi2 = 0.0;
    	for (int i=0; i<j.nplt; i++ ) chi2+=getChi2(i,0,j.ndims[i][0],0,j.ndims[i][1],0,j.ndims[i][2]); //loop over reaction channels
        numberOfCalls++;
        System.out.println(numberOfCalls+" "+chi2+" "+j.xnew[0]+" "+j.xnew[1]+" "+j.xnew[2]+" "+j.xnew[3]);
        endTime = System.currentTimeMillis();
		return chi2;
	}
	
	public double getChi2(int ifi, int w1, int w2, int c1, int c2, int f1, int f2) {
		
    	double[] dumd = new double[100];
    	double[] dum  = new double[100];
    	double[] dume = new double[100];
    	double[]  obs = new double[100];
    	double[] obs1 = new double[100];
		
    	int indx=0, iplt=0, k=0, ndf=0;
		double chi=0.0, chi2 = 0.0;
    	double xmx=-100;
		
//    	if (w1!=w2) iplt=2;
//     	if (c1!=c2) iplt=1; 
//    	if (f1!=f2) iplt=0;
    	
    	for (int iw=w1; iw<w2+1; iw++) {
    		int iww = iw*j.ndims[ifi][1]*j.ndims[ifi][0];
    		for (int ic=c1; ic<c2+1; ic++) {
    			int icc = ic*j.ndims[ifi][0];
    			for (int iff=f1; iff<f2+1; iff++) {
    				
    				 indx = iww+icc+iff;
 				 
//    				 if(xx[ifi][iplt][indx] < xmx) break;  
    				 
//    				 double  x = (float) xx[ifi][iplt][indx];
    				 double  y =  j.crs[ifi][indx];
    				 double ye = j.crse[ifi][indx]; 
    				 
    				 if(ye<0.00000001) continue;
     				 
    				 int i1 = (int) j.xx[ifi][6][indx]-1;
    				 int i2 = (int) j.xx[ifi][5][indx]-1;

    				 double dthe = j.xx[ifi][1][indx]; // could be either theta or cos(theta)

    				 j.loadpar(2); //new fit  
    				 j.jr.janrRun(j.xx[ifi][2][indx], 
						          j.xx[ifi][3][indx],
						          j.xx[ifi][4][indx],
						          dthe>1.0 ? Math.cos(Math.toRadians(dthe)):dthe,
						          j.xx[ifi][0][indx],
						          true);
				     double yv = j.robs[i1][i2]; 
				     double diff = (yv-y)/ye;
//	                 xmx=xplt[k];    	                 
	                 chi = diff * diff;
	                 chi2+=chi;
	                 k++;
    			}
    		}
    	}
    	ndf = k;
		return chi2;
	}
	
    public String getBenchmarkString(){
        StringBuilder str = new StringBuilder();
        double time = (double) (endTime-startTime);
        str.append(String.format("[fit-benchmark] Time = %.3f , Iterrations = %d"
                , time/1000.0,
                this.numberOfCalls));
        return str.toString();
    }
}

