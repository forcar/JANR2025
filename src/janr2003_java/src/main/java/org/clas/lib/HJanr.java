package org.clas.lib;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;

public class HJanr extends JanrMonitor {
	
	int[] c={8,3,6,7}, f={2,1,2,2}, w={5,18,8,10}; //default bins (coscm,phicm,w) displayed in DATA tab

    public HJanr(String name) {
        super(name);
        setRunNumber(0);
        dgmActive = true;
        setJanrTabNames(0,"DATA","MINUIT");
        setJanrTabNames(1,"FITC");
        useWCFButtons(true);
        use1234Buttons(true);
        useBinSliderPane(true);
        init();
        hjanr_init();
        createHistos(getRunNumber());
        hjanr_loadpar(3);
        plotHistos(getRunNumber(),0,0);
    }
    
    public void hjanr_init() { 	
    	int i=0;
    	hjanr_loadpar(1);
     	for (String name : j.fname) hjanr_input(i++,j.janrPath+name); //loop over all reactions
    }
    
    public void hjanr_input(int iF, String file) {
    	
        String[] lab1 = {"(#pi^0)","(#pi^+)"};
        String[] lab2 = {"  d#sigma","ALT^/","ALT","#Sigma"};
        
        double[]    xxk = new double[3];
        boolean[] iflag = new boolean[3];

        Arrays.fill(iflag, true);
        Arrays.fill(xxk, -1000.0);

        System.out.println("HJANR_INPUT: Reading file " + file);        

        try (BufferedReader br = new BufferedReader(new FileReader(file))) {
            String line;
            int k = 0;
            while ((line = br.readLine()) != null && k < j.MAX_K) { 
            	String[] tokens = line.trim().split("\\s+");
                if (tokens.length < 12) break;
                j.xx[iF][2][k] = Double.parseDouble(tokens[0]); //W
                j.xx[iF][3][k] = Double.parseDouble(tokens[1]); //q2
                j.xx[iF][4][k] = Double.parseDouble(tokens[2]); //eps
                j.xx[iF][1][k] = Double.parseDouble(tokens[3]); //coscm
                j.xx[iF][0][k] = Double.parseDouble(tokens[4]); //phicm
                j.crs[iF][k]   = Double.parseDouble(tokens[5]); //crs
                j.crse[iF][k]  = Double.parseDouble(tokens[6]); //crse
                
                j.xx[iF][5][k] = Double.parseDouble(tokens[10]); // 1:p-pi0 2:n-pi+
                j.xx[iF][6][k] = Double.parseDouble(tokens[11]); // 1:crs   2:alt'

                if (j.xx[iF][2][k] == 0.0) break;

                for (int i = 0; i < 3; i++) { //determine number of phi,costh,w bins
                    if (j.xx[iF][i][k] > xxk[i] && iflag[i]) {
                        j.ndims[iF][i] = j.ndims[iF][i] + 1;
                    } else if (j.xx[iF][i][k] < xxk[i]) {
                        iflag[i] = false;
                    }
                    xxk[i] = j.xx[iF][i][k];
                   }
                k++;
            }
        } catch (IOException e) {
            System.err.println("Error reading file: " + file);
            e.printStackTrace();
        }

        int lab2Index = (int) j.xx[iF][6][0] - 1;
        int lab1Index = (int) j.xx[iF][5][0] - 1;
        String lab2Str = (lab2Index >= 0 && lab2Index < lab2.length) ? lab2[lab2Index] : "";
        String lab1Str = (lab1Index >= 0 && lab1Index < lab1.length) ? lab1[lab1Index] : "";
        
        j.lab[iF] = lab2Str + lab1Str; //store observable label
    }
    
    @Override 
    public void createHistos(int run) {
    	histosExist = true; 
    	createDATA();
    	createFITC();
    }
    
    public void createDATA() {
    	dgm[0].add("DATA", j.nplt,3,0,0,0); 
    	for (int i=0; i<j.nplt; i++) {dgm[0].makeGE("WDAT" +i, 0, j.lab[i],"", i==0?    "W":"", 1,5);
    	                              dgm[0].makeGE("WSEL" +i,-2, "",      "", i==0?     "":"", 3,4,1,3);
                                      dgm[0].makeGE("WFIT1"+i,-2, "",      "", i==0?    "W":"", 1,1,1,1);
                                      dgm[0].makeGE("WFIT2"+i,-2, "",      "", i==0?    "W":"", 2,1,1,2);}
     	for (int i=0; i<j.nplt; i++) {dgm[0].makeGE("CDAT" +i, 0, "",      "", i==0?"COSCM":"", 1,5);
                                      dgm[0].makeGE("CSEL" +i,-2, "",      "", i==0?     "":"", 3,4,1,3);
    	                              dgm[0].makeGE("CFIT1"+i,-2, "",      "", i==0?"COSCM":"", 1,1,1,1);
                                      dgm[0].makeGE("CFIT2"+i,-2, "",      "", i==0?"COSCM":"", 2,1,1,2);}
        for (int i=0; i<j.nplt; i++) {dgm[0].makeGE("FDAT" +i, 0, "",      "", i==0?"PHICM":"", 1,5);
                                      dgm[0].makeGE("FSEL" +i,-2, "",      "", i==0?     "":"", 3,4,1,3);
       	                              dgm[0].makeGE("FFIT1"+i,-2, "",      "", i==0?"PHICM":"", 1,1,1,1);
       	                              dgm[0].makeGE("FFIT2"+i,-2, "",      "", i==0?"PHICM":"", 2,1,1,2);
    	                              }
    }
    
    public void createFITC() {
    	dgm[1].add("FITC", 1,1,0,0,0);   	
    	String titx = j.pname[0]+"="+String.format("%.2f",j.start_value[0])+" ";
    	for (int i=1; i<30; i++) titx+=" "+j.pname[i]+"="+String.format("%.2f",j.start_value[i])+"    ";
       	dgm[1].makeGE("COF1", -1, 0, 31, "", titx," ",1,5);     dgm[1].cc("COF1", false, false, 0.5f, 30.5f, -4f, 4f, 0f, 0f);
       	dgm[1].makeGE("SEL1", -2, 0, 31, "", titx," ",2,4,1); 
    }

    @Override       
    public void plotHistos(int run, int nt, int bin) {
    	if(!histosExist) return;
    	if(nt>=0) hjanr_bin_pick(nt, bin);
    	hjanr_plot_data();
    	dgm[0].drawGroup("DATA",0,0,0);    	
    	hjanr_plot_par(); 
    	dgm[1].drawGroup("FITC",0,0,0);
    }
    
    public void hjanr_pick() {
    	
    }
    
    public void hjanr_bin_pick(int nt, int bin) { //Pick new data points in (W,coscm,phicm)
    	switch (nt) {
    	case 0: binslider.setMaximum(j.ndims[0][2]-1); w[0] = bin; break; 
    	case 4: binslider.setMaximum(j.ndims[0][1]-1); c[0] = bin; break;
    	case 8: binslider.setMaximum(j.ndims[0][0]-1); f[0] = bin; break; 
    	case 1: binslider.setMaximum(j.ndims[1][2]-1); w[1] = bin; break; 
    	case 5: binslider.setMaximum(j.ndims[1][1]-1); c[1] = bin; break;
    	case 9: binslider.setMaximum(j.ndims[1][0]-1); f[1] = bin; break;
    	case 2: binslider.setMaximum(j.ndims[2][2]-1); w[2] = bin; break;
    	case 6: binslider.setMaximum(j.ndims[2][1]-1); c[2] = bin; break;
    	case 10:binslider.setMaximum(j.ndims[2][0]-1); f[2] = bin; break;
    	case 3: binslider.setMaximum(j.ndims[3][2]-1); w[3] = bin; break;
    	case 7: binslider.setMaximum(j.ndims[3][1]-1); c[3] = bin; break;
    	case 11:binslider.setMaximum(j.ndims[3][0]-1); f[3] = bin;
    	}
    }
 
    public void hjanr_plot_data() { 
    	for (int i=0; i<j.nplt; i++ ) hjanr_plot_fit(i,f[i],f[i],c[i],c[i],0,j.ndims[i][2],w[i]); //w
    	for (int i=0; i<j.nplt; i++ ) hjanr_plot_fit(i,f[i],f[i],0,j.ndims[i][1],w[i],w[i],c[i]); //coscm
    	for (int i=0; i<j.nplt; i++ ) hjanr_plot_fit(i,0,j.ndims[i][0],c[i],c[i],w[i],w[i],f[i]); //phicm
    	
    }
    
    public void hjanr_plot_fit(int ifi, int f1, int f2, int c1, int c2, int w1, int w2, int itg) {
    	
    	double[] dumd = new double[100];
    	double[] dum  = new double[100];
    	double[] dume = new double[100];
    	double[]  obs = new double[100];
    	double[] obs1 = new double[100];
    	double[] xmin = {-5,-1.1,1.08};
    	double[] xmax = {365,1.1,1.70};

    	int indx, iplt=0, k=0;
    	double ymin=0, ymax=0, xmx=-100;
    	
    	if (w1!=w2) iplt=2;
     	if (c1!=c2) iplt=1; 
    	if (f1!=f2) iplt=0;
    	
    	Arrays.fill(j.xplt,0.0);
    	Arrays.fill(dumd,  0.0);
    	Arrays.fill(dume,  0.0);
    	    	
    	for (int iw=w1; iw<w2+1; iw++) {
    		int iww = iw*j.ndims[ifi][1]*j.ndims[ifi][0];
    		for (int ic=c1; ic<c2+1; ic++) {
    			int icc = ic*j.ndims[ifi][0];
    			for (int iff=f1; iff<f2+1; iff++) {
    				 indx = iww+icc+iff;
    				 
    				 if(j.xx[ifi][iplt][indx] < xmx) break;    				 
    				 j.xplt[k] = (float) j.xx[ifi][iplt][indx];
    				 dumd[k] =  j.crs[ifi][indx];
    				 dume[k] = j.crse[ifi][indx]; 
    				 ymin = Math.min(ymin, dumd[k]);
    				 ymax = Math.max(ymax, dumd[k]);
    				 
    				 int i1 = (int) j.xx[ifi][6][indx]-1;
    				 int i2 = (int) j.xx[ifi][5][indx]-1;
    				 
    				 hjanr_loadpar(1); //default fit
    				 double dthe = j.xx[ifi][1][indx]; //theta_cm (deg)
    				 if (dthe>1.0) dthe=Math.cos(Math.toRadians(dthe)); //convert deg to coscm  
    				 
//    				 j.jr.test = ifi==2 && dthe<-0.8 && j.xx[ifi][2][indx]>1.49 && j.xx[ifi][2][indx]<1.51 && j.xx[ifi][0][indx]>66 && j.xx[ifi][0][indx]<68;              

    				 j.jr.janrRun(j.xx[ifi][2][indx], //w
    						      j.xx[ifi][3][indx], //q2
    						      j.xx[ifi][4][indx], //eps
    						      dthe,
    						      j.xx[ifi][0][indx], //phi
    						      true);
    				 obs[k]=j.robs[i1][i2];
    				 
                     hjanr_debug(j.jr.test);    				
    				 if (j.jr.test) {
    					 System.out.println(ifi+" "+j.xx[ifi][2][indx]+" "+
    					                            j.xx[ifi][3][indx]+" "+
    					                            j.xx[ifi][4][indx]+" "+
    					                            dthe+" "+
    					                            j.xx[ifi][0][indx]+" "+
    					                            dumd[k]+" "+dume[k]+" "+(float)obs[k]); 
    					                    
    					 System.out.println(ifi+" "+indx+" "+i1+" "+i2);
    				 } 

    				 hjanr_loadpar(2); //new fit  
    				 j.jr.janrRun(j.xx[ifi][2][indx], 
						          j.xx[ifi][3][indx],
						          j.xx[ifi][4][indx],
						          dthe,
						          j.xx[ifi][0][indx],
						          true);
				     obs1[k]=j.robs[i1][i2]; 				 
	                 xmx=j.xplt[k];
	                 k++;
    			}
    		}
    	}
    	
    	if (j.xx[ifi][6][0]==1 || j.xx[ifi][6][0]==5) {
    		ymin=0.0; ymax=ymax*1.2;
    	} else {
    		ymax=Math.max(Math.abs(ymin), Math.abs(ymax));
    		ymin=-ymax*1.5; ymax=ymax*1.5;
    	}
    	
    	if (iplt==2 && j.xx[1][1][1]>1) {
    		xmin[iplt]=-5 ; xmax[iplt]=180;
    	}
    	
    	hjanr_reset_data(iplt, ifi); 
    	
    	switch (iplt) {
    	case 2: for (int i=0; i<k; i++) dgm[0].getGE("WDAT"+ ifi).addPoint(j.xplt[i],   dumd[i],  0., dume[i]);
                                        dgm[0].getGE("WSEL" +ifi).addPoint(j.xplt[itg], dumd[itg],0., dume[itg]);
                for (int i=0; i<k; i++) dgm[0].getGE("WFIT1"+ifi).addPoint(j.xplt[i],    obs[i],  0., 0.);
                for (int i=0; i<k; i++) dgm[0].getGE("WFIT2"+ifi).addPoint(j.xplt[i],   obs1[i],  0., 0.);
    	break;
    	case 1: for (int i=0; i<k; i++) dgm[0].getGE("CDAT" +ifi).addPoint(j.xplt[i],   dumd[i],  0., dume[i]);
                                        dgm[0].getGE("CSEL" +ifi).addPoint(j.xplt[itg], dumd[itg],0., dume[itg]);
                for (int i=0; i<k; i++) dgm[0].getGE("CFIT1"+ifi).addPoint(j.xplt[i],    obs[i],  0., 0.);
                for (int i=0; i<k; i++) dgm[0].getGE("CFIT2"+ifi).addPoint(j.xplt[i],    obs1[i], 0., 0.);

    	break;
    	case 0: for (int i=0; i<k; i++) dgm[0].getGE("FDAT" +ifi).addPoint(j.xplt[i],   dumd[i],  0., dume[i]);
                                        dgm[0].getGE("FSEL" +ifi).addPoint(j.xplt[itg], dumd[itg],0., dume[itg]);
                for (int i=0; i<k; i++) dgm[0].getGE("FFIT1"+ifi).addPoint(j.xplt[i],    obs[i],  0., 0.);
                for (int i=0; i<k; i++) dgm[0].getGE("FFIT2"+ifi).addPoint(j.xplt[i],    obs1[i], 0., 0.);
    	}
    }
    
    public void hjanr_plot_par() {  
    	hjanr_reset_par();
    	String                   titx =    j.pname[0]+"="+String.format("%.2f",j.xnew[0])+" ";
    	for (int i=1; i<30; i++) titx+=" "+j.pname[i]+"="+String.format("%.2f",j.xnew[i])+"    ";
    	for (int i=0; i<30; i++) dgm[1].getGE("COF1").addPoint(1+i,   j.start_value[i],0,0);   	
    	for (int icof : cofList) dgm[1].getGE("SEL1").addPoint(1+icof, j.xnew[icof],     0,0); 
    	                         dgm[1].getGE("COF1").setTitleX(titx); 
    }
        
    public void hjanr_reset_data(int iplt, int ifi) {
    	switch (iplt) {
    	case 2: dgm[0].getGE("WDAT" +ifi).reset();
                dgm[0].getGE("WSEL" +ifi).reset();
                dgm[0].getGE("WFIT1"+ifi).reset();
                dgm[0].getGE("WFIT2"+ifi).reset();
    	break;
    	case 1: dgm[0].getGE("CDAT" +ifi).reset();
                dgm[0].getGE("CSEL" +ifi).reset();                
                dgm[0].getGE("CFIT1"+ifi).reset();
                dgm[0].getGE("CFIT2"+ifi).reset();

    	break;
    	case 0: dgm[0].getGE("FDAT" +ifi).reset();
                dgm[0].getGE("FSEL" +ifi).reset();                
                dgm[0].getGE("FFIT1"+ifi).reset();
                dgm[0].getGE("FFIT2"+ifi).reset();
    	} 	        
    }
    
    public void hjanr_reset_par() {
        dgm[1].getGE("COF1").reset();
        dgm[1].getGE("SEL1").reset();  
    }
       
    public void hjanr_debug(boolean test) {

         boolean test1 = false;
         if(test1) { 
        	 for (int jj=0; jj<2; jj++) {
        	 System.out.println(
        	 (float)j.H1[1][jj]+" "+
        	 (float)j.H2[1][jj]+" "+
        	 (float)j.H3[1][jj]+" "+
        	 (float)j.H4[1][jj]+" "+
        	 (float)j.H5[1][jj]+" "+
        	 (float)j.H6[1][jj]);
        	 }
         }
         if(test1) { 
        	 for (int jj=0; jj<2; jj++) {
        	 System.out.println(
        	 (float)j.phi_amp[0][jj][1]+" "+
        	 (float)j.phi_amp[1][jj][1]+" "+
        	 (float)j.phi_amp[2][jj][1]+" "+
        	 (float)j.phi_amp[3][jj][1]+" "+
        	 (float)j.phi_amp[4][jj][1]+" "+
        	 (float)j.phi_amp[5][jj][1]);
        	 }
         }
         if(test1) { 
        	 for (int jj=0; jj<2; jj++) {
        	 System.out.println(
        	 (float)j.ReacS1[jj][1]+" "+
        	 (float)j.ReacS2[jj][1]+" "+
        	 (float)j.ReacS3[jj][1]+" "+
        	 (float)j.ReacS4[jj][1]+" "+
           	 (float)j.ReacS5[jj][1]+" "+
           	 (float)j.ReacS6[jj][1]+" "+
           	 (float)j.ReacS7[jj][1]);
        	 }
         }
         if(test1) { 
        	 for (int jj=0; jj<2; jj++) {
        	 System.out.println(
	    	 (float)j.BS1[jj][1]+" "+j.BS1[jj][4]+" "+
	    	 (float)j.BS3[1]    +" "+j.BS3[4]    +" "+
	    	 (float)j.BE1[jj][1]+" "+j.BE1[jj][4]+" "+
	    	 (float)j.BE3[1]    +" "+j.BE3[4]);  
        	 }
         }         
         if(test1) { 
        	 for (int jj=0; jj<2; jj++) {
        	 System.out.println(
        	 (float)j.BreitS1[0][jj][1]+" "+(float)j.BreitS1[1][jj][1]+" "+(float)j.BreitS1[0][jj][4]+" "+(float)j.BreitS1[1][jj][4]+" "+
        	 (float)j.BreitE1[0][jj][1]+" "+(float)j.BreitE1[1][jj][1]+" "+(float)j.BreitE1[0][jj][4]+" "+(float)j.BreitE1[1][jj][4]);
        	 System.out.println((float) j.jr.Breit[jj]+" "+(float)j.as1[jj][0]+" "+(float)j.as1[jj][1]+" "+(float)j.ae1[jj][0]+" "+(float)j.ae1[jj][1]);
        	 System.out.println(
             (float)j.BreitS3[jj][1]+" "+(float)j.BreitS3[jj][11]+" "+(float)j.BreitE3[jj][1]+" "+(float)j.BreitE3[jj][11]+" "+
        	 (float)j.BreitS3[jj][4]+" "+(float)j.BreitS3[jj][8] +" "+(float)j.BreitE3[jj][4]+" "+(float)j.BreitE3[jj][8]);        	 
        	 }
         }  	
    }

    @Override
    public void hjanr_loadpar(int n) {
    	
        double[] xload = new double[36]; 
        
        switch (n) {
        case 1: for (int i = 0; i <= 35; i++)  xload[i] = j.start_value[i]; break;
        case 2: for (int i = 0; i <= 35; i++)  xload[i] = j.xnew[i];        break;
        case 3: for (int i = 0; i <= 35; i++) j.xnew[i] = j.start_value[i]; break;
        case 4: for (int i = 0; i <= 35; i++) j.start_value[i] = j.xnew[i];
        }

        j.cm1  = xload[0];
        j.ce1  = xload[1];
        j.cs1  = xload[2];
        j.cmp1 = xload[3];
        j.csp1 = xload[4];
        j.cep2 = xload[5];
        j.csp2 = xload[6];
        j.cmp3 = xload[7];
        j.cep3 = xload[8];
        j.csp3 = xload[9];
        j.cep4 = xload[10];
        j.csp4 = xload[11];
        j.cep5 = xload[12];
        j.csp5 = xload[13];
        j.cmp6 = xload[14];
        j.cep6 = xload[15];
        j.csp6 = xload[16];
        j.cmp7 = xload[17];
        j.cep7 = xload[18];
        j.csp7 = xload[19];
        j.cmp8 = xload[20];
        j.cep8 = xload[21];
        j.csp8 = xload[22];
        j.cmp9 = xload[23];
        j.cep9 = xload[24];
        j.cmp10 = xload[25];
        j.cep10 = xload[26];
        j.cmp11 = xload[27];
        j.cep11 = xload[28];
        j.cm2 = xload[29];
        j.ce2 = xload[30];
        j.cs2 = xload[31];
        j.cm3 = xload[32];
        j.ce3 = xload[33];
        j.cs3 = xload[34];
        j.cspa = xload[35];
    }
    
    public void txt1out() {
    	
    }
    
    public void txt2out() {
    	
    }

    public static void main(String[] args) { 
    	HJanr hj = new HJanr("JANR");
    };


}
