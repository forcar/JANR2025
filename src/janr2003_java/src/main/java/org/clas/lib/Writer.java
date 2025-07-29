package org.clas.lib;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;

public class Writer extends Constants {
	
	public Writer() {
		
	}
	
    static final int NW = 60;
    static final double[] MPION = {0.13498, 0.13957};   // pion masses
    static final double MP = 0.93827; 
	
    MultAna ma = new MultAna();
	JanrRun jr = new JanrRun();
	
	public void MultAna(double q2) {
        System.out.println("Janr03.Writer.MultAna("+(float)q2+")");
		ma.multAna(q2);
	}
	
	public void generateTables() {
		
        double q2 = 0.5 * (q2min + q2max);
        double mpa = MPION[0] * MPION[0];
        double mpb = MPION[1] * MPION[1];
        double mp2 = MP * MP;
        double eps = 1.0;
        
        System.out.printf("Janr03.Writer.generateTables: Tables at Q2=%.6f  Ebeam=%.6f%n", q2, ebeam);
	
        try (
                PrintWriter outPi0MF  = new PrintWriter(new BufferedWriter(new FileWriter(janrPath+"mf1g-jan-pi0.tbl")));
                PrintWriter outPiPMF  = new PrintWriter(new BufferedWriter(new FileWriter(janrPath+"mf1g-jan-pip.tbl")));
                PrintWriter outPi0RSP = new PrintWriter(new BufferedWriter(new FileWriter(janrPath+"rspf-jan-pi0.tbl")));
                PrintWriter outPiPRSP = new PrintWriter(new BufferedWriter(new FileWriter(janrPath+"rspf-jan-pip.tbl")))
            ) {
                for (int iw = 1; iw <= NW; iw++) {
                    double w = 1.1 + (iw - 1) * 0.01;

                    if (ebeam > 0) {
                        double nu = (w*w + q2 - mp2) / (2.0 * MP);
                        double ep = ebeam - nu;
                        double denom = 4.0 * ebeam * ep - q2;
                        eps = 1.0 / (1.0 + 2.0 * (q2 + nu*nu) / denom);
                    }

                    if (eps > 0) {
                        double termA = (w*w - mp2 + mpa) / (2.0 * w);
                        double qpi1 = Math.sqrt(termA*termA - mpa);
                        double termB = (w*w - mp2 + mpb) / (2.0 * w);
                        double qpi2 = Math.sqrt(termB*termB - mpb);
                        double ratio = qpi2 / qpi1;  // not used further here

                        for (int ic = 1; ic <= 21; ic++) {
                            double coscm = -1.0 + (ic - 1) * 0.1;
                            for (int ip = 1; ip <= 24; ip++) {
                                double phi = 7.5 + (ip - 1) * 15.0;

                                jr.janrRun(w, q2, eps, coscm, phi, true);

                                if ((iw % 2) == 1) {

                                    outPi0MF .printf("%10.3f%10.3f%10.3f%10.3f%10.3f%15.5e%15.5e%15.5e%15.5e%15.5e%n",
                                                     w, q2, eps, coscm, phi,
                                                     sigma[0], sigma[0]/10.0,
                                                     0.0, 0.0, 0.0);
                                    outPiPMF .printf("%10.3f%10.3f%10.3f%10.3f%10.3f%15.5e%15.5e%15.5e%15.5e%15.5e%n",
                                                     w, q2, eps, coscm, phi,
                                                     sigma[1], sigma[1]/10.0,
                                                     0.0, 0.0, 0.0);
                                }
                            }

                            outPi0RSP.printf("%10.3f%10.3f%10.3f%10.3f%15.5e%15.5e%15.5e%n",
                                             w, q2, eps, coscm,
                                             robs[6][0], robs[7][0], robs[8][0]);
                            outPiPRSP.printf("%10.3f%10.3f%10.3f%10.3f%15.5e%15.5e%15.5e%n",
                                             w, q2, eps, coscm,
                                             robs[6][1], robs[7][1], robs[8][1]);
                        }
                    }
                }

            }
            catch (IOException e) {
                e.printStackTrace();
            }	
	}

}
