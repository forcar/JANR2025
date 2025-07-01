package org.clas.lib;

import java.io.FileWriter;
import java.io.IOException;

public class JanrBosted2016 extends Constants {
	
    private static final double[] q2e16 = {1.72, 2.05, 2.44, 2.91, 3.48};
    private static final double[] q2eg4 = {
            0.0038, 0.0045, 0.0054, 0.0064, 0.0077,
            0.0092, 0.0110, 0.0131, 0.0156, 0.0187,
            0.0223, 0.0266, 0.0317, 0.0379, 0.0452,
            0.0540, 0.0645, 0.0770, 0.0919, 0.1100,
            0.1310, 0.1560, 0.1870, 0.2230, 0.2660,
            0.3170, 0.3790, 0.4520, 0.5400, 0.6450,
            0.7700, 0.9190, 1.1000, 1.3100, 1.5600,
            1.8700, 2.2300, 2.6600, 3.1700, 3.7900,
            4.52, 5.00
    };
    private static final double[] q2minexcl = {0.00, 0.00646, 0.0110, 0.0187, 0.0317, 0.0540, 0.0919, 0.156, 0.266, 0.452, 0.770, 1.31, 2.23, 3.79};
    private static final double[] q2maxexcl = {0.00646, 0.0110, 0.0187, 0.0317, 0.0540, 0.0919, 0.156, 0.266, 0.452, 0.770, 1.31, 2.23, 3.79, 6.45};
    private static final double[] cthmin = {-1., -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9};
    private static final double[] cthmax = {-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.7, 0.8, 0.9, 1.0};
    private static final double[] cthminc = {-0.20, 0.20, 0.44, 0.63, 0.78, 0.90};
    private static final double[] cthmaxc = {0.20, 0.44, 0.63, 0.78, 0.90, 1.0};

    public static void main(String[] args) { };
    
    public void janr_init() {
        String firun = System.getenv("QSQ");
        int iq2 = Integer.parseInt(firun.trim());
        double q2 = Math.min(4.99, 0.2 * iq2);

        System.out.printf("using q2=%.4f %02d%n", q2, iq2);
        String fname = String.format("janrpi0p%02d.dat", iq2);
        String fname2 = String.format("janrpipn%02d.dat", iq2);

        try (FileWriter fileWriter1 = new FileWriter(fname);
             FileWriter fileWriter2 = new FileWriter(fname2)) {

            janrInit();
            multAna(q2);
            readMultables();
            highMultables();

            double w = 1.5;
            double costh = 0.;
            double phid = 0.;

            for (int iq = 1; iq <= 1; iq++) {
                q2 = 0.2 * iq;
                janrIniPoint(q2);

                for (int iw = 1; iw <= 46; iw++) {
                    w = 1.1 + 0.02 * (iw - 1);

                    for (int it = 1; it <= 6; it++) {
                        costh = (cthminc[it - 1] + cthmaxc[it - 1]) / 2.;
                        phid = 90.;

                        double theta = Math.acos(costh);
                        theta = Math.toDegrees(theta);
                        double sinth = Math.sqrt(1. - costh * costh);
                        double phi = Math.toRadians(phid);

                        double s = w * w;
                        double E = s - mn22;
                        double Egamma = (E - q2) / w2;
                        double Epion = (E + mp22) / w2;
                        double q22 = Epion * Epion - mp22;

                        double k22 = (w - Egamma) * (w - Egamma) - mn22;
                        double qpion = Math.sqrt(q22);
                        double kgamma = Math.sqrt(k22);
                        double qk = qpion * kgamma * costh - Epion * Egamma;
                        double t = mp22 - q2 + 2.0 * qk;
                        double t_min = mp22 - q2 + 2.0 * (qpion * kgamma - Epion * Egamma);
                        double u = 2.0 * mn22 + mp22 - q2 - s - t;
                        double E2pion = (E + m2pion * m2pion) / w2;
                        double q2pion = E2pion * E2pion - m2pion * m2pion;
                        double Eeta = (E + meta * meta) / w2;
                        double qeta = Eeta * Eeta - meta * meta;

                        // Additional calculations and file writing can be done here
                    }
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void janrInit() {
        // Initialization code here
    }

    private static void multAna(double q2) {
        // Analysis code here
    }

    private static void readMultables() {
        // Read multables code here
    }

    private static void highMultables() {
        // High multables code here
    }

    private static void janrIniPoint(double q2) {
        // Initialize point code here
    }
}


