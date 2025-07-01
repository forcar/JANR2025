package org.clas.lib;

public class Fintg extends Constants {
	
    // Shared variables from common block
    private int lm;
    private int is;
    private int i;
    private int j;
    private double w;
    private double q2;
    
    // Constants and helper variables
    private static final int MAX_PHI_INDEX = 6;
    private static final int MAX_IS_INDEX = 3;
    private static final int MAX_J_INDEX = 2;
    
    JanrRePhi jrp = new JanrRePhi();
    
    public Fintg(double w, double q2, int lm, int is, int i, int j) {
        this.w = w;
        this.q2 = q2;
        this.lm = lm;
        this.is = is;
        this.i = i;
        this.j = j;
    }
    
    public double fintg(double x) {
        double[][][] PHis = new double[MAX_PHI_INDEX+1][MAX_IS_INDEX+1][MAX_J_INDEX+1];
        
        int ierr = 0;
        jrp.janrRePhi(w, q2, x, ierr);
        
        // Calculate PHis based on is value
        if (is == 3) {
            for (int ii = 1; ii <= 6; ii++) {
                PHis[ii][is][j] = PHisot[j][ii][3];
            }
        } else if (is == 2) {
            for (int ii = 1; ii <= 6; ii++) {
                PHis[ii][is][j] = PHisot[j][ii][1] + 2.0 * PHisot[j][ii][2];
            }
        } else if (is == 1) {
            for (int ii = 1; ii <= 6; ii++) {
                PHis[ii][is][j] = PHisot[j][ii][1] - PHisot[j][ii][2];
            }
        }
        
        double p = 0.0, p1 = 0.0, p2 = 0.0;
        
        // Switch based on lm value (Fortran uses 1-based indexing)
        switch (lm + 1) {
            case 1:
                p = 1.0;
                p1 = 0.0;
                p2 = 0.0;
                break;
            case 2:
                p = x;
                p1 = -Math.sqrt(1.0 - x * x);
                p2 = 0.0;
                break;
            case 3:
                p = 3.0 * x * x / 2.0 - 0.5;
                p1 = -3.0 * x * Math.sqrt(1.0 - x * x);
                p2 = 3.0 * (1.0 - x * x);
                break;
            case 4:
                p = (5.0 * x * x - 3.0) * x / 2.0;
                p1 = -Math.sqrt(1.0 - x * x) * (15.0 * x * x - 3.0) / 2.0;
                p2 = 15.0 * x * (1.0 - x * x);
                break;
            default:
                // Handle unexpected lm values
                break;
        }
        
        double fintg = 0.0;
        
        // Switch based on i value (Fortran uses 1-based indexing)
        switch (i) {
            case 1:
                // Integral 1 (I1)
                if (lm < 1) return fintg;
                double rim = p1 * Math.sqrt(1.0 - x * x) * (PHis[3][is][j] + x * PHis[4][is][j]);
                fintg = rim * (2.0 * lm + 1.0) / (2.0 * lm * (lm + 1.0));
                break;
            case 2:
                // Integral 2 (I2)
                rim = p * (2.0 * PHis[1][is][j] - 2.0 * x * PHis[2][is][j] + (1.0 - x * x) * PHis[4][is][j]);
                fintg = rim * (2.0 * lm + 1.0) / 2.0;
                break;
            case 3:
                // Integral 3 (I3)
                if (lm < 1) return fintg;
                rim = p1 * Math.sqrt(1.0 - x * x) * PHis[2][is][j];
                fintg = -rim * (2.0 * lm + 1.0) / (2.0 * lm * (lm + 1.0));
                break;
            case 4:
                // Integral 4 (I4)
                if (lm < 2) return fintg;
                rim = p2 * (1.0 - x * x) * PHis[4][is][j];
                fintg = rim * (2.0 * lm + 1.0) / (2.0 * lm * (lm + 1.0) * (lm - 1) * (lm + 2));
                break;
            case 5:
                // Integral 5 (I5)
                rim = p * (PHis[5][is][j] + x * PHis[6][is][j]);
                fintg = rim * (2.0 * lm + 1.0) / 2.0;
                break;
            case 6:
                // Integral 6 (I6)
                if (lm < 1) return fintg;
                rim = p1 * Math.sqrt(1.0 - x * x) * PHis[6][is][j];
                fintg = -rim * (2.0 * lm + 1.0) / (2.0 * lm * (lm + 1.0));
                break;
            default:
                // Handle unexpected i values
                break;
        }
        
        return fintg;
    }
    
    // Getters and setters for the shared variables
    public int getLm() { return lm; }
    public void setLm(int lm) { this.lm = lm; }
    
    public int getIs() { return is; }
    public void setIs(int is) { this.is = is; }
    
    public int getI() { return i; }
    public void setI(int i) { this.i = i; }
    
    public int getJ() { return j; }
    public void setJ(int j) { this.j = j; }
    
    public double getW() { return w; }
    public void setW(double w) { this.w = w; }
    
    public double getQ2() { return q2; }
    public void setQ2(double q2) { this.q2 = q2; }
}
