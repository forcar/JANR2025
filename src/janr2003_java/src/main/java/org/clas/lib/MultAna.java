package org.clas.lib;

import java.io.*;

public class MultAna extends Constants {
	
	int ierr;
   
	int lm,is,ii,jj; //common/ifintg/lm,is,ii,jj,w,q2
	double w,q2;     //common/ifintg/lm,is,ii,jj,w,q2
	
    double mp = 0.938;
	
    JanrRePhi jrp = new JanrRePhi(); 

	public MultAna() {

	}

    public void multAna(double qq) {
    	
        if (IntAccur == 0) IntAccur = 0.05f; 
    	
        ierr = 0;        
       
        int Nw = 101;
        double W_min = 1.1;
        double W_step = 0.01;

        for (jj = 0; jj < 1; jj++) {
        try {
            BufferedWriter[] luns = new BufferedWriter[100]; // Based on max lun=(lm+2)*10+is; lm=0..3, is=1..6 => lun=20..46
            int n=0;
            for (lm = 0; lm < 4; lm++) {
                for (is = 0; is < 6; is++) { 
                    int lun = (lm+2)*10+is;
                    int ifile = lm*6+is; 
                    String filename = janrPath+"multipols/" + fnamr[n]; 
                    luns[lun] = new BufferedWriter(new FileWriter(filename, false));
                    luns[lun].write(fnamr[ifile] + "\n"); n++;
                    luns[lun].write("\n");
                }
            }
            
            q2 = qq;

            for (int iw = 1; iw <= Nw; iw++) {
                w = W_min + (iw-1)*W_step;
                double El = (w*w-mp*mp)/2.0/mp;

                for (lm = 0; lm < 4; lm++) {
                    for (is = 0; is < 3; is++) {
                        for (ii = 0; ii < 6; ii++) {
                            Rint[lm][is][ii][jj] = gauss("Fintg", -1.0, 1.0, IntAccur);
                        }
                    }
                    for (is = 0; is < 3; is++) {
                        re_mult();
                    }
                    
                    int lun = (lm+2)*10;
                    
                    write201(luns[lun+0], w,El,Mlp,lm,jj);
                    write201(luns[lun+1], w,El,Mlm,lm,jj);
                    write201(luns[lun+2], w,El,Elp,lm,jj);
                    write201(luns[lun+3], w,El,Elm,lm,jj);
                    write201(luns[lun+4], w,El,Slp,lm,jj);
                    write201(luns[lun+5], w,El,Slm,lm,jj);
                }
            }

            for (int lm = 0; lm < 3; lm++) {
                for (int is = 0; is < 6; is++) {
                    int lun = (lm+2)*10+is;
                    if (luns[lun] != null) luns[lun].close();
                }
            }

        } catch (IOException e) {
            System.err.println("Error writing jj=0 data: " + e.getMessage());
        }
        }
        
        for (jj = 1; jj < 2; jj++) {
        try {
            BufferedWriter[] luns = new BufferedWriter[100];
            int n=0;
            for (lm = 0; lm < 4; lm++) {
                for (is = 0; is < 6; is++) {
                    int lun = (lm+2)*10+is;
                    int ifile = lm*6+is;
                    String filename = janrPath+"multipols/" + fnamrr[n];
                    luns[lun] = new BufferedWriter(new FileWriter(filename, false));
                    luns[lun].write(fnamrr[ifile] + "\n"); n++;
                    luns[lun].write("\n");
                }
            }

            q2 = qq;

            for (int iw = 1; iw <= Nw; iw++) {
                w = W_min + (iw-1)*W_step;
                double El = (w*w-mp*mp)/2.0/mp;

                for (lm = 0; lm < 4; lm++) {
                    for (is = 0; is < 3; is++) {
                        for (ii = 0; ii < 6; ii++) {
                            Rint[lm][is][ii][jj] = gauss("Fintg", -1.0, 1.0, IntAccur);
                        }
                    }
                    for (is = 0; is < 3; is++) {
                        re_mult();
                    }
                    
                    int lun = (lm+2)*10;
                    
                    write201(luns[lun+0], w,El,Mlp,lm,jj);
                    write201(luns[lun+1], w,El,Mlm,lm,jj);
                    write201(luns[lun+2], w,El,Elp,lm,jj);
                    write201(luns[lun+3], w,El,Elm,lm,jj);
                    write201(luns[lun+4], w,El,Slp,lm,jj);
                    write201(luns[lun+5], w,El,Slm,lm,jj);
                }
            }

            for (int lm = 0; lm < 3; lm++) {
                for (int is = 0; is < 6; is++) {
                    int lun = (lm+2)*10+is;
                    if (luns[lun] != null) luns[lun].close();
                }
            }
        } catch (IOException e) {
            System.err.println("Error writing jj=1 data: " + e.getMessage());
        }
        }

        q2 = qq;
        jj = 1;
        
        try {
            BufferedWriter[] luns = new BufferedWriter[100]; 
            int n=0;
            for (lm = 4; lm < 6; lm++) {
                for (is = 0; is < 6; is++) { 
                    int lun = (lm-4)*10+is;
                    int ifile = (lm-4)*6+is; 
                    String filename = janrPath+"multipols/" + fnamhigh[n];
                    luns[lun] = new BufferedWriter(new FileWriter(filename, false));
                    luns[lun].write(fnamhigh[ifile] + "\n"); n++;
                    luns[lun].write("\n");
                }
            }

            for (int iw = 1; iw <= Nw; iw++) {
                w = W_min + (iw-1)*W_step;
                double El = (w*w-mp*mp)/2.0/mp;

                for (lm = 4; lm < 6; lm++) {
                    for (is = 0; is < 3; is++) {
                        for (ii = 0; ii < 6; ii++) {
                            Rint[lm][is][ii][jj] = gauss("Finthigh", -1.0, 1.0, IntAccur);
                        }
                    }
                    for (is = 0; is < 3; is++) {
                        re_mult();
                    }
                    
                    int lun = (lm-4)*10;
                    
                    write201(luns[lun+0], w,El,Mlp,lm,jj);
                    write201(luns[lun+1], w,El,Mlm,lm,jj);
                    write201(luns[lun+2], w,El,Elp,lm,jj);
                    write201(luns[lun+3], w,El,Elm,lm,jj);
                    write201(luns[lun+4], w,El,Slp,lm,jj);
                    write201(luns[lun+5], w,El,Slm,lm,jj);
                }
            }

            for (int lm = 4; lm < 6; lm++) {
                for (int is = 0; is < 6; is++) {
                    int lun = (lm-4)*10+is;
                    if (luns[lun] != null) luns[lun].close();
                }
            }

        } catch (IOException e) {
            System.err.println("Error writing jj=0 data: " + e.getMessage());
        }
        }

    private static void write201(BufferedWriter out, double W, double El, double[][][] arr, int lm, int jj) throws IOException {
        String line = String.format("  %7.3f %7.3f  %12.5f %12.5f %12.5f\n", W, El, arr[lm][0][jj], arr[lm][1][jj], arr[lm][2][jj]);
        out.write(line);
        out.flush();
    }
    
    interface Fint {
        double apply(double x);
    }
    
    class ffImpl implements Fint {
    	@Override
    	public double apply(double x) {
    		return ff(x);
    	}
    }
    
    class FintgImpl implements Fint {
    	@Override
    	public double apply(double x) {
    		return Fintg(x);
    	}
    }
    
    class FinthighImpl implements Fint {
    	@Override
    	public double apply(double x) {
    		return Finthigh(x);
    	}
    }
    
    static double ff(double x) { //for testing
        return Math.exp(- x * x / 2) / Math.sqrt(2 * Math.PI);
    }
    
    public float gauss(String fun, double a, double b, double eps) {
    	Fint f=null;
    	if(fun=="ff")       f = new ffImpl();
    	if(fun=="Fintg")    f = new FintgImpl();
    	if(fun=="Finthigh") f = new FinthighImpl();
        return dgauss(f,a,b,eps);
    }
      
  //https://github.com/apc-llc/cernlib/blob/ef83476ab09eed4972aff95685fcb27c24e6a18c/2005/src/packlib/kernlib/kernnum/obsolete/d103fort/gauss.F#L10 (FORTRAN)
  //https://codingfleet.com/code-converter/java/ (DeepSeek V3) (JAVA)
    
    private static final double[] W = {
            0.10122853629037625915, 0.22238103445337470435,
            0.31370664587788728733, 0.36268378337836198296,
            0.02715245941175409485, 0.06225352393864789286,
            0.09515851168249278480, 0.12462897125553387205,
            0.14959598881657673208, 0.16915651939500253818,
            0.18260341504492358886, 0.18945061045505496849
        };

    private static final double[] X = {
            0.96028985649753623168, 0.79666647741362673959,
            0.52553240991632898581, 0.18343464249564980493,
            0.98940093499164993259, 0.94457502307323257607,
            0.86563120238783174388, 0.75540440835500303389,
            0.61787624440264374844, 0.45801677765722738634,
            0.28160355077925891323, 0.09501250983763744018
    };

    public static float dgauss(Fint f, double a, double b, double eps) {
    	
        if (b == a) {
            return 0.0f;
        }

        double gauss = 0.0;
        double s8,s16;
        final double constVal = 0.005 / (b - a);
        double bb = a;

        while (true) {
            double aa = bb;
            bb = b;

            while (true) {
                double c1 = 0.5 * (bb + aa);
                double c2 = 0.5 * (bb - aa);

                s8 = 0.0;
                for (int i = 0; i < 4; i++) {
                    double u = c2 * X[i];
                    s8 += W[i] * (f.apply(c1 + u) + f.apply(c1 - u));
                }
                s8 *= c2;

                s16 = 0.0;
                for (int i = 4; i < 12; i++) {
                    double u = c2 * X[i];
                    s16 += W[i] * (f.apply(c1 + u) + f.apply(c1 - u));
                }
                s16 *= c2;

                if (Math.abs(s16 - s8) <= eps * (1.0 + Math.abs(s16))) {
                    break;
                }

                bb = c1;
                if (1.0 + Math.abs(constVal * c2) == 1.0) {
                    System.err.println("FUNCTION DGAUSS ... TOO HIGH ACCURACY REQUIRED");
                    return 0.0f;
                }
            }

            gauss += s16;
            if (bb == b) {
                break;
            }
        }

        return (float) gauss;
    }    
    
    public double Fintg(double x) {
    	
    	int ierr=0;
    	
    	double fintg;

        int MAX_PHI_INDEX = 6;  
        int MAX_IS_INDEX  = 3;  
        int MAX_J_INDEX   = 2;  
        
        double[][][] PHis = new double[MAX_PHI_INDEX+1][MAX_IS_INDEX+1][MAX_J_INDEX+1];
        
        jrp.janrRePhi(w, q2, x, ierr);
            
        if (is == 2) {
            for (int i = 0; i < 6; i++) {
                PHis[i][is][jj] = PHisot[jj][i][2];
            }
        } else if (is == 1) {
            for (int i = 0; i < 6; i++) {
                PHis[i][is][jj] = PHisot[jj][i][0] + 2.0 * PHisot[jj][i][1];
            }
        } else if (is == 0) {
            for (int i = 0; i < 6; i++) {
                PHis[i][is][jj] = PHisot[jj][i][0] - PHisot[jj][i][1];
            }
        }
            
        double p = 0.0, p1 = 0.0, p2 = 0.0;
            
        switch (lm) {
            case 0:
                p = 1.0;
                p1 = 0.0;
                p2 = 0.0;
                break;
            case 1:
                p = x;
                p1 = -Math.sqrt(1.0 - x * x);
                p2 = 0.0;
                break;
            case 2:
                p = 3.0 * x * x / 2.0 - 0.5;
                p1 = -3.0 * x * Math.sqrt(1.0 - x * x);
                p2 = 3.0 * (1.0 - x * x);
                break;
            case 3:
                p = (5.0 * x * x - 3.0) * x / 2.0;
                p1 = -Math.sqrt(1.0 - x * x) * (15.0 * x * x - 3.0) / 2.0;
                p2 = 15.0 * x * (1.0 - x * x);
                break;
            default:
                break;
            }
            
        fintg = 0.0;
            
        switch (ii) {
            case 0:
                if (lm < 1) return fintg;
                double rim = p1 * Math.sqrt(1.0 - x * x) * (PHis[2][is][jj] + x * PHis[3][is][jj]);
                fintg = rim * (2.0 * lm + 1.0) / (2.0 * lm * (lm + 1.0));
                break;
            case 1:
                rim = p * (2.0 * PHis[0][is][jj] - 2.0 * x * PHis[1][is][jj] + (1.0 - x * x) * PHis[3][is][jj]);
                fintg = rim * (2.0 * lm + 1.0) / 2.0;
                break;
            case 2:
                if (lm < 1) return fintg;
                rim = p1 * Math.sqrt(1.0 - x * x) * PHis[1][is][jj];
                fintg = -rim * (2.0 * lm + 1.0) / (2.0 * lm * (lm + 1.0));
                break;
            case 3:
                if (lm < 2) return fintg;
                rim = p2 * (1.0 - x * x) * PHis[3][is][jj];
                fintg = rim * (2.0 * lm + 1.0) / (2.0 * lm * (lm + 1.0) * (lm - 1) * (lm + 2));
                break;
            case 4:
                rim = p * (PHis[4][is][jj] + x * PHis[5][is][jj]);
                fintg = rim * (2.0 * lm + 1.0) / 2.0;
                break;
            case 5:
                if (lm < 1) return fintg;
                rim = p1 * Math.sqrt(1.0 - x * x) * PHis[5][is][jj];
                fintg = -rim * (2.0 * lm + 1.0) / (2.0 * lm * (lm + 1.0));
                break;
            default:
                break;
            }
            
            return fintg;          	
    }

    public double Finthigh(double x) {
    	
        double p0, p1, p2, p3, p4, p5, p6;
        double[][]      rim = new double[6][3]; 
        double[][]     PHis = new double[6][3]; 
        double[][]      pim = new double[6][3]; 
        double[][][]    pip = new double[6][6][3]; // was 6,5,2 in Inna's code, 
        double[][] PHisotop = new double[6][3]; 

        double extB, extR;
        int ierr = 0;
        
        jrp.janrRePhi(w, q2, x, ierr);

        extB = 1.0 / (1.0 + Math.pow(w * w - 1.16, 2));
        extR = extB * Math.pow(w * w - 1.16, 2);

        for (int i = 0; i < 6; i++) {
            for (int ir = 0; ir < 3; ir++) {
                PHisotop[i][ir] = PHisot[0][i][ir] * extB + PHisot[1][i][ir] * extR;
            }
        }

        // is=0 for (3/2), is=1 for (1/2), is=2 for (0)

        if (is == 2) {
            // (0)
            for (int i = 0; i < 6; i++) {
                PHis[i][is] = PHisotop[i][2];
            }
        } else if (is == 1) {
            // (1/2)
            for (int i = 0; i < 6; i++) {
                PHis[i][is] = PHisotop[i][0] + 2.0 * PHisotop[i][1];
            }
        } else if (is == 0) {
            // (3/2)
            for (int i = 0; i < 6; i++) {
                PHis[i][is] = PHisotop[i][0] - PHisotop[i][1];
            }
        }

        p0 = 1.0 - x * x;
        // P_4
        p1 = (3.0 - 30.0 * x * x + 35.0 * Math.pow(x, 4)) / 8.0;
        // P_5
        p2 = (15.0 * x - 70.0 * Math.pow(x, 3) + 63.0 * Math.pow(x, 5)) / 8.0;
        // P'_4
        p3 = -7.5 * x + 17.5 * Math.pow(x, 3);
        // P'_5
        p4 = 1.875 - 26.25 * x * x + 39.375 * Math.pow(x, 4);
        // P''_4
        p5 = -7.5 + 52.5 * x * x;
        // P''_5
        p6 = -52.5 * x + 157.5 * Math.pow(x, 3);

        if (is >= 0 && is <= 2) {
            rim[0][is] = PHis[2][is] + x * PHis[3][is];
            rim[1][is] = 2.0 * PHis[0][is] - 2.0 * x * PHis[1][is] + p0 * PHis[3][is];
            rim[2][is] = PHis[1][is];
            rim[3][is] = PHis[3][is];
            rim[4][is] = PHis[4][is] + x * PHis[5][is];
            rim[5][is] = PHis[5][is];

            pim[0][0] = -9.0 * p0 * p3 / 45.0;
            pim[0][1] = -11.0 * p0 * p4 / 60.0;
            pim[1][0] = 4.5 * p1;
            pim[1][1] = 5.5 * p2;
            pim[2][0] = -pim[0][0];
            pim[2][1] = -pim[0][1];
            pim[3][0] = p0 * p0 * p5 / 80.0;
            pim[3][1] = 11.0 * p0 * p0 * p6 / 1680.0;
            pim[4][0] = pim[1][0];
            pim[4][1] = pim[1][1];
            pim[5][0] = pim[2][0];
            pim[5][1] = pim[2][1];

            if (lm == 5) {
                for (int i = 0; i < 6; i++) {
                    pip[i][lm][is] = pim[i][1] * rim[i][is];
                }
            } else if (lm == 4) {
                for (int i = 0; i < 6; i++) {
                    pip[i][lm][is] = pim[i][0] * rim[i][is];
                }
            }
        }

        double finthigh = 0.0;

        switch (ii) {
            case 0:
                // Integral 1 (I1)
                finthigh = pip[0][lm][is];
                break;
            case 1:
                // Integral 2 (I2)
                finthigh = pip[1][lm][is];
                break;
            case 2:
                // Integral 3 (I3)
                finthigh = pip[2][lm][is];
                break;
            case 3:
                // Integral 4 (I4)
                finthigh = pip[3][lm][is];
                break;
            case 4:
                // Integral 5 (I5)
                finthigh = pip[4][lm][is];
                break;
            case 5:
                // Integral 6 (I6)
                finthigh = pip[5][lm][is];
                break;
            default:
                // Invalid i value
                throw new IllegalArgumentException("Parameter ii out of range (1..6), got i=" + ii);
        }

        return finthigh;	
    }
    
    public void re_mult() {
    	
        int lmIdx = lm;
        int isIdx = is;
        int jjIdx = jj;

        if (lm >= 0) {
            Elp[lmIdx][isIdx][jjIdx] = (float) (
                (-2 * lm * Rint[lmIdx][isIdx][0][jjIdx]
                         + Rint[lmIdx][isIdx][1][jjIdx]
                 +2 * lm * Rint[lmIdx][isIdx][2][jjIdx]
         - lm * (lm - 1) * Rint[lmIdx][isIdx][3][jjIdx])
                         /(2 * (lm + 1) * (2 * lm + 1))
            );
            Slp[lmIdx][isIdx][jjIdx] = (float) (
                          (Rint[lmIdx][isIdx][4][jjIdx]
                    - lm * Rint[lmIdx][isIdx][5][jjIdx])
                         /     (lm + 1) / (2 * lm + 1)
            );
            Mlp[lmIdx][isIdx][jjIdx] = 0.0f;
            Mlm[lmIdx][isIdx][jjIdx] = 0.0f;
            Elm[lmIdx][isIdx][jjIdx] = 0.0f;
            Slm[lmIdx][isIdx][jjIdx] = 0.0f;
        }

        if (lm >= 1) {
            Mlp[lmIdx][isIdx][jjIdx] = (float) (
                  (2 *               Rint[lmIdx][isIdx][0][jjIdx]
                                   + Rint[lmIdx][isIdx][1][jjIdx]
                      + 2 * lm     * Rint[lmIdx][isIdx][2][jjIdx]
           + (lm - 1) * (lm + 2)   * Rint[lmIdx][isIdx][3][jjIdx])
                                   / (2  * (lm + 1) * (2 * lm + 1))
            );
            Mlm[lmIdx][isIdx][jjIdx] = (float) (
                  (-2 *              Rint[lmIdx][isIdx][0][jjIdx]
                                   - Rint[lmIdx][isIdx][1][jjIdx]
                  + 2 * (lm + 1)   * Rint[lmIdx][isIdx][2][jjIdx]
           - (lm - 1) * (lm + 2)   * Rint[lmIdx][isIdx][3][jjIdx])
                                   / (2 * lm * (2 * lm + 1))
            );
            Slm[lmIdx][isIdx][jjIdx] = (float) (
                                    (Rint[lmIdx][isIdx][4][jjIdx]
                 +      (lm + 1)   * Rint[lmIdx][isIdx][5][jjIdx])
                                   / (lm * (2 * lm + 1))
            );
        }

        if (lm >= 2) {
            Elm[lmIdx][isIdx][jjIdx] = (float) (
                   (2 * (lm + 1)   * Rint[lmIdx][isIdx][0][jjIdx]
                                   + Rint[lmIdx][isIdx][1][jjIdx]
                  - 2 * (lm + 1)   * Rint[lmIdx][isIdx][2][jjIdx]
           - (lm + 1) * (lm + 2)   * Rint[lmIdx][isIdx][3][jjIdx])
                                   / (2 * lm * (2 * lm + 1))
            );
        }
    }

    public static void main(String[] args) {
    	MultAna ma = new MultAna();
        double result = ma.gauss("ff", 0, 1, 1e-6);
        System.out.println("Integration result: " + result); 
    }
}
