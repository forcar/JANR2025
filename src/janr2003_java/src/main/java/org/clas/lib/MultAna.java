package org.clas.lib;

import java.io.*;

public class MultAna extends Constants {
	
	int lm,is,i,j; //common/ifintg/lm,is,i,j,w,q2
	double w,q2;   //common/ifintg/lm,is,i,j,w,q2
	
    JanrRePhi jrp = new JanrRePhi(); 

	public MultAna() {

	}

    public void multAna(double qq) {
    	
        int ierr = 0;        

        String[] fnamr = {
            "rm0p.dat","rm0m.dat","re0p.dat","re0m.dat","rs0p.dat","rs0m.dat",
            "rm1p.dat","rm1m.dat","re1p.dat","re1m.dat","rs1p.dat","rs1m.dat",
            "rm2p.dat","rm2m.dat","re2p.dat","re2m.dat","rs2p.dat","rs2m.dat",
            "rm3p.dat","rm3m.dat","re3p.dat","re3m.dat","rs3p.dat","rs3m.dat"
        };
        String[] fnamrr = {
            "rrm0p.dat","rrm0m.dat","rre0p.dat","rre0m.dat",
            "rrs0p.dat","rrs0m.dat",
            "rrm1p.dat","rrm1m.dat","rre1p.dat","rre1m.dat","rrs1p.dat","rrs1m.dat",
            "rrm2p.dat","rrm2m.dat","rre2p.dat","rre2m.dat","rrs2p.dat","rrs2m.dat",
            "rrm3p.dat","rrm3m.dat","rre3p.dat","rre3m.dat","rrs3p.dat","rrs3m.dat"
        };
        String[] fnamhigh = {
            "rm4p.dat","rm4m.dat","re4p.dat","re4m.dat","rs4p.dat","rs4m.dat",
            "rm5p.dat","rm5m.dat","re5p.dat","re5m.dat","rs5p.dat","rs5m.dat"
        };

        double mp = 0.938;

        cm1  = start_value[0];
        ce1  = start_value[1];
        cs1  = start_value[2];
        cmp1 = start_value[3];
        csp1 = start_value[4];
        cep2 = start_value[5];
        csp2 = start_value[6];
        cmp3 = start_value[7];
        cep3 = start_value[8];
        csp3 = start_value[9];
        cep4 = start_value[10];
        csp4 = start_value[11];
        cep5 = start_value[12];
        csp5 = start_value[13];
        cmp6 = start_value[14];
        cep6 = start_value[15];
        csp6 = start_value[16];
        cmp7 = start_value[17];
        cep7 = start_value[18];
        csp7 = start_value[19];
        cmp8 = start_value[20];
        cep8 = start_value[21];
        csp8 = start_value[22];
        cmp9 = start_value[23];
        cep9 = start_value[24];
        cmp10 = start_value[25];
        cep10 = start_value[26];
        cmp11 = start_value[27];
        cep11 = start_value[28];
        cspa  = start_value[29];


        if (IntAccur == 0) IntAccur = 5; // Fortran: 0.05; but IntAccur is int in Java, so using 5 (0.05*100). Adjust as needed.

        System.out.println("MULT_ANA: IntAccur = " + IntAccur);

        System.out.print(String.format("%25s", "Value of parameters => 1-30 "));
        
        System.out.println();
        for (int i = 0; i < parms; i++) {
            System.out.print(String.format("%9.5f", start_value[i]));
            if ((i+1)%6 == 0) System.out.println();
        }
        System.out.println();

        System.out.println("MULT_ANA: Evaluating background at Q2= " + qq);

        int Nw = 101;
        double W_min = 1.1;
        double W_step = 0.01;

        System.out.println("MULT_ANA: Number of W points = " + Nw);
        System.out.println("MULT_ANA: W_min,W_step = " + W_min + "," + W_step);

        for (int jj = 1; jj <= 1; jj++) {
        try {
            BufferedWriter[] luns = new BufferedWriter[40]; // Based on max lun=(lm+2)*10+is; lm=0..3, is=1..6 => lun=20..46
    
            for (int lm = 0; lm <= 3; lm++) {
                for (int is = 1; is <= 6; is++) {
                    int lun = (lm+2)*10+is;
                    int ifile = lm*6+(is-1); 
                    String filename = janrPath+"multipols/" + fnamr[ifile];
                    luns[lun] = new BufferedWriter(new FileWriter(filename, false));
                    luns[lun].write(fnamr[ifile] + "\n");
                    luns[lun].write("\n");
                }
            }
            

            double q2 = qq;

            for (int iw = 1; iw <= Nw; iw++) {
                double W = W_min + (iw-1)*W_step;
                double El = (Math.pow(W,2)-Math.pow(mp,2))/2.0/mp;

                for (int lm = 0; lm <= 3; lm++) {
                    for (int is = 1; is <= 3; is++) {
                        for (int ii = 1; ii <= 6; ii++) {
                            Rint[lm][is][ii][jj] = gauss("Fintg", -1.0, 1.0, IntAccur);
                        }
                    }
                    for (int is = 1; is <= 3; is++) {
                        re_mult();
                    }
                    int lun = (lm+2)*10;
                    write201(luns[lun+1], W,El,Mlp,lm,jj);
                    write201(luns[lun+2], W,El,Mlm,lm,jj);
                    write201(luns[lun+3], W,El,Elp,lm,jj);
                    write201(luns[lun+4], W,El,Elm,lm,jj);
                    write201(luns[lun+5], W,El,Slp,lm,jj);
                    write201(luns[lun+6], W,El,Slm,lm,jj);
                }
            }

            for (int lm = 0; lm <= 3; lm++) {
                for (int is = 1; is <= 6; is++) {
                    int lun = (lm+2)*10+is;
                    if (luns[lun] != null) luns[lun].close();
                }
            }

        } catch (IOException e) {
            System.err.println("Error writing jj=1 data: " + e.getMessage());
        }
        }
        
        for (int  jj = 2; jj <= 2; jj++) {
        try {
            BufferedWriter[] luns = new BufferedWriter[40];

            for (int lm = 0; lm <= 3; lm++) {
                for (int is = 1; is <= 6; is++) {
                    int lun = (lm+2)*10+is;
                    int ifile = lm*6+(is-1);
                    String filename = "multipols/" + fnamrr[ifile];
                    luns[lun] = new BufferedWriter(new FileWriter(filename, false));
                    luns[lun].write(fnamrr[ifile] + "\n");
                    luns[lun].write("\n");
                }
            }

            double Q2 = qq;

            for (int iw = 1; iw <= Nw; iw++) {
                double W = W_min + (iw-1)*W_step;
                double El = (Math.pow(W,2)-Math.pow(mp,2))/2.0/mp;

                for (int lm = 0; lm <= 3; lm++) {
                    for (int is = 1; is <= 3; is++) {
                        for (int ii = 1; ii <= 6; ii++) {
                            Rint[lm][is][ii][jj] = gauss("Fintg", -1.0, 1.0, IntAccur);
                        }
                    }
                    for (int is = 1; is <= 3; is++) {
                        re_mult();
                    }
                    int lun = (lm+2)*10;
                    write201(luns[lun+1], W,El,Mlp,lm,jj);
                    write201(luns[lun+2], W,El,Mlm,lm,jj);
                    write201(luns[lun+3], W,El,Elp,lm,jj);
                    write201(luns[lun+4], W,El,Elm,lm,jj);
                    write201(luns[lun+5], W,El,Slp,lm,jj);
                    write201(luns[lun+6], W,El,Slm,lm,jj);
                }
            }

            for (int lm = 0; lm <= 3; lm++) {
                for (int is = 1; is <= 6; is++) {
                    int lun = (lm+2)*10+is;
                    if (luns[lun] != null) luns[lun].close();
                }
            }
        } catch (IOException e) {
            System.err.println("Error writing jj=2 data: " + e.getMessage());
        }
        }

        try {
        for (int lm = 4; lm <= 5; lm++) {
            for (int is = 1; is <= 6; is++) {
                int lun = (lm+2)*10+is;
                int ifile = (lm-4)*6+(is-1);
                String filename = "multipols/" + fnamhigh[ifile];
                BufferedWriter out = new BufferedWriter(new FileWriter(filename, false));
                out.write(fnamhigh[ifile] + "\n");
                out.write("\n");
                out.close();
            }
        }
        } catch (IOException e) {
            System.err.println("Error writing lm=4,5 data: " + e.getMessage());
        }
        
        double Q2 = qq;
        int jj = 1;

        try {
        for (int iw = 1; iw <= Nw; iw++) {
            double W = W_min + (iw-1)*W_step;
            double El = (Math.pow(W,2)-Math.pow(mp,2))/2.0/mp;
            for (int lm = 4; lm <= 5; lm++) {
                for (int is = 1; is <= 3; is++) {
                    for (int ii = 1; ii <= 6; ii++) {
                        Rint[lm][is][ii][jj] = gauss("Finthigh", -1.0, 1.0, IntAccur);
                    }
                }
                for (int is = 1; is <= 3; is++) {
                    re_mult();
                }
                int lun = (lm+2)*10;
                write201(new BufferedWriter(new FileWriter("multipols/" + fnamhigh[(lm-4)*6+1])), W,El,Mlp,lm,jj);
                write201(new BufferedWriter(new FileWriter("multipols/" + fnamhigh[(lm-4)*6+2])), W,El,Mlm,lm,jj);
                write201(new BufferedWriter(new FileWriter("multipols/" + fnamhigh[(lm-4)*6+3])), W,El,Elp,lm,jj);
                write201(new BufferedWriter(new FileWriter("multipols/" + fnamhigh[(lm-4)*6+4])), W,El,Elm,lm,jj);
                write201(new BufferedWriter(new FileWriter("multipols/" + fnamhigh[(lm-4)*6+5])), W,El,Slp,lm,jj);
                write201(new BufferedWriter(new FileWriter("multipols/" + fnamhigh[(lm-4)*6+6])), W,El,Slm,lm,jj);
            }
        }
        
        } catch (IOException e) {
            System.err.println("Error writing jj=1 data: " + e.getMessage());
        }
   
        for (int lm = 4; lm <= 5; lm++) {
            for (int is = 1; is <= 6; is++) {
                int lun = (lm+2)*10+is;
            }
        }
    }

    // --- Write Fortran FORMAT(2x,2F7.3,2x,3F12.5) ---
    private static void write201(BufferedWriter out, double W, double El, float[][][] arr, int lm, int jj) throws IOException {
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
        int MAX_IS_INDEX = 3;
        int MAX_J_INDEX = 2;
        
        double[][][] PHis = new double[MAX_PHI_INDEX+1][MAX_IS_INDEX+1][MAX_J_INDEX+1];
            
        jrp.janrRePhi(w, q2, x, ierr);
            
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
                break;
            }
            
        fintg = 0.0;
            
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
                    break;
            }
            
            return fintg;          	
    }

    public double Finthigh(double x) {
    	
        double p0, p1, p2, p3, p4, p5, p6;
        double[][]      rim = new double[6][3];
        double[][]     PHis = new double[6][3];
        double[][]      pim = new double[6][2];
        double[][][]    pip = new double[6][5][2];
        double[][] PHisotop = new double[6][3];

        double extB, extR;
        int ierr = 0;
        
        jrp.janrRePhi(w, q2, x, ierr);

        extB = 1.0 / (1.0 + Math.pow(w * w - 1.16, 2));
        extR = extB * Math.pow(w * w - 1.16, 2);

        for (int ii = 0; ii < 6; ii++) {
            for (int ir = 0; ir < 3; ir++) {
                PHisotop[ii][ir] = PHisot[0][ii][ir] * extB + PHisot[1][ii][ir] * extR;
            }
        }

        // is=1 for (3/2), is=2 for (1/2), is=3 for (0)

        if (is == 3) {
            // (0), Fortran is=3 → Java is==3
            for (int ii = 0; ii < 6; ii++) {
                PHis[ii][is - 1] = PHisotop[ii][2];
            }
        } else if (is == 2) {
            // (1/2)
            for (int ii = 0; ii < 6; ii++) {
                PHis[ii][is - 1] = PHisotop[ii][0] + 2.0 * PHisotop[ii][1];
            }
        } else if (is == 1) {
            // (3/2)
            for (int ii = 0; ii < 6; ii++) {
                PHis[ii][is - 1] = PHisotop[ii][0] - PHisotop[ii][1];
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

        // rim(1,is) = PHis(3,is)+x*PHis(4,is)
        // Fortran: rim(1,is) → Java: rim[0][is-1]
        if (is >= 1 && is <= 3) {
            rim[0][is - 1] = PHis[2][is - 1] + x * PHis[3][is - 1];
            rim[1][is - 1] = 2.0 * PHis[0][is - 1] - 2.0 * x * PHis[1][is - 1] + p0 * PHis[3][is - 1];
            rim[2][is - 1] = PHis[1][is - 1];
            rim[3][is - 1] = PHis[3][is - 1];
            rim[4][is - 1] = PHis[4][is - 1] + x * PHis[5][is - 1];
            rim[5][is - 1] = PHis[5][is - 1];

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

            // For lm==5 or lm==4, Fortran is 1-based (lm=4 or 5), Java index is lm-1
            if (lm == 5) {
                for (int ii = 0; ii < 6; ii++) {
                    pip[ii][lm - 1][is - 1] = pim[ii][1] * rim[ii][is - 1];
                }
            } else if (lm == 4) {
                for (int ii = 0; ii < 6; ii++) {
                    pip[ii][lm - 1][is - 1] = pim[ii][0] * rim[ii][is - 1];
                }
            }
        }

        double finthigh = 0.0;

        switch (i) {
            case 1:
                // Integral 1 (I1)
                finthigh = pip[0][lm - 1][is - 1];
                break;
            case 2:
                // Integral 2 (I2)
                finthigh = pip[1][lm - 1][is - 1];
                break;
            case 3:
                // Integral 3 (I3)
                finthigh = pip[2][lm - 1][is - 1];
                break;
            case 4:
                // Integral 4 (I4)
                finthigh = pip[3][lm - 1][is - 1];
                break;
            case 5:
                // Integral 5 (I5)
                finthigh = pip[4][lm - 1][is - 1];
                break;
            case 6:
                // Integral 6 (I6)
                finthigh = pip[5][lm - 1][is - 1];
                break;
            default:
                // Invalid i value
                throw new IllegalArgumentException("Parameter i out of range (1..6), got i=" + i);
        }

        return finthigh;	
    }
    
    public void re_mult() {
    	
        int jj = j;

        int lmIdx = lm;
        int isIdx = is;
        int jjIdx = jj;

        if (lm >= 0) {
            Elp[lmIdx][isIdx][jjIdx] = (float) (
                (-2.0 * lm * Rint[lmIdx][isIdx][1][jjIdx]
                           + Rint[lmIdx][isIdx][2][jjIdx]
                + 2.0 * lm * Rint[lmIdx][isIdx][3][jjIdx]
           - lm * (lm - 1) * Rint[lmIdx][isIdx][4][jjIdx])
                / (2.0 * (lm + 1.0) * (2.0 * lm + 1.0))
            );
            Slp[lmIdx][isIdx][jjIdx] = (float) (
                (Rint[lmIdx][isIdx][5][jjIdx]
          - lm * Rint[lmIdx][isIdx][6][jjIdx])
                / (lm + 1.0) / (2.0 * lm + 1.0)
            );
            Mlp[lmIdx][isIdx][jjIdx] = 0.0f;
            Mlm[lmIdx][isIdx][jjIdx] = 0.0f;
            Elm[lmIdx][isIdx][jjIdx] = 0.0f;
            Slm[lmIdx][isIdx][jjIdx] = 0.0f;
        }

        if (lm >= 1) {
            Mlp[lmIdx][isIdx][jjIdx] = (float) (
                (2.0 *                 Rint[lmIdx][isIdx][1][jjIdx]
                                     + Rint[lmIdx][isIdx][2][jjIdx]
                          + 2.0 * lm * Rint[lmIdx][isIdx][3][jjIdx]
           + (lm - 1.0) * (lm + 2.0) * Rint[lmIdx][isIdx][4][jjIdx])
                / (2.0 * (lm + 1.0) * (2.0 * lm + 1.0))
            );
            Mlm[lmIdx][isIdx][jjIdx] = (float) (
                (-2.0 *                 Rint[lmIdx][isIdx][1][jjIdx]
                                      - Rint[lmIdx][isIdx][2][jjIdx]
                     + 2.0 * (lm + 1) * Rint[lmIdx][isIdx][3][jjIdx]
            - (lm - 1.0) * (lm + 2.0) * Rint[lmIdx][isIdx][4][jjIdx])
                / (2.0 * lm * (2.0 * lm + 1.0))
            );
            Slm[lmIdx][isIdx][jjIdx] = (float) (
                (Rint[lmIdx][isIdx][5][jjIdx]
                 - (lm + 1.0) * Rint[lmIdx][isIdx][6][jjIdx])
                / (lm * (2.0 * lm + 1.0))
            );
        }

        if (lm >= 2) {
            Elm[lmIdx][isIdx][jjIdx] = (float) (
                (2.0 *      (lm + 1.0) * Rint[lmIdx][isIdx][1][jjIdx]
                                       + Rint[lmIdx][isIdx][2][jjIdx]
                    - 2.0 * (lm + 1.0) * Rint[lmIdx][isIdx][3][jjIdx]
             - (lm + 1.0) * (lm + 2.0) * Rint[lmIdx][isIdx][4][jjIdx])
                / (2.0 * lm * (2.0 * lm + 1.0))
            );
        }
    }

    public static void main(String[] args) {
    	MultAna ma = new MultAna();
        double result = ma.gauss("ff", 0, 1, 1e-6);
        System.out.println("Integration result: " + result); 
    }
}
