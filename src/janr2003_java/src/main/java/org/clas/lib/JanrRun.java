package org.clas.lib;

public class JanrRun extends Constants {

    private static double wsave = -1.0;
    private static double qsave = -1.0;
    private static double csave = -2.0;
    private static double phid  =  0.0;
    
    JanrIniPoint jini = new JanrIniPoint();
    POLINT          p = new POLINT();
    boolean      test = false;
    
    public double[] Breit = new double[2];
    public double    rv[] = new double[1];
        
    public JanrRun() {
    	
    }

    public void janrRun(double wi,    
                        double Q2i,
                        double epsilon,
                        double costhi,
                        double phidi,
                        boolean iupd) {

        double phi = (float) (phidi * pi / 180.0);
        double q2,q22,k22, qpion, qkratio, kgamma, E2pion, Eeta; 
        double s=0, w=0, theta=0, costh=0;

        boolean okWsave = (wi   == wsave);
        boolean okQsave = (Q2i  == qsave);
        boolean okCsave = (costhi == csave);
        boolean skip    = okWsave && okCsave && okQsave;
   
//        if (!okQsave) jini.janrIniPoint(Q2i);

        wsave = wi;
        qsave = Q2i;
        csave = costhi;
        phid  = phidi;
        
        if (!( !iupd && skip )) {  	
            w     = wi;
            q2    = Q2i;            
            costh = costhi;
            
            jini.janrIniPoint(q2);
            
            theta = (float) Math.acos(costh);
            theta = (float) (theta * 180.0 / pi);
            
            s       = w * w;
            w2      = 2 * w;
            E       = s - mn22;
            Egamma  = (E - q2) / w2;
            Epion   = (E + mp22) / w2;
            q22     = Epion * Epion - mp22;
            k22     = (w - Egamma) * (w - Egamma) - mn22;
            qpion   = Math.sqrt(Math.max(0.0, q22));
            qkratio = qpion * w2 / E;
            kgamma  = Math.sqrt(Math.max(0.0, k22));
            qk      = qpion * kgamma * costh - Epion * Egamma;
            t       = (mp22 - q2 + 2.0 * qk);
            u       = 2.0 * mn22 + mp22 - q2 - s - t;
            E2pion  = (E + m2pion*m2pion) / w2;
            q2pion  = E2pion * E2pion - m2pion*m2pion;
            Eeta    = (E + meta*meta) / w2;
            qeta    = Eeta * Eeta - meta*meta;

            phiAmpl(s, costh, q2, kgamma, qpion);
            helAmpl(q2, theta, kgamma);
        } else {
            // We still need s and q2 for the cross_section call
            s  = wsave * wsave;
            q2 = qsave;
        }

        crossSections(s, q2, phi, epsilon);

    }

    public void phiAmpl(double s,
                        double costh,
                        double q2,
                        double kgamma,
                        double qpion) {

        double[][] E1 = new double[2][2];
        double[]   E3 = new double[2]; //what was this?
        double p1, p2, p3, p4, p5;
        
        multAmpl(s, q2, kgamma, qpion);
        highMult(s, q2, costh);
        
        // Calculate derivatives of Legendre polynomials
        // p1=P2', p2=P3', p3=P4', p4=P3", p5=P4"
        
        p1 = 3.0f * costh;
        p2 = (15.0f * costh * costh - 3.0f) / 2.0f;
        p3 = 2.5f * (7.0f * costh * costh * costh - 3.0f * costh);
        p4 = 15.0f * costh;
        p5 = 7.5f * (7.0f * costh * costh - 1.0f);
        
        // Calculate intermediate amplitudes for transition from 
        // multipole amplitudes to helicity amplitudes and cross sections.
        // Multipole ampl. are normalized according to CGLN Phys.Rev. 106 (1957) 1345.
        // They are equal to mult.ampl. of Devenish&Lyth/8W/pi
        // mult.ampl. (1) => (0+)
        // mult.ampl. (2) => (1-)
        // mult.ampl. (3) => (1+)
        // mult.ampl. (4) => (2-)
        // mult.ampl. (5) => (2+)
        // mult.ampl. (6) => (3-)
        // mult.ampl. (7) => (3+)
        
        for (int i = 0; i < 2; i++) { // Re(0) Im(1)
            for (int j = 0; j < 2; j++) { // pi0(0) pi+(1)
      	
                E1[i][j] = ReacM6[i][j] + ReacE6[i][j];
                
                phi_amp[0][i][j] = ReacE1[i][j] 
                    + p1 *        (ReacM3[i][j] + ReacE3[i][j])
                    + p3 *  (3.0 * ReacM7[i][j] + ReacE7[i][j])
                    +        3.0 * ReacM4[i][j] + ReacE4[i][j]
                    + p1 *  (4.0 * ReacM6[i][j] + ReacE6[i][j])
                    + p2 *  (2.0 * ReacM5[i][j] + ReacE5[i][j]);
                
                phi_amp[1][i][j] = 2.0 * ReacM3[i][j] +        ReacM2[i][j]
                    + p1 *        (2.0 * ReacM4[i][j] + 3.0f * ReacM5[i][j])
                    + p2 *         3.0 * ReacM6[i][j]
                    + p2 *         4.0 * ReacM7[i][j];
                
                phi_amp[2][i][j] = 3.0 * (ReacE3[i][j] - ReacM3[i][j] + E1[i][j]) 
                    + p5 *               (ReacE7[i][j] - ReacM7[i][j])
                    + p4 *               (ReacE5[i][j] - ReacM5[i][j]);
                
                phi_amp[3][i][j] = -3.0 * (ReacM4[i][j] + ReacE4[i][j])
                    +               3.0 * (ReacM5[i][j] - ReacE5[i][j])
//                  - p4 *                 ReacM6[i][j] + ReacE6[i][j]    //no parentheses
                    - p4 *                     E1[i][j]                   //lcs 8.29.2025 
                    + p4 *                (ReacM7[i][j] - ReacE7[i][j]);
                
                phi_amp[4][i][j] = ReacS1[i][j] + 2.0f * p1 * ReacS3[i][j]
                    + 4.0 * p3 *   ReacS7[i][j] + 3.0f * p2 * ReacS5[i][j]
                    - 2.0 *        ReacS4[i][j] - 3.0f * p1 * ReacS6[i][j];
                
                phi_amp[5][i][j] = ReacS2[i][j] + 2.0f * p1 * ReacS4[i][j]
                    + 3.0f * p2 *  ReacS6[i][j] - 2.0f *      ReacS3[i][j]
                    - 4.0f * p2 *  ReacS7[i][j] - 3.0f * p1 * ReacS5[i][j];
            }
        }
        
        for (int i = 0; i < 6; i++) {
            for (int j = 0; j < 2; j++) {
                phi_amp[i][0][j] += phih[i][j]; //from HighMult
            }
        }
    }
    
    public void multAmpl(double s, double q2, double kgamma, double qpion) {
        
        double[][]   Br1 = new double[2][9];
        double[][]   Br3 = new double[2][8];
        double[][]   EBr1 = new double[2][9]; //EBAC phase
        double[][]   EBr3 = new double[2][9]; //EBAC phase
        double[][][] BrM1 = new double[2][2][9];
        double[][][] BrE1 = new double[2][2][9];
        double[][][] BrS1 = new double[2][2][9];
        double[][]   BrM3 = new double[2][8];
        double[][]   BrE3 = new double[2][8];
        double[][]   BrS3 = new double[2][8];
        
        double RMP33, IMP33, coeff;
        double XIN, XIN1, elab, WW, Suppress;
        double YOUTSimr, YOUTSimi;
        double[] mpion = new double[2];
        double[] qpi = new double[2];
        double ratio;
        
        breitAmpl(s, q2, kgamma, qpion);
        background(s, q2);
        backHigh(s);
        
        WW = Math.sqrt(s);
        mpion[0] = 0.13498;
        mpion[1] = 0.13957;
        
        for (int i = 0; i < 2; i++) {
            qpi[i] = Math.sqrt(Math.pow(((s-mn22 + mpion[i]*mpion[i])/2.0/WW), 2) - mpion[i]*mpion[i]);
        }
        
        for (int j = 0; j < 9; j++) {
            for (int i = 0; i < 2; i++) {
                EBr1[i][j] = 1.0;
                if (j < 8) EBr3[i][j] = 1.0; //not used
            }
        }
        
        // Used below to correct pi+ multipoles for pion mass ratio
        // All calculations in isospin channels use pi0 mass by default
        
        ratio = qpi[1] / qpi[0];
        
        XIN = WW;
        elab = (s - mn * mn) / 2.0 / mn;
        
        if (WW <= 1.34) {
            YOUTSimr = interp(WSa, simr, 13, 3, XIN);
            YOUTSimi = interp(WSa, simi, 13, 3, XIN);
            
            for (int j = 0; j < 2; j++) {
                BreitE1[1][j][1] = BreitE1[0][j][1] * YOUTSimr / (1.0 - YOUTSimi);
                BreitS1[1][j][1] = BreitS1[0][j][1] * YOUTSimr / (1.0 - YOUTSimi);
            }
        }
        
        // Breit-Wigner amplitudes I=1/2
        for (int i = 0; i < 2; i++) {         // Re(0) Im(1)
            for (int j = 0; j < 2; j++) {     // pA(1/2) (j=0) nA(1/2) (j=1)
                BrM1[i][j][5] = BreitM1[i][j][5] + BreitM1[i][j][12];     
                BrE1[i][j][5] = BreitE1[i][j][5] + BreitE1[i][j][12];
                BrS1[i][j][5] = BreitS1[i][j][5] + BreitS1[i][j][12];

                BrM1[i][j][6] = BreitM1[i][j][6] + BreitM1[i][j][10];   
                BrE1[i][j][6] = BreitE1[i][j][6] + BreitE1[i][j][10];
                BrS1[i][j][6] = BreitS1[i][j][6] + BreitS1[i][j][10];

                BrE1[i][j][1] = BreitE1[i][j][1] + BreitE1[i][j][4]; //[1] is wrong     
                BrS1[i][j][1] = BreitS1[i][j][1] + BreitS1[i][j][4]; //[1] is wrong

                BrM1[i][j][0] = BreitM1[i][j][0] + BreitM1[i][j][9];    
                BrS1[i][j][0] = BreitS1[i][j][0] + BreitS1[i][j][9];

                BrM1[i][j][2] = BreitM1[i][j][2] + BreitM1[i][j][3] + BreitM1[i][j][11];  
                BrE1[i][j][2] = BreitE1[i][j][2] + BreitE1[i][j][3] + BreitE1[i][j][11];
                BrS1[i][j][2] = BreitS1[i][j][2] + BreitS1[i][j][3] + BreitS1[i][j][11];

                BrM1[i][j][8] = BreitM1[i][j][8];                        
                BrE1[i][j][8] = BreitE1[i][j][8];
                BrS1[i][j][8] = BreitS1[i][j][8];

                BrM1[i][j][7] = BreitM1[i][j][7];                        
                BrE1[i][j][7] = BreitE1[i][j][7];
                BrS1[i][j][7] = BreitS1[i][j][7];
            }
        }
        
        // 1+ P33 M1+
        XIN1 = XIN;
        RMP33 = interp(WSa, mp33r, 51, 4, XIN1);
        IMP33 = interp(WSa, mp33i, 51, 4, XIN1);
        
        // 1- P11 1(1)
        Br1[0][0] = interp(WSa, p11r, 51, 3, XIN);
        Br1[1][0] = interp(WSa, p11i, 51, 3, XIN);
        
        if (nebac > 0) {
//            EBr1[0][0] = interp(EWSa, ep11r, 48, 3, XIN);
//            EBr1[1][0] = interp(EWSa, ep11i, 48, 3, XIN);
        }
        
        // 0+ S11 1(2)
        Br1[0][1] = interp(WSa, s11r, 51, 3, XIN);
        Br1[1][1] = interp(WSa, s11i, 51, 3, XIN);
        
        // 2- D13 1(3)
        Br1[0][2] = interp(WSa, d13r, 51, 3, XIN);
        Br1[1][2] = interp(WSa, d13i, 51, 3, XIN);
        
        Br1[0][3] = 0.0;
        Br1[1][3] = 0.0;
        Br1[0][4] = 0.0;
        Br1[1][4] = 0.0;
        Br1[0][7] = 0.0;
        Br1[1][7] = 0.0;
        
        // 3- F15 1(6)
        Br1[0][5] = interp(WSa, f15r, 51, 3, XIN);
        Br1[1][5] = interp(WSa, f15i, 51, 3, XIN);
        
        // 1+ P13 1(7)
        Br1[0][6] = interp(WSa, p13r, 51, 3, XIN);
        Br1[1][6] = interp(WSa, p13i, 51, 3, XIN);
        
        // 2+ D15 1(9)
        Br1[0][8] = interp(WSa, d15r, 51, 3, XIN);
        Br1[1][8] = interp(WSa, d15i, 51, 3, XIN);
        
        // 1+ P33 3(1)
        Br3[0][0] = interp(WSa, p33r, 51, 3, XIN) - 0.05;
        Br3[1][0] = interp(WSa, p33i, 51, 3, XIN) * 1.017;
        
        // 2- D33 3(2)
        Br3[0][1] = interp(WSa, d33r, 51, 3, XIN);
        Br3[1][1] = interp(WSa, d33i, 51, 3, XIN);
        
        // 3+ F37 3(3)
        Br3[0][2] = interp(WSa, f37r, 51, 3, XIN);
        Br3[1][2] = interp(WSa, f37i, 51, 3, XIN);
        
        Br3[0][3] = 0.0;
        Br3[1][3] = 0.0;
        
        // 0+ S31 3(5)
        Br3[0][4] = interp(WSa, s31r, 51, 3, XIN);
        Br3[1][4] = interp(WSa, s31i, 51, 3, XIN);
        
        Br3[0][5] = 0.0;
        Br3[1][5] = 0.0;
        
        Br3[0][6] = 0.0;
        Br3[1][6] = 0.0;
        
        // 1- P31 3(8)
        Br3[0][7] = interp(WSa, p31r, 51, 3, XIN);
        Br3[1][7] = interp(WSa, p31i, 51, 3, XIN);
        
        //n1=1 (background on) n2=1 (resonance on) EBr1=EBAC phase
        for (int j = 0; j < 9; j++) {
            for (int i = 0; i < 2; i++) { // pA(1/2) (i=0) nA(1/2) (i=1)           	
                MM1[0][i][j] = n1 * BM1[i][j] * (1.0 - Br1[1][j]) + n2 * BrM1[0][i][j] * EBr1[0][j];
                ME1[0][i][j] = n1 * BE1[i][j] * (1.0 - Br1[1][j]) + n2 * BrE1[0][i][j] * EBr1[0][j];
                MS1[0][i][j] = n1 * BS1[i][j] * (1.0 - Br1[1][j]) + n2 * BrS1[0][i][j] * EBr1[0][j];
                MM1[1][i][j] = n1 * BM1[i][j] *        Br1[0][j]  + n2 * BrM1[1][i][j] * EBr1[1][j];
                ME1[1][i][j] = n1 * BE1[i][j] *        Br1[0][j]  + n2 * BrE1[1][i][j] * EBr1[1][j];
                MS1[1][i][j] = n1 * BS1[i][j] *        Br1[0][j]  + n2 * BrS1[1][i][j] * EBr1[1][j];
            }
        }
        
        // Breit-Wigner amplitudes I=3/2
        for (int i = 0; i < 2; i++) { // Re(0) Im(1)
            BrM3[i][0] = BreitM3[i][0] + BreitM3[i][3] + BreitM3[i][9];
            BrE3[i][0] = BreitE3[i][0] + BreitE3[i][3] + BreitE3[i][9];
            BrS3[i][0] = BreitS3[i][0] + BreitS3[i][3] + BreitS3[i][9];

            BrM3[i][1] = BreitM3[i][1] + BreitM3[i][11];
            BrE3[i][1] = BreitE3[i][1] + BreitE3[i][11];
            BrS3[i][1] = BreitS3[i][1] + BreitS3[i][11];

            BrM3[i][2] = BreitM3[i][2];
            BrE3[i][2] = BreitE3[i][2];
            BrS3[i][2] = BreitS3[i][2];

            BrE3[i][4] = BreitE3[i][4] + BreitE3[i][8];
            BrS3[i][4] = BreitS3[i][4] + BreitS3[i][8]; // BrS3[*][4] should be zero

            BrM3[i][5] = BreitM3[i][5] + BreitM3[i][12];
            BrE3[i][5] = BreitE3[i][5] + BreitE3[i][12];
            BrS3[i][5] = BreitS3[i][5] + BreitS3[i][12];

            BrM3[i][6] = BreitM3[i][6];
            BrE3[i][6] = BreitE3[i][6];
            BrS3[i][6] = BreitS3[i][6];

            BrM3[i][7] = BreitM3[i][7] + BreitM3[i][10];
            BrS3[i][7] = BreitS3[i][7] + BreitS3[i][10];
        }
        
        // Construct Re(0) Im(1) parts of (1+) multipole amplitudes for D(1232)
        MM3[0][0] = n1 * BM3[0] * (1.0 - Br3[1][0]) + n2 * BrM3[0][0];
        ME3[0][0] = n1 * BE3[0] * (1.0 - Br3[1][0]) + n2 * BrE3[0][0];
        MS3[0][0] = n1 * BS3[0] * (1.0 - Br3[1][0]) + n2 * BrS3[0][0];
        MM3[1][0] = n1 * BM3[0] *        Br3[0][0]  + n2 * BrM3[1][0];
        ME3[1][0] = n1 * BE3[0] *        Br3[0][0]  + n2 * BrE3[1][0];
        MS3[1][0] = n1 * BS3[0] *        Br3[0][0]  + n2 * BrS3[1][0];
        
        // Note RMP33 and IMP33 are previously defined shape for M1+ which
        // was determined from fits to Q2=0 VPI M1+ data in PRC 67, 015209 (2003).
        // Shape follows description of modified B-W and Regge modified backgrounds
        // described in Eq. 16-17 of paper.  resM is GM* parameterization defined 
        // in janr_ini_point.F   
        
        MM3[0][0] = resM * RMP33;
        MM3[1][0] = resM * IMP33;
        
        // Empirical adjustment to shape of D(1232) away from resonance position
        // These parameters are intended to be used for Q2=0.4 fits.
        // This commented code is omitted in this java version
        
        if (WW <= 1.2283) {
            for (int i = 0; i < 2; i++) {
                MM3[i][0] = MM3[i][0] * Math.pow(WW / 1.2283, cm2);
                ME3[i][0] = ME3[i][0] * Math.pow(WW / 1.2283, ce2);
                MS3[i][0] = MS3[i][0] * Math.pow(WW / 1.2283, cs2);
            }
        } else {
            for (int i = 0; i < 2; i++) {
                MM3[i][0] = MM3[i][0] * Math.pow(WW / 1.2283, cm3);
                ME3[i][0] = ME3[i][0] * Math.pow(WW / 1.2283, ce3);
                MS3[i][0] = MS3[i][0] * Math.pow(WW / 1.2283, cs3);
            }
        }
        
        // Construct Re(0) Im(1) parts of other multipole amplitudes A(3/2)
        for (int j = 1; j < 8; j++) {       	 
            MM3[0][j] = n1 * BM3[j] * (1.0 - Br3[1][j]) + n2 * BrM3[0][j];
            ME3[0][j] = n1 * BE3[j] * (1.0 - Br3[1][j]) + n2 * BrE3[0][j];
            MS3[0][j] = n1 * BS3[j] * (1.0 - Br3[1][j]) + n2 * BrS3[0][j]; 
            MM3[1][j] = n1 * BM3[j] *        Br3[0][j]  + n2 * BrM3[1][j];
            ME3[1][j] = n1 * BE3[j] *        Br3[0][j]  + n2 * BrE3[1][j];
            MS3[1][j] = n1 * BS3[j] *        Br3[0][j]  + n2 * BrS3[1][j]; 
        }
        
        coeff = 1.0; // units mFm/10 ?
        
        for (int j = 0; j < 9; j++) {
            for (int i = 0; i < 2; i++) {
                for (int k = 0; k < 2; k++) {
                    MM1[k][i][j] = coeff * MM1[k][i][j];
                    ME1[k][i][j] = coeff * ME1[k][i][j];
                    MS1[k][i][j] = coeff * MS1[k][i][j];
                }
            }
        }
        
        for (int j = 0; j < 8; j++) {
            for (int k = 0; k < 2; k++) {
                MM3[k][j] = coeff * MM3[k][j];
                ME3[k][j] = coeff * ME3[k][j];
                MS3[k][j] = coeff * MS3[k][j];
            }
        }
        
        // ReacM1(i,j)=0+	i=0,1 for Re,Im;	j=0,1 for pi0,pi+	 
        // ReacM2(i,j)=1-
        // ReacM3(i,j)=1+
        // ReacM4(i,j)=2-
        // ReacM5(i,j)=2+
        // ReacM6(i,j)=3-
        // ReacM7(i,j)=3+
        
        // Reaction multipoles
        for (int i = 0; i < 2; i++) { //i=0: Re i=1: Im
            ReacE1[i][0] =                   ME1[i][0][1] + ME3[i][4] / 1.5;
            ReacE1[i][1] = Math.sqrt(2.0) * (ME1[i][0][1] - ME3[i][4] / 3.0);
            ReacS1[i][0] =                   MS1[i][0][1] + MS3[i][4] / 1.5;          //bad value for i=0
            ReacS1[i][1] = Math.sqrt(2.0) * (MS1[i][0][1] - MS3[i][4] / 3.0);
            
            ReacM2[i][0] =                   MM1[i][0][0] + MM3[i][7] / 1.5;
            ReacM2[i][1] = Math.sqrt(2.0) * (MM1[i][0][0] - MM3[i][7] / 3.0) * ratio;
            ReacS2[i][0] =                   MS1[i][0][0] + MS3[i][7] / 1.5;          //bad value for i=0
            ReacS2[i][1] = Math.sqrt(2.0) * (MS1[i][0][0] - MS3[i][7] / 3.0) * ratio; //bad value for i=1
            
            ReacM3[i][0] =                   MM1[i][0][6] + MM3[i][0] / 1.5;
            ReacM3[i][1] = Math.sqrt(2.0) * (MM1[i][0][6] - MM3[i][0] / 3.0) * ratio;
            ReacE3[i][0] =                   ME1[i][0][6] + ME3[i][0] / 1.5;
            ReacE3[i][1] = Math.sqrt(2.0) * (ME1[i][0][6] - ME3[i][0] / 3.0) * ratio;
            ReacS3[i][0] =                   MS1[i][0][6] + MS3[i][0] / 1.5;
            ReacS3[i][1] = Math.sqrt(2.0) * (MS1[i][0][6] - MS3[i][0] / 3.0) * ratio;
            
            ReacM4[i][0] =                   MM1[i][0][2] + MM3[i][1] / 1.5;
            ReacM4[i][1] = Math.sqrt(2.0) * (MM1[i][0][2] - MM3[i][1] / 3.0) * ratio * ratio;
            ReacE4[i][0] =                   ME1[i][0][2] + ME3[i][1] / 1.5;
            ReacE4[i][1] = Math.sqrt(2.0) * (ME1[i][0][2] - ME3[i][1] / 3.0) * ratio * ratio;
            ReacS4[i][0] =                   MS1[i][0][2] + MS3[i][1] / 1.5;
            ReacS4[i][1] = Math.sqrt(2.0) * (MS1[i][0][2] - MS3[i][1] / 3.0) * ratio * ratio;
            
            ReacM5[i][0] =                   MM1[i][0][8] + MM3[i][6] / 1.5;
            ReacM5[i][1] = Math.sqrt(2.0) * (MM1[i][0][8] - MM3[i][6] / 3.0) * ratio * ratio;
            ReacE5[i][0] =                   ME1[i][0][8] + ME3[i][6] / 1.5;
            ReacE5[i][1] = Math.sqrt(2.0) * (ME1[i][0][8] - ME3[i][6] / 3.0) * ratio * ratio;
            ReacS5[i][0] =                   MS1[i][0][8] + MS3[i][6] / 1.5;
            ReacS5[i][1] = Math.sqrt(2.0) * (MS1[i][0][8] - MS3[i][6] / 3.0) * ratio * ratio;
            
            ReacM6[i][0] =                   MM1[i][0][5] + MM3[i][5] / 1.5;
            ReacM6[i][1] = Math.sqrt(2.0) * (MM1[i][0][5] - MM3[i][5] / 3.0) * ratio * ratio * ratio;
            ReacE6[i][0] =                   ME1[i][0][5] + ME3[i][5] / 1.5;
            ReacE6[i][1] = Math.sqrt(2.0) * (ME1[i][0][5] - ME3[i][5] / 3.0) * ratio * ratio * ratio;
            ReacS6[i][0] =                   MS1[i][0][5] + MS3[i][5] / 1.5;
            ReacS6[i][1] = Math.sqrt(2.0) * (MS1[i][0][5] - MS3[i][5] / 3.0) * ratio * ratio * ratio;
            
            ReacM7[i][0] =                   MM1[i][0][7] + MM3[i][2] / 1.5;
            ReacM7[i][1] = Math.sqrt(2.0) * (MM1[i][0][7] - MM3[i][2] / 3.0) * ratio * ratio * ratio;
            ReacE7[i][0] =                   ME1[i][0][7] + ME3[i][2] / 1.5;
            ReacE7[i][1] = Math.sqrt(2.0) * (ME1[i][0][7] - ME3[i][2] / 3.0) * ratio * ratio * ratio;
            ReacS7[i][0] =                   MS1[i][0][7] + MS3[i][2] / 1.5;
            ReacS7[i][1] = Math.sqrt(2.0) * (MS1[i][0][7] - MS3[i][2] / 3.0) * ratio * ratio * ratio;
        }
    }
    
    public double interp(double[] x, double[] y, int npts, int nterms, double xin) {
        int i1 = 0, i2 = 0, k, j, ix, imax, ixmax;
        double denom, deltax, sum, prod;
        double[] delta = new double[10];
        double[] a = new double[10];
        double yout = 0;
        
        // Find the appropriate range of points for interpolation
        for (int i = 0; i < npts; i++) {
            if (xin < x[i]) {
                i1 = i - nterms / 2;
                if (i1 < 0) {
                    i1 = 0;
                }
                break;
            } else if (xin == x[i]) {
                return y[i];  // Exact match found
            }
            
            // If we reach the end without finding x[i] > xin
            if (i == npts - 1) {
                i1 = npts - nterms;
            }
        }
        
        // Adjust i1 and i2 to stay within bounds
        i2 = i1 + nterms - 1;
        if (i2 >= npts) {
            i2 = npts - 1;
            i1 = i2 - nterms + 1;
            if (i1 < 0) {
                i1 = 0;
                nterms = i2 - i1 + 1;
            }
        }
        
        // Calculate delta values
        denom = x[i1 + 1] - x[i1];
        deltax = (xin - x[i1]) / denom;
        for (int i = 0; i < nterms; i++) {
            ix = i1 + i;
            delta[i] = (x[ix] - x[i1]) / denom;
        }
        
        // Calculate coefficients
        a[0] = y[i1];
        for (k = 1; k < nterms; k++) {
            prod = 1.0f;
            sum = 0.0f;
            imax = k;
            ixmax = i1 + imax;
            
            for (int i = 0; i < imax; i++) {
                j = k - i - 1;
                prod = prod * (delta[k] - delta[j]);
                sum = sum - a[j] / prod;
            }
            a[k] = sum + y[ixmax] / prod;
        }
        
        // Calculate the interpolated value
        sum = a[0];
        for (j = 1; j < nterms; j++) {
            prod = 1.0f;
            imax = j;
            for (int i = 0; i < imax; i++) {
                prod = prod * (deltax - delta[i]);
            }
            sum = sum + a[j] * prod;
        }
        
        yout = sum;
        return yout;
    }
    
    private void breitAmpl(double s, double q2, double kgamma, double qpion) {
    	
//        double[] Breit = new double[2];
//        double    rv[] = new double[1];
        
        // Isospin 3/2 resonances (i=1 to 6)
        for (int i = 0; i < 6; i++) {
            Breit = getBreit(s, q2, kgamma, qpion, i+1, rv);   
            for (int j = 0; j < 2; j++) { // real (j=0), imag (j=1)
                BreitE3[j][i] = Breit[j] * ae3[i];
                BreitM3[j][i] = Breit[j] * am3[i];
                BreitS3[j][i] = Breit[j] * as3[i] * rv[0]; //as3[4,8] should be zero in janrIniRes
           }

        }
        
        // Zero out remaining elements (i=7 to 13)
//        for (int i = 6; i < 13; i++) {
        for (int i = 7; i < 13; i++) {
            for (int j = 0; j < 2; j++) {
                BreitE3[j][i] = 0.0f;
                BreitM3[j][i] = 0.0f;
                BreitS3[j][i] = 0.0f;
            }
        }
        
        // Isospin 1/2 resonances (i=1 to 10)
        for (int i = 0; i < 10; i++) {
            int ii = i + 21; 
            Breit = getBreit(s, q2, kgamma, qpion, ii, rv); 
            for (int k = 0; k < 2; k++) { // pA(1/2) (k=0), nA(1/2) (k=1)
                for (int j = 0; j < 2; j++) { // real (j=0), imag (j=1)
                    BreitE1[j][k][i] = Breit[j] * ae1[k][i];
                    BreitM1[j][k][i] = Breit[j] * am1[k][i];
                    BreitS1[j][k][i] = Breit[j] * as1[k][i] * rv[0];
//                    if (test && (i==0||i==1||i==2)) {
//                    	System.out.println(j+" "+k+" "+i+" "+(float)BreitE1[j][k][i]+" "+(float)Breit[j]+" "+ae1[k][i]);
//                    	System.out.println(j+" "+k+" "+i+" "+(float)BreitS1[j][k][i]+" "+(float)Breit[j]+" "+as1[k][i]+" "+rv[0]);
//                    }
                }
            }
        }
        
        // Zero out remaining elements (i=11 to 13)
//        for (int i = 10; i < 13; i++) {
        for (int i = 11; i < 13; i++) {
            for (int k = 0; k < 2; k++) {
                for (int j = 0; j < 2; j++) {
                    BreitE1[j][k][i] = 0.0f;
                    BreitM1[j][k][i] = 0.0f;
                    BreitS1[j][k][i] = 0.0f;
                }
            }
        }
    } 
    
    private double[] getBreit(double s, double q2, double kgamma, double qpion, int i, double[] rv) {
        
    	double wthr,breta,eeta,qreta,reta;
    	double ww,q22,k22,kkr,x2,w2pion;
    	double m22,m2,en1,en2,kr22,qr22,kr,qr;
    	double e2pion,qr2pion,r2pion,rpion,qxl,qx2l;
    	double gpion,geta,ginel,gtotal,qxgamma,ggamma;
    	double den,breitw,f1,f2,a,jr;
    	
    	int ll,lp;
    	
        a     = 1.0;
        wthr  = 1.0;
        breta = 0.0;
        reta  = 0;
        
        ww  = Math.sqrt(s);
        m22 = mres[i]*mres[i];
        m2  = 2 * mres[i];
        
//        double mn22 = Math.pow(mn, 2);
//        double mp22 = Math.pow(mp, 2);
//        double m2pion = 0.13957; // Assuming pion mass if not provided
        
        // Eta channel calculation
        if (i == 22) {
            breta = 0.55;
            if (ww < (mn + meta)) wthr = 0.0;
            eeta  = (m22 + meta*meta - mn22) / m2;
            qreta = eeta*eeta - meta*meta;
//            qeta  = Math.sqrt(Math.max(0, eeta*eeta - meta*meta));
            reta  = Math.abs(qeta / qreta);
            reta  = Math.sqrt(reta);
        }
        
        q22 = qpion*qpion;
        k22 = kgamma*kgamma;
        
        ll = (int) langul[i];
        lp = (int) lprime[i];
        jr = jres[i] * 0.5;
        x2 = Math.pow(xres[i], 2);
        
        if (ww <= (mn + 2.0 * mp)) q2pion = 0.0;

        kkr = (m22 - mn22) / m2;
        
        en1 = (m22 + mn22 + q2) / m2;
        en2 = (m22 + mn22 - mp22) / m2;
        
        kr22 = en1*en1 - mn22;
        qr22 = en2*en2 - mn22;
        
        kr = Math.sqrt(kr22);
        qr = Math.sqrt(qr22);
        
        rv[0] = kgamma / kr;

        rpion = qpion / qr;
        
        e2pion = (m22 + Math.pow(m2pion, 2) - mn22) / m2;
        qr2pion = Math.pow(e2pion, 2) - Math.pow(m2pion, 2);
        r2pion = Math.pow(q2pion / qr2pion, ll + 2);
        
        qxl  = Math.pow((x2 + qr22)    / (x2 + q22),    ll);
        qx2l = Math.pow((x2 + qr2pion) / (x2 + q2pion), ll + 2);
        
        gpion = eta[i] * gres[i] * Math.pow(rpion, 2.0 * ll + 1.0) * qxl;
        geta = breta * gres[i] * reta;
        ginel = (1.0 - eta[i] - breta) * gres[i] * r2pion * qx2l;
        
        gtotal = gpion + geta * wthr + ginel;
        if (i == 1) gtotal = gpion / Math.sqrt(s / m22) + geta * wthr + ginel;
        
        qxgamma = Math.pow((x2 + kr22) / (x2 + k22), lp);
        ggamma  = gres[i] * Math.pow(rv[0], 2*lp+1) * qxgamma;
        
        den = Math.pow(m22 - s, 2) + Math.pow(mres[i] * gtotal, 2);
        
        breitw = mres[i] * Math.sqrt(gpion * ggamma / (eta[i] * rv[0] * rpion)) / den; 
        
        f1 = kkr * mn * eta[i];
        f2 = (2*jr+1) * pi * qr * mres[i] * gres[i];
        
        if (i != 1) a = Math.sqrt(f1 / f2) / 50.6;
        
        double[] breit = new double[2];
        breit[0] = a * breitw * (m22 - s);
        breit[1] = a * breitw * mres[i] * gtotal;
        
        return breit;    	
    }

    private void background(double s, double q2) {
    	
        double[][] exte_r1 = new double[2][9];
        double[][] exte_b1 = new double[2][9];
        double[][] extm_r1 = new double[2][9];
        double[][] extm_b1 = new double[2][9];
        double[][] exts_r1 = new double[2][9];
        double[][] exts_b1 = new double[2][9];
        
        double[] exte_r3 = new double[8];
        double[] exte_b3 = new double[8];
        double[] extm_r3 = new double[8];
        double[] extm_b3 = new double[8];
        double[] exts_r3 = new double[8];
        double[] exts_b3 = new double[8];
        
    	
        double[][]  ae1 = new double[2][9];
        double[][] cre1 = new double[2][9];
        double[][] cbe1 = new double[2][9];
        double[][]  am1 = new double[2][9];
        double[][] crm1 = new double[2][9];
        double[][] cbm1 = new double[2][9];
        
        double[]  ae3 = new double[8];
        double[] cre3 = new double[8];
        double[] cbe3 = new double[8];
        double[]  am3 = new double[8];
        double[] crm3 = new double[8];
        double[] cbm3 = new double[8];
        
        double[]   as3 = new double[8];
        double[][] as1 = new double[2][9];
        
    	
        backBorn(s, q2);
        backRegge(s, q2);
        
        // Initialize transverse amplitudes
        // M1-
        am1[0][0] = 1.0;
        crm1[0][0] = 0.0;
        cbm1[0][0] = 0.0;
        am3[7] = 0.05;
        crm3[7] = 0.0;
        cbm3[7] = 0.0;
        
        // E0+
        ae1[0][1] = 0.05;
        cre1[0][1] = 0.0;
        cbe1[0][1] = 0.0;
        ae3[4] = 0.05;
        cre3[4] = 0.0;
        cbe3[4] = 0.0;
        
        // M1+
        am1[0][6] = 0.05;
        crm1[0][6] = 0.0;
        cbm1[0][6] = 0.0;
        am3[0] = 1.0;
        crm3[0] = 0.0;
        cbm3[0] = 0.0;
        
        // E1+
        ae1[0][6] = 0.05;
        cre1[0][6] = 0.0;
        cbe1[0][6] = 0.0;
        ae3[0] = 1.0;
        cre3[0] = 0.0;
        cbe3[0] = 0.0;
        
        // M2-
        am1[0][2] = 1.0;
        crm1[0][2] = 3.0;
        cbm1[0][2] = 1.0;
        am3[1] = 1.0;
        crm3[1] = 0.0;
        cbm3[1] = 0.0;
        
        // E2-
        ae1[0][2] = 1.0;
        cre1[0][2] = 0.0;
        cbe1[0][2] = -0.9;
        ae3[1] = 1.0;
        cre3[1] = 0.0;
        cbe3[1] = 1.0;
        
        // M2+
        am1[0][8] = 1.0;
        crm1[0][8] = 0.0;
        cbm1[0][8] = 0.0;
        am3[6] = 1.0;
        crm3[6] = 0.0;
        cbm3[6] = 0.0;
        
        // E2+
        ae1[0][8] = 1.0;
        cre1[0][8] = 0.0;
        cbe1[0][8] = 0.0;
        ae3[6] = 1.0;
        cre3[6] = 0.0;
        cbe3[6] = 0.0;
        
        // M3-
        am1[0][5] = 1.0;
        crm1[0][5] = 3.0;
        cbm1[0][5] = 1.0;
        am3[5] = 1.0;
        crm3[5] = 0.0;
        cbm3[5] = 0.0;
        
        // E3-
        ae1[0][5] = 1.0;
        cre1[0][5] = 2.0;
        cbe1[0][5] = 1.0;
        ae3[5] = 1.0;
        cre3[5] = 1.0;
        cbe3[5] = 0.0;
        
        // M3+
        am1[0][7] = 1.0;
        crm1[0][7] = 0.0;
        cbm1[0][7] = 0.0;
        am3[2] = 1.0;
        crm3[2] = 2.0;
        cbm3[2] = 1.0;
        
        // E3+
        ae1[0][7] = 1.0;
        cre1[0][7] = 0.0;
        cbe1[0][7] = 0.0;
        ae3[2] = 1.0;
        cre3[2] = 0.0;
        cbe3[2] = 0.0;

        // Calculate transverse amplitudes
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 9; j++) {
                exte_r1[i][j] = ae1[i][j] * Math.pow(s / 4, cre1[i][j])     * Math.pow(s - 1.16, 2);
                exte_r1[i][j] = exte_r1[i][j] / (ae1[i][j]                  * Math.pow(s - 1.16, 2) + 1.0);
                exte_b1[i][j] = Math.pow(1.16 / s, cbe1[i][j]) / (ae1[i][j] * Math.pow(s - 1.16, 2) + 1.0);
                
                extm_r1[i][j] = am1[i][j] * Math.pow(s / 4, crm1[i][j])     * Math.pow(s - 1.16, 2);
                extm_r1[i][j] = extm_r1[i][j] / (am1[i][j]                  * Math.pow(s - 1.16, 2) + 1.0);
                extm_b1[i][j] = Math.pow(1.16 / s, cbm1[i][j]) / (am1[i][j] * Math.pow(s - 1.16, 2) + 1.0);
            }
        }

        for (int j = 0; j < 8; j++) {
            exte_r3[j] = ae3[j] * Math.pow(s / 4, cre3[j])     * Math.pow(s - 1.16, 2);
            exte_r3[j] = exte_r3[j] / (ae3[j]                  * Math.pow(s - 1.16, 2) + 1.0);
            exte_b3[j] = Math.pow(1.16 / s, cbe3[j]) / (ae3[j] * Math.pow(s - 1.16, 2) + 1.0);
            
            extm_r3[j] = am3[j] * Math.pow(s / 4, crm3[j])     * Math.pow(s - 1.16, 2);
            extm_r3[j] = extm_r3[j] / (am3[j]                  * Math.pow(s - 1.16, 2) + 1.0);
            extm_b3[j] = Math.pow(1.16 / s, cbm3[j]) / (am3[j] * Math.pow(s - 1.16, 2) + 1.0);
        }

        // Combine transverse amplitudes
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 9; j++) {
                BM1[i][j] = extm_r1[i][j] * RM1[i][j] + extm_b1[i][j] * BbM1[i][j];
                BE1[i][j] = exte_r1[i][j] * RE1[i][j] + exte_b1[i][j] * BbE1[i][j];
            }
        }

        for (int j = 0; j < 8; j++) {
            BM3[j] = extm_r3[j] * RM3[j] + extm_b3[j] * BbM3[j];
            BE3[j] = exte_r3[j] * RE3[j] + exte_b3[j] * BbE3[j];
        }

        // Initialize longitudinal amplitudes
        // S1-
        as1[0][0] = 0.2407;
        as3[7] = 0.5525;
        
        // S0+
        as1[0][1] = 0.0;
        as3[4] = 0.0;
        
        // S1+
        as1[0][6] = 0.0;
        as3[0] = 0.0;
        
        // S2-
        as1[0][2] = 1.0;
        as3[1] = 0.0;
        
        // S2+
        as1[0][8] = 0.762;
        as3[6] = 0.2254;
        
        // S3-
        as1[0][5] = 1.0;
        as3[5] = 0.63688;
        
        // S3+
        as1[0][7] = 1.0;
        as3[2] = 1.0;

        // Calculate longitudinal amplitudes
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 9; j++) {
                exts_r1[i][j] =        as1[i][j] * Math.pow(s - 1.16, 2) / (as1[i][j] * Math.pow(s - 1.16, 2) + 1.0);
                exts_b1[i][j] = 1.0 / (as1[i][j] * Math.pow(s - 1.16, 2)                                      + 1.0);
                BS1[i][j] = exts_r1[i][j] * RS1[i][j] + exts_b1[i][j] * BbS1[i][j];
            }
        }

        for (int j = 0; j < 8; j++) {
            exts_r3[j] =        as3[j] * Math.pow(s - 1.16, 2) /           (as3[j]    * Math.pow(s - 1.16, 2) + 1.0);
            exts_b3[j] = 1.0 / (as3[j] * Math.pow(s - 1.16, 2)                                                + 1.0);
            BS3[j] = exts_r3[j] * RS3[j] + exts_b3[j] * BbS3[j];
        }    	
    }
    
    private void backBorn(double s, double q2) {
    	
    	int kpower=8,jw,minp,maxp,mwp=1;
    	
        double[] tmult = new double[8];
        double[]    wx = new double[8];

        double[][] bem = new double[4][3];
        double[][] bep = new double[4][3];
        double[][] bmm = new double[4][3];
        double[][] bmp = new double[4][3];
        double[][] bsm = new double[4][3];
        double[][] bsp = new double[4][3];
    	
        double wpoint = Math.sqrt(s);

        for (int iw = 0; iw < Maxmpoints; iw++) {
            if (wpoint < wtab[iw] && mwp == 1) {
                mwp = iw;
                break;
            }
        }

        minp = mwp - 3;
        
        if (minp <= 0) {
            minp = 1;
            maxp = 8;
        } else {
            maxp = mwp + 4;
            if (maxp > Maxmpoints) {
                maxp = Maxmpoints;
                minp = Maxmpoints - 7;
            }
        }
        
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
            	
                // E"l"m
                jw = 0;
                
                for (int iw = minp - 1; iw < maxp; iw++) {
                    tmult[jw] = aem[iw][i][j];
                    wx[jw] = wtab[iw];
                    jw++;
                }                
                bem[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // E"l"p
                jw = 0;
                for (int iw = minp - 1; iw < maxp; iw++) {
                    tmult[jw] = aep[iw][i][j];
                    jw++;
                } 
                bep[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // M"l"m
                jw = 0;
                for (int iw = minp - 1; iw < maxp; iw++) {
                    tmult[jw] = amm[iw][i][j];
                    jw++;
                } 
                bmm[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // M"l"p
                jw = 0;
                for (int iw = minp - 1; iw < maxp; iw++) {
                    tmult[jw] = amp[iw][i][j];
                    jw++;
                } 
                bmp[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // S"l"m
                jw = 0;
                for (int iw = minp - 1; iw < maxp; iw++) {
                    tmult[jw] = asm[iw][i][j];
                    jw++;
                }
                bsm[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // S"l"p
                jw = 0;
                for (int iw = minp - 1; iw < maxp; iw++) {
                    tmult[jw] = asp[iw][i][j];
                    jw++;
                }
                bsp[i][j] = p.polint(tmult, wx, kpower, wpoint);

            }
        }

        mixingAmpl(s, q2);

        // Background multipole amplitudes with I=3/2        
        BbM3[0] = bmp[1][2];
        BbM3[1] = bmm[2][2];
        BbM3[2] = bmp[3][2];
        BbM3[3] = 0.0;
        BbM3[4] = 0.0;
        BbM3[5] = bmm[3][2];
        BbM3[6] = bmp[2][2];
        BbM3[7] = bmm[1][2] + M1[0];

        BbE3[0] = bep[1][2];
        BbE3[1] = bem[2][2];
        BbE3[2] = bep[3][2];
        BbE3[3] = 0.0;
        BbE3[4] = bep[0][2] + E0[0];
        BbE3[5] = bem[3][2];
        BbE3[6] = bep[2][2];
        BbE3[7] = bem[1][2];

        BbS3[0] = bsp[1][2];
        BbS3[1] = bsm[2][2];
        BbS3[2] = bsp[3][2];
        BbS3[3] = 0.0;
        BbS3[4] = bsp[0][2] + S0[0];
        BbS3[5] = bsm[3][2];
        BbS3[6] = bsp[2][2];
        BbS3[7] = bsm[1][2] + S1[0];

        // Background multipole amplitudes with I=1/2 for combinations (1/2) and (0)        
        BbM1[0][0] = bmm[1][1] + M1[0];
        BbM1[0][1] = 0.0;
        BbM1[0][2] = bmm[2][1];
        BbM1[0][3] = 0.0;
        BbM1[0][4] = 0.0;
        BbM1[0][5] = bmm[3][1];
        BbM1[0][6] = bmp[1][1];
        BbM1[0][7] = bmp[3][1];
        BbM1[0][8] = bmp[2][1];

        BbM1[1][0] = bmm[1][0] + M1[1];
        BbM1[1][1] = 0.0;
        BbM1[1][2] = bmm[2][0];
        BbM1[1][3] = 0.0;
        BbM1[1][4] = 0.0;
        BbM1[1][5] = bmm[3][0];
        BbM1[1][6] = bmp[1][0];
        BbM1[1][7] = bmp[3][0];
        BbM1[1][8] = bmp[2][0];

        BbE1[0][0] = bem[1][1];
        BbE1[0][1] = bep[0][1] + E0[0];
        BbE1[0][2] = bem[2][1];
        BbE1[0][3] = 0.0;
        BbE1[0][4] = 0.0;
        BbE1[0][5] = bem[3][1];
        BbE1[0][6] = bep[1][1];
        BbE1[0][7] = bep[3][1];
        BbE1[0][8] = bep[2][1];

        BbE1[1][0] = bem[1][0];
        BbE1[1][1] = bep[0][0] + E0[1];
        BbE1[1][2] = bem[2][0];
        BbE1[1][3] = 0.0;
        BbE1[1][4] = 0.0;
        BbE1[1][5] = bem[3][0];
        BbE1[1][6] = bep[1][0];
        BbE1[1][7] = bep[3][0];
        BbE1[1][8] = bep[2][0];

        BbS1[0][0] = bsm[1][1] + S1[0];
        BbS1[0][1] = bsp[0][1] + S0[0];
        BbS1[0][2] = bsm[2][1];
        BbS1[0][3] = 0.0;
        BbS1[0][4] = 0.0;
        BbS1[0][5] = bsm[3][1];
        BbS1[0][6] = bsp[1][1];
        BbS1[0][7] = bsp[3][1];
        BbS1[0][8] = bsp[2][1];

        BbS1[1][0] = bsm[1][0] + S1[1];
        BbS1[1][1] = bsp[0][0] + S0[1];
        BbS1[1][2] = bsm[2][0];
        BbS1[1][3] = 0.0;
        BbS1[1][4] = 0.0;
        BbS1[1][5] = bsm[3][0];
        BbS1[1][6] = bsp[1][0];
        BbS1[1][7] = bsp[3][0];
        BbS1[1][8] = bsp[2][0];

        // Background multipole amplitudes with I=1/2 for combinations p(1/2) and n(1/2)        
        for (int j = 0; j < 9; j++) {
            BbM1[0][j] = BbM1[0][j] / 3.0 + BbM1[1][j];
            BbM1[1][j] = 2.0 * BbM1[1][j] - BbM1[0][j];

            BbE1[0][j] = BbE1[0][j] / 3.0 + BbE1[1][j];
            BbE1[1][j] = 2.0 * BbE1[1][j] - BbE1[0][j];

            BbS1[0][j] = BbS1[0][j] / 3.0 + BbS1[1][j];
            BbS1[1][j] = 2.0 * BbS1[1][j] - BbS1[0][j];
        }
   	
    }
       
    public void mixingAmpl(double s, double q2) {

        double W, on1, on2, op1, om1, op2, om2;
        double Lpv, Lmix;
        double[] B1 = new double[2];

        W = Math.sqrt(s);
        on1 = W - Egamma;
        on2 = W - Epion;
        op1 = on1 + mn;
        om1 = on1 - mn;
        op2 = on2 + mn;
        om2 = on2 - mn;
        
        Lmix = 0.55 + cspa;
        Lpv = Lmix * Lmix / (Lmix * Lmix + Epion * Epion - mp22);
        B1[0] = Lpv * 3.134 * (F2p - F2n) / mn / W / 2.0;
        B1[1] = Lpv * 3.134 * (F2p + F2n) / mn / W / 2.0;

        // Loop for (+) and (0) amplitudes. (-) = 0.
        for (int i = 0; i < 2; i++) {
            E0[i] =  Math.sqrt(op1 * op2) * B1[i] * (W - mn);
            M1[i] = -Math.sqrt(om1 * om2) * B1[i] * (W + mn);
            S0[i] =  Math.sqrt(om1 * op2) * B1[i] * op1;
            S1[i] = -Math.sqrt(op1 * om2) * B1[i] * om1;
        }
    }
    
    public void backRegge(double s, double q2) {
    	
        int kpower=8,jw,minp,maxp,mwp = 1;  
        
        double[] tmult = new double[8];
        double[] wx = new double[8];

        double[][] bem = new double[4][3];
        double[][] bep = new double[4][3];
        double[][] bmm = new double[4][3];
        double[][] bmp = new double[4][3];
        double[][] bsm = new double[4][3];
        double[][] bsp = new double[4][3];
        
        double wpoint = Math.sqrt(s);
        
        for (int iw = 0; iw < Maxmpoints; iw++) {
            if (wpoint < wtab[iw] && mwp == 1) {
                mwp = iw;
                break;
            }
        }
        
        minp = mwp - 3;
        
        if (minp <= 0) {
            minp = 1;
            maxp = 8;
        } else {
            maxp = mwp + 4;
            if (maxp > Maxmpoints) {
                maxp = Maxmpoints;
                minp = Maxmpoints - 7;
            }
        }

        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                // E"l"m            	
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = raem[iw][i][j];
                    wx[jw] = wtab[iw];
                    jw++;
                }
                bem[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // E"l"p
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = raep[iw][i][j];
                    jw++;
                }                
                bep[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // M"l"m
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = ramm[iw][i][j];
                    jw++;
                }                
                bmm[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // M"l"p
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = ramp[iw][i][j];
                    jw++;
                }                
                bmp[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // S"l"m
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = rasm[iw][i][j];
                    jw++;
                }                
                bsm[i][j] = p.polint(tmult, wx, kpower, wpoint);

                // S"l"p
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = rasp[iw][i][j];
                    jw++;
                }               
                bsp[i][j] = p.polint(tmult, wx, kpower, wpoint);
            }
        }

        // Background multipole amplitudes with I=3/2
        RM3[0] = bmp[1][2];
        RM3[1] = bmm[2][2];
        RM3[2] = bmp[3][2];
        RM3[3] = 0.0;
        RM3[4] = 0.0;
        RM3[5] = bmm[3][2];
        RM3[6] = bmp[2][2];
        RM3[7] = bmm[1][2];

        RE3[0] = bep[1][2];
        RE3[1] = bem[2][2];
        RE3[2] = bep[3][2];
        RE3[3] = 0.0;
        RE3[4] = bep[0][2];
        RE3[5] = bem[3][2];
        RE3[6] = bep[2][2];
        RE3[7] = bem[1][2];

        RS3[0] = bsp[1][2];
        RS3[1] = bsm[2][2];
        RS3[2] = bsp[3][2];
        RS3[3] = 0.0;
        RS3[4] = bsp[0][2];
        RS3[5] = bsm[3][2];
        RS3[6] = bsp[2][2];
        RS3[7] = bsm[1][2];

        // Background multipole amplitudes with I=1/2 for combinations (1/2) and (0)
        RM1[0][0] = bmm[1][1];
        RM1[0][1] = 0.0;
        RM1[0][2] = bmm[2][1];
        RM1[0][3] = 0.0;
        RM1[0][4] = 0.0;
        RM1[0][5] = bmm[3][1];
        RM1[0][6] = bmp[1][1];
        RM1[0][7] = bmp[3][1];
        RM1[0][8] = bmp[2][1];

        RM1[1][0] = bmm[1][0];
        RM1[1][1] = 0.0;
        RM1[1][2] = bmm[2][0];
        RM1[1][3] = 0.0;
        RM1[1][4] = 0.0;
        RM1[1][5] = bmm[3][0];
        RM1[1][6] = bmp[1][0];
        RM1[1][7] = bmp[3][0];
        RM1[1][8] = bmp[2][0];

        RE1[0][0] = bem[1][1];
        RE1[0][1] = bep[0][1];
        RE1[0][2] = bem[2][1];
        RE1[0][3] = 0.0;
        RE1[0][4] = 0.0;
        RE1[0][5] = bem[3][1];
        RE1[0][6] = bep[1][1];
        RE1[0][7] = bep[3][1];
        RE1[0][8] = bep[2][1];

        RE1[1][0] = bem[1][0];
        RE1[1][1] = bep[0][0];
        RE1[1][2] = bem[2][0];
        RE1[1][3] = 0.0;
        RE1[1][4] = 0.0;
        RE1[1][5] = bem[3][0];
        RE1[1][6] = bep[1][0];
        RE1[1][7] = bep[3][0];
        RE1[1][8] = bep[2][0];

        RS1[0][0] = bsm[1][1];
        RS1[0][1] = bsp[0][1];
        RS1[0][2] = bsm[2][1];
        RS1[0][3] = 0.0;
        RS1[0][4] = 0.0;
        RS1[0][5] = bsm[3][1];
        RS1[0][6] = bsp[1][1];
        RS1[0][7] = bsp[3][1];
        RS1[0][8] = bsp[2][1];

        RS1[1][0] = bsm[1][0];
        RS1[1][1] = bsp[0][0];
        RS1[1][2] = bsm[2][0];
        RS1[1][3] = 0.0;
        RS1[1][4] = 0.0;
        RS1[1][5] = bsm[3][0];
        RS1[1][6] = bsp[1][0];
        RS1[1][7] = bsp[3][0];
        RS1[1][8] = bsp[2][0];

        // Background multipole amplitudes with I=1/2 for combinations p(1/2) and n(1/2)
        for (int j = 0; j < 9; j++) {
            RM1[0][j] = RM1[0][j] / 3.0 + RM1[1][j];
            RM1[1][j] = 2.0 * RM1[1][j] - RM1[0][j];
            RE1[0][j] = RE1[0][j] / 3.0 + RE1[1][j];
            RE1[1][j] = 2.0 * RE1[1][j] - RE1[0][j];
            RS1[0][j] = RS1[0][j] / 3.0 + RS1[1][j];
            RS1[1][j] = 2.0 * RS1[1][j] - RS1[0][j];
        }
    }
    
    public void highMult(double s, double q2, double costh) {
    	   	
        backHigh(s);

        // Polynomial coefficients
        double p1 = -52.5 * costh + 157.5 * Math.pow(costh, 3); // P''_5
        double p2 = 15.0 * costh;                               // P''_3
        double p3 = -7.5 + 52.5 * Math.pow(costh, 2);           // P''_4
        double p4 = 1.875 - 26.25 * Math.pow(costh, 2) + 39.375 * Math.pow(costh, 4); // P'_5
        double p5 = -1.5 + 7.5 * Math.pow(costh, 2);            // P'_3
        double p6 = -7.5 * costh + 17.5 * Math.pow(costh, 3);   // P'_4
        double p7 = 13.125 - 236.25 * Math.pow(costh, 2) + 433.125 * Math.pow(costh, 4); // P''_6
        double p8 = 13.125 * costh - 78.75 * Math.pow(costh, 3) + 86.625 * Math.pow(costh, 5); // P'_6

        // Loops for p(1/2), n(1/2), 3/2 
        for (int i = 0; i < 3; i++) {
            phih[0][i] =  (4.0 * M4p[i] + E4p[i]) * p4 + (5.0 * M5p[i] + E5p[i]) * p8;
            phih[0][i] += (5.0 * M4m[i] + E4m[i]) * p5 + (6.0 * M5m[i] + E5m[i]) * p6;
            phih[1][i] =  (5.0 * M4p[i] + 4.0 * M4m[i]) * p6 + (6.0 * M5p[i] + 5.0 * M5m[i]) * p4;
            phih[2][i] =  (E4p[i] - M4p[i]) * p1 + (E5p[i] - M5p[i]) * p7;
            phih[2][i] += (M4m[i] + E4m[i]) * p2 + (M5m[i] + E5m[i]) * p3;
            phih[3][i] =  (M4p[i] - E4p[i]        - M4m[i] - E4m[i]) * p3;
            phih[3][i] += (M5p[i] - E5p[i]        - M5m[i] - E5m[i]) * p1;
            phih[4][i] =  5.0 * S4p[i] * p4 - 4.0 * S4m[i] * p5 + 6.0 * S5p[i] * p8 - 5.0 * S5m[i] * p6;
            phih[5][i] = -5.0 * S4p[i] * p6 + 4.0 * S4m[i] * p6 - 6.0 * S5p[i] * p4 + 5.0 * S5m[i] * p4;
        }

        // Construct phih for reactions pi0, pi+
        for (int i = 0; i < 6; i++) {
            phih[i][0] = phih[i][0] + phih[i][2] / 1.5;
            phih[i][1] = Math.sqrt(2.0) * (phih[i][0] - phih[i][2]);
        }
    }    	
    
    public void backHigh(double s) {
    	
    	int kpower=8,jw,minp,maxp,mwp=1;
    	
        double[] tmult = new double[kpower];
        double[] wx = new double[kpower];
    	
        double[][] bem = new double[2][3];
        double[][] bep = new double[2][3];
        double[][] bmm = new double[2][3];
        double[][] bmp = new double[2][3];
        double[][] bsm = new double[2][3];
        double[][] bsp = new double[2][3];

        double wpoint = Math.sqrt(s);

        for (int iw = 0; iw < Maxmpoints; iw++) {
            if (wpoint < wtab[iw] && mwp == 1) {
                mwp = iw;
                break;
            }
        }
        
        minp = mwp - 3;

        if (minp <= 0) {
            minp = 1;
            maxp = 8;
        } else {
            maxp = mwp + 4;
            if (maxp > Maxmpoints) {
                maxp = Maxmpoints;
                minp = Maxmpoints - 7;
            }
        }

        for (int i = 4; i < 6; i++) {
            for (int j = 0; j < 3; j++) {
                // E"l"m
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = aem[iw][i][j];
                    wx[jw] = wtab[iw];
                    jw++;
                }
                bem[i - 4][j] = p.polint(tmult, wx, kpower, wpoint);

                // E"l"p
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = aep[iw][i][j];
                    jw++;
                }
                bep[i - 4][j] = p.polint(tmult, wx, kpower, wpoint);

                // M"l"m
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = amm[iw][i][j];
                    jw++;
                }
                bmm[i - 4][j] = p.polint(tmult, wx, kpower, wpoint);

                // M"l"p
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = amp[iw][i][j];
                    jw++;
                }
                bmp[i - 4][j] = p.polint(tmult, wx, kpower, wpoint);

                // S"l"m
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = asm[iw][i][j];
                    jw++;
                }
                bsm[i - 4][j] = p.polint(tmult, wx, kpower, wpoint);

                // S"l"p
                jw = 0;
                for (int iw = minp-1; iw < maxp; iw++) {
                    tmult[jw] = asp[iw][i][j];
                    jw++;
                }
                bsp[i - 4][j] = p.polint(tmult, wx, kpower, wpoint);
            }
        }

        // Construct background multipole amplitudes with I=1/2
        E4m[0] =  bem[0][1] / 3.0 + bem[0][0];
        E4m[1] = -bem[0][1] / 3.0 + bem[0][0];
        E4p[0] =  bep[0][1] / 3.0 + bep[0][0];
        E4p[1] = -bep[0][1] / 3.0 + bep[0][0];
        M4m[0] =  bmm[0][1] / 3.0 + bmm[0][0];
        M4m[1] = -bmm[0][1] / 3.0 + bmm[0][0];
        M4p[0] =  bmp[0][1] / 3.0 + bmp[0][0];
        M4p[1] = -bmp[0][1] / 3.0 + bmp[0][0];
        S4m[0] =  bsm[0][1] / 3.0 + bsm[0][0];
        S4m[1] = -bsm[0][1] / 3.0 + bsm[0][0];
        S4p[0] =  bsp[0][1] / 3.0 + bsp[0][0];
        S4p[1] = -bsp[0][1] / 3.0 + bsp[0][0];
        E5m[0] =  bem[1][1] / 3.0 + bem[1][0];
        E5m[1] = -bem[1][1] / 3.0 + bem[1][0];
        E5p[0] =  bep[1][1] / 3.0 + bep[1][0];
        E5p[1] = -bep[1][1] / 3.0 + bep[1][0];
        M5m[0] =  bmm[1][1] / 3.0 + bmm[1][0];
        M5m[1] = -bmm[1][1] / 3.0 + bmm[1][0];
        M5p[0] =  bmp[1][1] / 3.0 + bmp[1][0];
        M5p[1] = -bmp[1][1] / 3.0 + bmp[1][0];
        S5m[0] =  bsm[1][1] / 3.0 + bsm[1][0];
        S5m[1] = -bsm[1][1] / 3.0 + bsm[1][0];
        S5p[0] =  bsp[1][1] / 3.0 + bsp[1][0];
        S5p[1] = -bsp[1][1] / 3.0 + bsp[1][0];

        // Construct background multipole amplitudes with I=3/2
        E4m[2] = bem[0][2];
        E5m[2] = bem[1][2];
        E4p[2] = bep[0][2];
        E5p[2] = bep[1][2];
        M4m[2] = bmm[0][2];
        M5m[2] = bmm[1][2];
        M4p[2] = bmp[0][2];
        M5p[2] = bmp[1][2];
        S4m[2] = bsm[0][2];
        S5m[2] = bsm[1][2];
        S4p[2] = bsp[0][2];
        S5p[2] = bsp[1][2];
    } 
    
    public void helAmpl(double q2, double theta, double kgamma) {
    	
        double teta, qq, sx, cx, sx2;
        double[][][] PHreac = new double[6][3][2];

        teta = theta * pi / 180.0;
        qq = Math.sqrt(q2) / kgamma;
        sx = Math.sin(teta);
        cx = Math.cos(teta);
        sx2 = sx * sx;

        // Intermediate amplitudes for different reactions
        for (int i = 0; i < 6; i++) {  
            for (int j = 0; j < 2; j++) { // Re(0) Im(1)
                // Reaction g + p => p + pi0
                PHreac[i][0][j] = phi_amp[i][j][0];  

                // Reaction g + p => n + pi+
                PHreac[i][1][j] = phi_amp[i][j][1];  
                // PHreac[i][2][j] is unused
            }
        }

        // Helicity Amplitudes
        for (int i = 0; i < 2; i++) { // p+pi0, n+pi+
            for (int j = 0; j < 2; j++) { // Re(0) Im(1)
                // h1 = f+ + - =  f- - +
                H1[i][j] = -sx * (PHreac[2][i][j] + cx * PHreac[3][i][j]) / a2;

                // h2 = f- + + = -f+ - -
                H2[i][j] = 2.0 * PHreac[0][i][j] - 2.0 * cx * PHreac[1][i][j] + sx2 * PHreac[3][i][j];
                H2[i][j] = -H2[i][j] / a2;

                // h3 = f+ - + = f- + -
                H3[i][j] = -sx2 * PHreac[3][i][j] / a2;

                // h4 = f+ + + = f- - -
                H4[i][j] = 2.0 * PHreac[1][i][j] + PHreac[2][i][j] + cx * PHreac[3][i][j];
                H4[i][j] = sx * H4[i][j] / a2;

                // h5 = f- - 0 = -f+ + 0
                H5[i][j] = qq * (PHreac[4][i][j] + cx * PHreac[5][i][j]);

                // h6 = f+ - 0 =  f- + 0
                H6[i][j] = qq * sx * PHreac[5][i][j];
            }
        }
    }  
    
    public void responseFunctions(double s, double azimut, double epsilon) {
    	
        double mnuc = 0.93827;
        double[] mpion = new double[2];
        mpion[0] = 0.13498;
        mpion[1] = 0.13957;
        
        // Loop over pi-zero and pi+
        for (int i = 0; i < 2; i++) {

            double s_mnuc2 = s - mnuc * mnuc;
            double mpion2 = mpion[i] * mpion[i];

            double term1 = (s - mnuc * mnuc + mpion2) / 2.0;
            double qpi = Math.sqrt((term1 * term1) / s - mpion2);
            double qkratio = qpi / (s_mnuc2 / 2.0 / Math.sqrt(s));

            sigmaT[i] = 0.0;
            sigmaL[i] = 0.0;
            sigmaTT[i] = 0.0;
            sigmaTL[i] = 0.0;
            sigmaTLP[i] = 0.0;
            sigmaTy[i] = 0.0;
            sigmaLy[i] = 0.0;
            sigmaTTx[i] = 0.0;
            sigmaTTy[i] = 0.0;
            sigmaTTz[i] = 0.0;
            sigmaTLx[i] = 0.0;
            sigmaTLy[i] = 0.0;
            sigmaTLz[i] = 0.0;
            sigmaTTPx[i] = 0.0;
            sigmaTTPz[i] = 0.0;
            sigmaTLPx[i] = 0.0;
            sigmaTLPy[i] = 0.0;
            sigmaTLPz[i] = 0.0;

            // Loop over j (real and imaginary parts)
            for (int j = 0; j < 2; j++) {
                sigmaT[i] += Math.pow(H1[i][j], 2) + Math.pow(H2[i][j], 2) +
                             Math.pow(H3[i][j], 2) + Math.pow(H4[i][j], 2);
                sigmaL[i] += Math.pow(H5[i][j], 2) + Math.pow(H6[i][j], 2);
                sigmaTT[i] -= H4[i][j] * H1[i][j] - H3[i][j] * H2[i][j];

                sigmaTL[i] += H5[i][j] * (H1[i][j] - H4[i][j]) +
                              H6[i][j] * (H2[i][j] + H3[i][j]);

                sigmaTLPx[i] += H5[i][j] * (-H3[i][j] + H2[i][j]) +
                                H6[i][j] * (H1[i][j] + H4[i][j]);
                sigmaTLPy[i] -= H5[i][j] * (H3[i][j] + H2[i][j]) +
                                H6[i][j] * (H1[i][j] - H4[i][j]);
                sigmaTLPz[i] += H5[i][j] * (H1[i][j] + H4[i][j]) +
                                H6[i][j] * (H3[i][j] - H2[i][j]);

                sigmaTTPx[i] += H1[i][j] * H2[i][j] + H3[i][j] * H4[i][j];
                sigmaTTPz[i] += Math.pow(H1[i][j], 2) - Math.pow(H2[i][j], 2) +
                                Math.pow(H3[i][j], 2) - Math.pow(H4[i][j], 2);
            }

            // Cross terms between real and imaginary parts
            sigmaTLP[i] += H5[i][1] * (H1[i][0] - H4[i][0]) - H5[i][0] * (H1[i][1] - H4[i][1]) +
                           H6[i][1] * (H2[i][0] + H3[i][0]) - H6[i][0] * (H2[i][1] + H3[i][1]);

            sigmaTy[i] += H1[i][1] * H2[i][0] - H1[i][0] * H2[i][1] +
                          H3[i][1] * H4[i][0] - H3[i][0] * H4[i][1];

            sigmaLy[i] += H5[i][1] * H6[i][0] - H5[i][0] * H6[i][1];

            sigmaTLx[i] += H5[i][0] * (H3[i][1] - H2[i][1]) - H5[i][1] * (H3[i][0] - H2[i][0]) -
                           H6[i][0] * (H1[i][1] + H4[i][1]) + H6[i][1] * (H1[i][0] + H4[i][0]);

            sigmaTLy[i] -= H5[i][0] * (H3[i][1] + H2[i][1]) - (-H5[i][1] * (H3[i][0] + H2[i][0])) +
                           H6[i][0] * (H1[i][1] - H4[i][1]) -   H6[i][1] * (H1[i][0] - H4[i][0]);

            sigmaTLz[i] -= H5[i][0] * (H1[i][1] + H4[i][1]) - (-H5[i][1] * (H1[i][0] + H4[i][0])) +
                           H6[i][0] * (H2[i][1] - H3[i][1]) -   H6[i][1] * (H2[i][0] - H3[i][0]);

            sigmaTTx[i] += H1[i][0] * H3[i][1] - H1[i][1] * H3[i][0] -
                           H2[i][0] * H4[i][1] + H2[i][1] * H4[i][0];

            sigmaTTy[i] -= H1[i][0] * H3[i][1] - (-H1[i][1] * H3[i][0]) -
                           H2[i][0] * H4[i][1] +   H2[i][1] * H4[i][0];

            sigmaTTz[i] -= H1[i][0] * H4[i][1] - (-H1[i][1] * H4[i][0]) -
                           H2[i][0] * H3[i][1] +   H2[i][1] * H3[i][0];

            // Now apply the scaling factors
            sigmaT[i]    *= qkratio / 2.0;
            sigmaL[i]    *= qkratio;
            sigmaTT[i]   *= qkratio;
            sigmaTL[i]   *= qkratio / Math.sqrt(2.0);
            sigmaTLP[i]  *= qkratio / Math.sqrt(2.0);
            sigmaTy[i]   *= qkratio;
            sigmaLy[i]   *= qkratio * 2.0;
            sigmaTLx[i]  *= qkratio / Math.sqrt(2.0);
            sigmaTLy[i]  *= qkratio / Math.sqrt(2.0);
            sigmaTLz[i]  *= qkratio / Math.sqrt(2.0);
            sigmaTTx[i]  *= qkratio;
            sigmaTTy[i]  *= qkratio;
            sigmaTTz[i]  *= qkratio;
            sigmaTLPx[i] *= qkratio / Math.sqrt(2.0);
            sigmaTLPy[i] *= qkratio / Math.sqrt(2.0);
            sigmaTLPz[i] *= qkratio / Math.sqrt(2.0);
            sigmaTTPx[i] *= qkratio;
            sigmaTTPz[i] *= qkratio / 2.0;

            // Final calculation
            sigma[i] = sigmaT[i]
                     + epsilon * sigmaL[i]
                     + epsilon * sigmaTT[i] * Math.cos(2.0 * azimut)
                     + Math.sqrt(2.0 * epsilon * (1.0 + epsilon)) * sigmaTL[i] * Math.cos(azimut);
        }

    } 
    
    public void crossSections(double s, double q2, double phi, double eps) {

   	
    	double vl, vtt, vlt, vltp;
        double qkratio, qpi, sinth;
        double sigma_tt, sigma_lt, sigma_ltp;
            	
        double    mnuc = 0.93827;
        double[] mpion = {0.13498, 0.13957};
        
        vl   = eps;
        vtt  = eps;
        vlt  = Math.sqrt(2 * eps * (1 + eps));
        vltp = Math.sqrt(2 * eps * (1 - eps));
        
        boolean stop = false;
        
        float fw = (float) Math.sqrt(s);

        for (int i = 0; i < 2; i++) { //i=0: p-pi0 i=1: n-pi+
            qpi = Math.sqrt(Math.pow((s - mnuc * mnuc + mpion[i] * mpion[i]) / 2.0, 2) / s - mpion[i] * mpion[i]);
            qkratio = qpi / ((s - mnuc * mnuc) / 2. / Math.sqrt(s));

            sigmaT[i] = 0.0;
            sigmaL[i] = 0.0;
            sigmaTT[i] = 0.0;
            sigmaTL[i] = 0.0;
            sigmaTLP[i] = 0.0;

            for (int j = 0; j < 2; j++) { // Loop over j (real and imaginary parts)
                sigmaT[i]  = sigmaT[i] + Math.pow(H1[i][j], 2) + Math.pow(H2[i][j], 2)
                                       + Math.pow(H3[i][j], 2) + Math.pow(H4[i][j], 2);

                sigmaL[i]  = sigmaL[i] + Math.pow(H5[i][j], 2) + Math.pow(H6[i][j], 2);

                sigmaTT[i] = sigmaTT[i] - H4[i][j] * H1[i][j] + H3[i][j] * H2[i][j];

                sigmaTL[i] = sigmaTL[i] + H5[i][j] * (H1[i][j] - H4[i][j])
                                        + H6[i][j] * (H2[i][j] + H3[i][j]);
            }
//            if(fw>=1.11&&fw<=1.12) System.out.println(Math.toDegrees(phi)+" "+Math.sqrt(s));

            // The following assumes that the second part is outside the above loop
            // and that the arrays H1, H4, H5, H2, H3, H6 have at least two columns (j=0,1)
            
//            if(Math.toDegrees(phi)< 7.6) System.out.println(i+" "+Math.sqrt(s)+" "+(float)sigmaT[i]+" "+(float)sigmaL[i]+" "+(float)sigmaTT[i]+" "+(float)sigmaTL[i]);

            sigmaTLP[i] = sigmaTLP[i]
                + H5[i][1] * (H1[i][0] - H4[i][0]) - H5[i][0] * (H1[i][1] - H4[i][1])
                + H6[i][1] * (H2[i][0] + H3[i][0]) - H6[i][0] * (H2[i][1] + H3[i][1]);

            sigmaT[i]   = sigmaT[i]   * qkratio / 2.;
            sigmaL[i]   = sigmaL[i]   * qkratio;
            sigmaTT[i]  = sigmaTT[i]  * qkratio;
            sigmaTL[i]  = sigmaTL[i]  * qkratio / Math.sqrt(2.);
            sigmaTLP[i] = sigmaTLP[i] * qkratio / Math.sqrt(2.);

            sigma_tt  = vtt  * sigmaTT[i]  * Math.cos(2 * phi);
            sigma_ltp = vltp * sigmaTLP[i] * Math.sin(phi);
            sigma_lt  = vlt  * sigmaTL[i]  * Math.cos(phi);
            
            sigma[i] = sigmaT[i]
                     + vl * sigmaL[i]
                     + sigma_tt
                     + sigma_lt;

            if (q2 > 0)  robs[0][i] = sigma[i];
            if (q2 == 0) robs[0][i] = sigmaT[i];

            robs[1][i] = sigma_ltp / sigma[i];
            robs[2][i] = sigma_lt  / sigma[i];
            robs[3][i] = -sigmaTT[i] / sigmaT[i];
            robs[4][i] = sigmaT[i];
            robs[5][i] = sigmaT[i] + vl * sigmaL[i];
            robs[6][i] = sigmaTT[i];
            robs[7][i] = sigmaTL[i];
            robs[8][i] = sigmaTLP[i];
        }    	
    }
	
/*
    // Assume polint is defined as:
    // Given y[], x[], n, x0, returns y0 via a_point[0]
    public static double polint(double[] y, double[] x, int n, double x0) {
    	double a_point;
        // Implementation required!
        // For now, simple Lagrange interpolation
        double y0 = 0.0;
        for (int i = 0; i < n; i++) {
            double term = y[i];
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    term *= (x0 - x[j]) / (x[i] - x[j]);
                }
            }
            y0 += term;
        }
        return y0;
    }
*/
    
}
