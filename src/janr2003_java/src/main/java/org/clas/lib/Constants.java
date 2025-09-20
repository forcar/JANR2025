package org.clas.lib;

import java.util.ArrayList;
import java.util.List;

public class Constants {
	
	    public static Boolean debug=false;
        //String janrPath = FileSystemView.getFileSystemView().getHomeDirectory().toString()+"/janr/";
	    public static String janrPath = "run/";
	    
	    public static String infile   = janrPath+"janr_input_e1b";
	    public static String resfile  = janrPath+"res-q0.4.inp";
	    public static String parifile = janrPath+"par_e1b.inp";
	    public static String resifile = janrPath+"res-q0.4.inp";
	
        public static final int MAX_PHI_INDEX = 6;
        public static final int MAX_IS_INDEX = 3;
        public static final int MAX_J_INDEX = 2;

        
 // MULTIPOL        
        public static final String[] fnamr = {
                "rm0p.dat","rm0m.dat","re0p.dat","re0m.dat","rs0p.dat","rs0m.dat",
                "rm1p.dat","rm1m.dat","re1p.dat","re1m.dat","rs1p.dat","rs1m.dat",
                "rm2p.dat","rm2m.dat","re2p.dat","re2m.dat","rs2p.dat","rs2m.dat",
                "rm3p.dat","rm3m.dat","re3p.dat","re3m.dat","rs3p.dat","rs3m.dat"
            };
        
        public static final String[] fnamhigh = {
                "rm4p.dat","rm4m.dat","re4p.dat","re4m.dat","rs4p.dat","rs4m.dat",
                "rm5p.dat","rm5m.dat","re5p.dat","re5m.dat","rs5p.dat","rs5m.dat"
            };
        
        public static final String[] fnamrr = {
                "rrm0p.dat","rrm0m.dat","rre0p.dat","rre0m.dat",
                "rrs0p.dat","rrs0m.dat",
                "rrm1p.dat","rrm1m.dat","rre1p.dat","rre1m.dat","rrs1p.dat","rrs1m.dat",
                "rrm2p.dat","rrm2m.dat","rre2p.dat","rre2m.dat","rrs2p.dat","rrs2m.dat",
                "rrm3p.dat","rrm3m.dat","rre3p.dat","rre3m.dat","rrs3p.dat","rrs3m.dat"
            };
	
// AMPLITUDES
	    public static final int MAX1res = 13;
	    public static final int MAX3res = 13;
	    public static final int MAX13res = 40;

	    public static double am3[] = new double[MAX3res];
	    public static double ae3[] = new double[MAX3res];
	    public static double as3[] = new double[MAX3res];
	    public static double am1[][] = new double[2][MAX1res];
	    public static double ae1[][] = new double[2][MAX1res];
	    public static double as1[][] = new double[2][MAX1res];
	    public static double aa3[] = new double[50];
	    public static double aa1[] = new double[50];
	    public static double sa1[] = new double[50];
	    public static double BreitE3[][] = new double[2][MAX3res];
	    public static double BreitE1[][][] = new double[2][2][MAX1res];
	    public static double BreitM3[][] = new double[2][MAX3res];
	    public static double BreitS3[][] = new double[2][MAX3res];
	    public static double BreitM1[][][] = new double[2][2][MAX1res];
	    public static double BreitS1[][][] = new double[2][2][MAX1res];
	    public static double resM;

	
// BACK_AMP
	    public static double[]  BM3 = new double[8];
	    public static double[]  BE3 = new double[8];
	    public static double[]  BS3 = new double[8];
	    public static double[] BbM3 = new double[8];
	    public static double[] BbE3 = new double[8];
	    public static double[] BbS3 = new double[8];
	    public static double[]  RM3 = new double[8];
	    public static double[]  RE3 = new double[8];
	    public static double[]  RS3 = new double[8];

	    public static double[][] BM1 = new double[2][9];
	    public static double[][] BE1 = new double[2][9];
	    public static double[][] BS1 = new double[2][9];
	    public static double[][] BbM1 = new double[2][9];
	    public static double[][] BbE1 = new double[2][9];
	    public static double[][] BbS1 = new double[2][9];
	    public static double[][] RM1 = new double[2][9];
	    public static double[][] RE1 = new double[2][9];
	    public static double[][] RS1 = new double[2][9];

	    public static double[] E0 = new double[2];
	    public static double[] M1 = new double[2];
	    public static double[] S0 = new double[2];
	    public static double[] S1 = new double[2];

	
// BORN_TERMS
		public static double [][][] br = new double[2][6][3];
		
// CROSS_SEC 
		/*
		public static double[] sigmaT  = new double[3];
		public static double[] sigmaL  = new double[3];
		public static double[] sigmaTT = new double[3];
		public static double[] sigmaTL = new double[3];
		public static double[] sigma   = new double[3];
		public static double[] sigmaP  = new double[3];
		public static double[] sigmaI  = new double[3];
*/
		
// DELTA_DAT 
		public static double[] xM  = new double[74];
		public static double[] xE  = new double[74];
		public static double[] xS  = new double[74];
		public static double[] xOM = new double[74];
		public static double[] xOE = new double[74];
		public static double[] xOS = new double[74];
		public static double[] xH  = new double[74];
		public static double[] yM  = new double[74];
		public static double[] yE  = new double[74];
		public static double[] yS  = new double[74];
		public static double[] yOM = new double[74];
		public static double[] yOE = new double[74];
		public static double[] yOS = new double[74];
		public static double[] yH  = new double[74];
		
// EBAC
		

// ETTAPAR 
	   public static double[][] ett = {
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0}};
		
// FITPAR
	    public static final int MAXpar = 99;
	    public static final int f_out = 49;
	    
	    public static String[]       pname = new String[MAXpar];
	    public static double[] start_value = new double[MAXpar];
	    public static double[]   step_size = new double[MAXpar];
	    public static double[]     low_lim = new double[MAXpar];
	    public static double[]      up_lim = new double[MAXpar];
	    public static int[]       par_stat = new int[MAXpar];
	    public static int[]          plist = new int[MAXpar];
	    public static String[]    res_name = new String[MAXpar];
	    
	    public static int pi_type;
	    public static int pi_type2fit;
	    public static int obsv2fit;
	    public static int ndf;
	    public static int parms;
	    public static int lmini;
	    public static int lmino;
	    public static int lminf;
	    
	    public static double cm1, ce1, cs1;
	    public static double cm2, ce2, cs2;
	    public static double cm3, ce3, cs3;
	    public static double cmp1, csp1, cep2, csp2, cep5, csp5;
	    public static double cmp3, cep3, csp3;
	    public static double cmp4, cep4, csp4;
	    public static double cmp6, cep6, csp6;
	    public static double cmp7, cep7, csp7;
	    public static double cmp8, cep8, csp8;
	    public static double cmp9, cep9;
	    public static double cmp10,cep10;
	    public static double cmp11,cep11;
	    public static double cspa, a_pip;
	    
	    public static double ChiMax;
	    public static double fix_param;
	    public static double tolorance;
	    public static double maxcalls;
	    public static double errorup;
	
// GPAR 
		public static double w2,E,Egamma,Epion,m2pion,q2pion,meta,qeta;
		public static double pi,a2,mn,mp,mn22,mp22,mn2,u,t,qk,F2p,F2n;
	
// HELICITY 
		public static double [][]  H1 = new double[3][2];
		public static double [][]  H2 = new double[3][2];
		public static double [][]  H3 = new double[3][2];
		public static double [][]  H4 = new double[3][2];
		public static double [][]  H5 = new double[3][2];
		public static double [][]  H6 = new double[3][2];
		public static double [][] sig = new double[3][2];
	
//  HIGHMUL 
		public static double [][] phih = new double[6][3];
		public static double[]     E4m = new double[3];
		public static double[]     E4p = new double[3];
		public static double[]     M4m = new double[3];
		public static double[]     M4p = new double[3];
		public static double[]     S4m = new double[3];
		public static double[]     S4p = new double[3];
		public static double[]     E5m = new double[3];
		public static double[]     E5p = new double[3];
		public static double[]     M5m = new double[3];
		public static double[]     M5p = new double[3];
		public static double[]     S5m = new double[3];
		public static double[]     S5p = new double[3];
	
// IMA 
		public static double [][] ba = new double[6][3];
			
// IFINTG
//		public static double lm,is,i,j,w,q2;
		
// INPAR
		public static String 	fitopt;
		public static String 	grfopt;
		public static double 	n1;
		public static double 	n2;
		public static int   	nebac;
		
// INT_FUNC
		public static double[][] h = new double[6][3];
		public static int jk,km;
		
// JANR_DATA
	    public static final int MAXpoints = 100000;
	    
	    public static int[][] ibad = new int[2][100];
	    public static int nbad;
	    public static int datform;	    
	    public static List<String> fname = new ArrayList<>();
	    public static String parofile;
	    public static int[] point_stat = new int[MAXpoints];
// HJANR_HCRS
	    public static int MAX_K=8200, MAX_I=3, MAX_OBS=4, MAX_J=7;
	    public static int nplt;
	    public static int[][]   ndims = new    int[MAX_OBS][MAX_I];  
	    public static double[][][] xx = new double[MAX_OBS][MAX_J][MAX_K]; 
	    public static double[][]  crs = new double[MAX_OBS][MAX_K];
	    public static double[][] crse = new double[MAX_OBS][MAX_K];
	    public static String[]    lab = new String[4];
	    
// HJANR_HLIM
	    private static int[] c = new int[4];
	    private static int[] w = new int[4];
	    private static int[] f = new int[4];
	    
// HJANR_HPLT	    
	    public static double[]    xplt = new double[100];
	    public static double[][] xpltt = new double[100][30];
	    public static int nt1;
	    
//HJANR_VAL	    
	    public static double[] xnew = new double[36];
	    
// JANR_TAB
	    public static int maxevents;
	    public static int maxiter;
	    public static int randomseed;

	    public static final int MAX_Q2_GRID = 20;
	    public static final int MAX_W_GRID = 100;
	    public static final int MAX_COS_GRID = 20;
	    public static final int MAX_PHI_GRID = 30;

	    public static double wmin;
	    public static double wmax;
	    public static double wstep;
	    public static double q2min=0;
	    public static double q2max=1;
	    public static double q2step;
	    public static double cosmin;
	    public static double cosmax;
	    public static double cosstep;
	    public static double phimin;
	    public static double phimax;
	    public static double phistep;
	    public static double ebeam;
	    public static double IntAccur;

	    public static double[][][][] csTab = new double[MAX_Q2_GRID][MAX_W_GRID][MAX_COS_GRID][MAX_PHI_GRID];

	    public static String tabfile;
	    public static String inffile;

// MULTAMPL
	    public static final int Maxmpoints = 101;

	    public static double[] wtab = new double[Maxmpoints];
	    
	    public static double[][][] aem = new double[Maxmpoints][6][3];
	    public static double[][][] aep = new double[Maxmpoints][6][3];
	    
	    public static double[][][] asm = new double[Maxmpoints][6][3];
	    public static double[][][] asp = new double[Maxmpoints][6][3];
	    
	    public static double[][][] amm = new double[Maxmpoints][6][3];
	    public static double[][][] amp = new double[Maxmpoints][6][3];
	    
	    public static double[][][] raem = new double[Maxmpoints][4][3];
	    public static double[][][] raep = new double[Maxmpoints][4][3];
	    
	    public static double[][][] rasm = new double[Maxmpoints][4][3];
	    public static double[][][] rasp = new double[Maxmpoints][4][3];
	    
	    public static double[][][] ramm = new double[Maxmpoints][4][3];
	    public static double[][][] ramp = new double[Maxmpoints][4][3];
	    
// MULTIPOLES	    
	    public static final int ROWS = 2;
	    public static final int COLS = 2;

	    public static double[][] ReacM2 = new double[ROWS][COLS];
	    public static double[][] ReacM3 = new double[ROWS][COLS];
	    public static double[][] ReacM4 = new double[ROWS][COLS];
	    public static double[][] ReacM5 = new double[ROWS][COLS];
	    public static double[][] ReacM6 = new double[ROWS][COLS];
	    public static double[][] ReacM7 = new double[ROWS][COLS];

	    public static double[][] ReacE1 = new double[ROWS][COLS];
	    public static double[][] ReacE3 = new double[ROWS][COLS];
	    public static double[][] ReacE4 = new double[ROWS][COLS];
	    public static double[][] ReacE5 = new double[ROWS][COLS];
	    public static double[][] ReacE6 = new double[ROWS][COLS];
	    public static double[][] ReacE7 = new double[ROWS][COLS];

	    public static double[][] ReacS1 = new double[ROWS][COLS];
	    public static double[][] ReacS2 = new double[ROWS][COLS];
	    public static double[][] ReacS3 = new double[ROWS][COLS];
	    public static double[][] ReacS4 = new double[ROWS][COLS];
	    public static double[][] ReacS5 = new double[ROWS][COLS];
	    public static double[][] ReacS6 = new double[ROWS][COLS];
	    public static double[][] ReacS7 = new double[ROWS][COLS];
	    
	    public static double[][] MM3 = new double[2][8]; 
	    public static double[][] ME3 = new double[2][8]; 
	    public static double[][] MS3 = new double[2][8]; 
	    public static double[][][] MM1 = new double[2][2][9]; 
	    public static double[][][] ME1 = new double[2][2][9]; 
	    public static double[][][] MS1 = new double[2][2][9]; 

	    
	    public static final int N51 = 51;
	    public static final int N13 = 13;
	    	    
	    public static double[] WSa = new double[N51];
	    
	    public static double[] mp33r = new double[N51];
	    public static double[] mp33i = new double[N51];
	    
	    public static double[] s11r = new double[N51];
	    public static double[] s11i = new double[N51];
	    
	    public static double[] s31r = new double[N51];
	    public static double[] s31i = new double[N51];
	    
	    public static double[] p11r = new double[N51];
	    public static double[] p11i = new double[N51];
	    
	    public static double[] p13r = new double[N51];
	    public static double[] p13i = new double[N51];
	    
	    public static double[] p31r = new double[N51];
	    public static double[] p31i = new double[N51];
	    
	    public static double[] p33r = new double[N51];
	    public static double[] p33i = new double[N51];
	    
	    public static double[] d13r = new double[N51];
	    public static double[] d13i = new double[N51];
	    
	    public static double[] d15r = new double[N51];
	    public static double[] d15i = new double[N51];
	    
	    public static double[] d33r = new double[N51];
	    public static double[] d33i = new double[N51];
	    
	    public static double[] f15r = new double[N51];
	    public static double[] f15i = new double[N51];
	    
	    public static double[] f37r = new double[N51];
	    public static double[] f37i = new double[N51];
	    
	    public static double[] simr = new double[N13];
	    public static double[] simi = new double[N13];
	    
	    public static double[] pimr = new double[N13];
	    public static double[] pimi = new double[N13];
		
// PHIAMPL
	    public static double[][][] phi_amp = new double[6][2][3];
	    
// REMULT
	    public static double[][][][] Rint = new double[6][3][6][2];  
	    public static double[][][]    Mlp = new double[6][3][2];
	    public static double[][][]    Mlm = new double[6][3][2];
	    public static double[][][]    Elp = new double[6][3][2];
	    public static double[][][]    Elm = new double[6][3][2];
	    public static double[][][]    Slp = new double[6][3][2];
	    public static double[][][]    Slm = new double[6][3][2];
	    
// REPHI 
	    public static double[][][] PHisot = new double[2][6][3];  
	    
// RESFUNC
	    public static double[] sigmaT = new double[2];
	    public static double[] sigmaL = new double[2];
	    public static double[] sigmaTT = new double[2];
	    public static double[] sigmaTL = new double[2];
	    public static double[] sigma = new double[2];
	    public static double[] sigmaTLP = new double[2];
	    public static double[] sigmaTy = new double[2];
	    public static double[] sigmaLy = new double[2];
	    public static double[] sigmaTTx = new double[2];
	    public static double[] sigmaTTy = new double[2];
	    public static double[] sigmaTTz = new double[2];
	    public static double[] sigmaTLx = new double[2];
	    public static double[] sigmaTLy = new double[2];
	    public static double[] sigmaTLz = new double[2];
	    public static double[] sigmaTTPx = new double[2];
	    public static double[] sigmaTTPz = new double[2];
	    public static double[] sigmaTLPx = new double[2];
	    public static double[] sigmaTLPy = new double[2];
	    public static double[] sigmaTLPz = new double[2];
	    public static double[][] robs = new double[10][2];
	    
// RESONANCES
	    public static final int Nres = 26;
	    public static final int MAXres = 34;
	    
	    public static double[] mres = new double[MAXres];
	    public static double[] gres = new double[MAXres];
	    public static double[] xres = new double[MAXres];
	    public static double[]  eta = new double[MAXres];
	    
	    public static int[] langul = new int[MAXres];
	    public static int[] lprime = new int[MAXres];
	    public static int[]   jres = new int[MAXres];
	    public static int[]   ires = new int[MAXres];
	    public static int[]   lres = new int[MAXres];
	    	
}
