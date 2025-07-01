package org.clas.lib;

import java.util.Arrays;

public class Constants {
	
        public static final int MAX_PHI_INDEX = 6;
        public static final int MAX_IS_INDEX = 3;
        public static final int MAX_J_INDEX = 2;
	
// AMPLITUDES
	    public static final int MAX1res = 13;
	    public static final int MAX3res = 13;
	    public static final int MAX13res = 40;

	    public double am3[] = new double[MAX3res];
	    public double ae3[] = new double[MAX3res];
	    public double as3[] = new double[MAX3res];
	    public double am1[][] = new double[2][MAX1res];
	    public double ae1[][] = new double[2][MAX1res];
	    public double as1[][] = new double[2][MAX1res];
	    public double aa3[] = new double[50];
	    public double aa1[] = new double[50];
	    public double sa1[] = new double[50];
	    public double BreitE3[][] = new double[2][MAX3res];
	    public double BreitE1[][][] = new double[2][2][MAX1res];
	    public double BreitM3[][] = new double[2][MAX3res];
	    public double BreitS3[][] = new double[2][MAX3res];
	    public double BreitM1[][][] = new double[2][2][MAX1res];
	    public double BreitS1[][][] = new double[2][2][MAX1res];
	    public double resM;

	
// BACK_AMP
	    public float[] BM3 = new float[8];
	    public float[] BE3 = new float[8];
	    public float[] BS3 = new float[8];
	    public float[] BbM3 = new float[8];
	    public float[] BbE3 = new float[8];
	    public float[] BbS3 = new float[8];
	    public float[] RM3 = new float[8];
	    public float[] RE3 = new float[8];
	    public float[] RS3 = new float[8];

	    public float[][] BM1 = new float[2][9];
	    public float[][] BE1 = new float[2][9];
	    public float[][] BS1 = new float[2][9];
	    public float[][] BbM1 = new float[2][9];
	    public float[][] BbE1 = new float[2][9];
	    public float[][] BbS1 = new float[2][9];
	    public float[][] RM1 = new float[2][9];
	    public float[][] RE1 = new float[2][9];
	    public float[][] RS1 = new float[2][9];

	    public float[] E0 = new float[2];
	    public float[] M1 = new float[2];
	    public float[] S0 = new float[2];
	    public float[] S1 = new float[2];

	
// BORN_TERMS
		public float [][][] br = new float[2][6][3];
		
// CROSS_SEC 
		/*
		public float[] sigmaT  = new float[3];
		public float[] sigmaL  = new float[3];
		public float[] sigmaTT = new float[3];
		public float[] sigmaTL = new float[3];
		public float[] sigma   = new float[3];
		public float[] sigmaP  = new float[3];
		public float[] sigmaI  = new float[3];
*/
		
// DELTA_DAT 
		public float[] xM  = new float[74];
		public float[] xE  = new float[74];
		public float[] xS  = new float[74];
		public float[] xOM = new float[74];
		public float[] xOE = new float[74];
		public float[] xOS = new float[74];
		public float[] xH  = new float[74];
		public float[] yM  = new float[74];
		public float[] yE  = new float[74];
		public float[] yS  = new float[74];
		public float[] yOM = new float[74];
		public float[] yOE = new float[74];
		public float[] yOS = new float[74];
		public float[] yH  = new float[74];

// ETTAPAR 
	   public double[][] ett = {
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0},
	            {1.0, 1.0, 1.0}};
		
// FITPAR
	    public static final int MAXpar = 99;
	    public static final int f_out = 49;
	    
	    public static String[]      pname = new String[MAXpar];
	    public static float[] start_value = new float[MAXpar];
	    public static float[]   step_size = new float[MAXpar];
	    public static float[]     low_lim = new float[MAXpar];
	    public static float[]      up_lim = new float[MAXpar];
	    public static int[]      par_stat = new int[MAXpar];
	    public static int[]         plist = new int[MAXpar];
	    
	    public static int pi_type;
	    public static int pi_type2fit;
	    public static int obsv2fit;
	    public static int ndf;
	    public static int parms;
	    public static int lmini;
	    public static int lmino;
	    public static int lminf;
	    
	    double cm1, ce1, cs1;
	    double cm2, ce2, cs2;
	    double cm3, ce3, cs3;
	    double cmp1, csp1, cep2, csp2, cep5, csp5;
	    double cmp3, cep3, csp3;
	    double cmp4, cep4, csp4;
	    double cmp6, cep6, csp6;
	    double cmp7, cep7, csp7;
	    double cmp8, cep8, csp8;
	    double cmp9, cep9;
	    double cmp10,cep10;
	    double cmp11,cep11;
	    double cspa, a_pip;
	    
	    public static float ChiMax;
	    public static float fix_param;
	    public static float tolorance;
	    public static float maxcalls;
	    public static float errorup;
	
// GPAR 
		public double w2,E,Egamma,Epion,m2pion,q2pion,meta,qeta;
		public double pi,a2,mn,mp,mn22,mp22,mn2,u,t,qk,F2p,F2n;
	
// HELICITY 
		public float [][]  h1 = new float[3][2];
		public float [][]  h2 = new float[3][2];
		public float [][]  h3 = new float[3][2];
		public float [][]  h4 = new float[3][2];
		public float [][]  h5 = new float[3][2];
		public float [][]  h6 = new float[3][2];
		public float [][] sig = new float[3][2];
	
//  HIGHMUL 
		public float [][] phih = new float[6][3];
		public float[]     E4m = new float[3];
		public float[]     E4p = new float[3];
		public float[]     M4m = new float[3];
		public float[]     M4p = new float[3];
		public float[]     S4m = new float[3];
		public float[]     S4p = new float[3];
		public float[]     E5m = new float[3];
		public float[]     E5p = new float[3];
		public float[]     M5m = new float[3];
		public float[]     M5p = new float[3];
		public float[]     S5m = new float[3];
		public float[]     S5p = new float[3];
	
// IMA 
		public float [][] ba = new float[6][3];
			
// INPAR
		public String 	fitopt;
		public String 	grfopt;
		public float 	n1;
		public float 	n2;
		public int   	nebac;
		
// INT_FUNC
		public float[][] h = new float[6][3];
		public int jk,km;
		
// JANR_DATA
	    public static final int MAXpoints = 100000;
	    
	    public static int[][] ibad = new int[2][100];
	    public static int nbad;
	    public static int datform;
	    
	    public static String[] fname = new String[4];
	    public static String parofile;
	    public static int[] point_stat = new int[MAXpoints];

// JANR_TAB
	    public int maxevents;
	    public int maxiter;
	    public int randomseed;

	    public static final int MAX_Q2_GRID = 20;
	    public static final int MAX_W_GRID = 100;
	    public static final int MAX_COS_GRID = 20;
	    public static final int MAX_PHI_GRID = 30;

	    public float wmin;
	    public float wmax;
	    public float wstep;
	    public float q2min;
	    public float q2max;
	    public float q2step;
	    public float cosmin;
	    public float cosmax;
	    public float cosstep;
	    public float phimin;
	    public float phimax;
	    public float phistep;
	    public float ebeam;
	    public float IntAccur;

	    public double[][][][] csTab = new double[MAX_Q2_GRID][MAX_W_GRID][MAX_COS_GRID][MAX_PHI_GRID];

	    public String tabfile;
	    public String inffile;

// MULTAMPL
	    public static final int Maxmpoints = 101;

	    // Declare arrays with the specified dimensions
	    public double[] Wtab = new double[Maxmpoints];
	    
	    public double[][][] AEm = new double[Maxmpoints][6][3];
	    public double[][][] AEp = new double[Maxmpoints][6][3];
	    
	    public double[][][] ASm = new double[Maxmpoints][6][3];
	    public double[][][] ASp = new double[Maxmpoints][6][3];
	    
	    public double[][][] AMm = new double[Maxmpoints][6][3];
	    public double[][][] AMp = new double[Maxmpoints][6][3];
	    
	    public double[][][] RAEm = new double[Maxmpoints][4][3];
	    public double[][][] RAEp = new double[Maxmpoints][4][3];
	    
	    public double[][][] RASm = new double[Maxmpoints][4][3];
	    public double[][][] RASp = new double[Maxmpoints][4][3];
	    
	    public double[][][] RAMm = new double[Maxmpoints][4][3];
	    public double[][][] RAMp = new double[Maxmpoints][4][3];
	    
// MULTIPOLES
	    
	    public static final int ROWS = 2;
	    public static final int COLS = 2;

	    public double[][] ReacM2 = new double[ROWS][COLS];
	    public double[][] ReacM3 = new double[ROWS][COLS];
	    public double[][] ReacM4 = new double[ROWS][COLS];
	    public double[][] ReacM5 = new double[ROWS][COLS];
	    public double[][] ReacM6 = new double[ROWS][COLS];
	    public double[][] ReacM7 = new double[ROWS][COLS];

	    public double[][] ReacE1 = new double[ROWS][COLS];
	    public double[][] ReacE3 = new double[ROWS][COLS];
	    public double[][] ReacE4 = new double[ROWS][COLS];
	    public double[][] ReacE5 = new double[ROWS][COLS];
	    public double[][] ReacE6 = new double[ROWS][COLS];
	    public double[][] ReacE7 = new double[ROWS][COLS];

	    public double[][] ReacS1 = new double[ROWS][COLS];
	    public double[][] ReacS2 = new double[ROWS][COLS];
	    public double[][] ReacS3 = new double[ROWS][COLS];
	    public double[][] ReacS4 = new double[ROWS][COLS];
	    public double[][] ReacS5 = new double[ROWS][COLS];
	    public double[][] ReacS6 = new double[ROWS][COLS];
	    public double[][] ReacS7 = new double[ROWS][COLS];
	    
	    public static final int N51 = 51;
	    public static final int N13 = 13;
	    	    
	    public double[] wsa = new double[N51];
	    
	    public double[] mp33r = new double[N51];
	    public double[] mp33i = new double[N51];
	    
	    public double[] s11r = new double[N51];
	    public double[] s11i = new double[N51];
	    
	    public double[] s31r = new double[N51];
	    public double[] s31i = new double[N51];
	    
	    public double[] p11r = new double[N51];
	    public double[] p11i = new double[N51];
	    
	    public double[] p13r = new double[N51];
	    public double[] p13i = new double[N51];
	    
	    public double[] p31r = new double[N51];
	    public double[] p31i = new double[N51];
	    
	    public double[] p33r = new double[N51];
	    public double[] p33i = new double[N51];
	    
	    public double[] d13r = new double[N51];
	    public double[] d13i = new double[N51];
	    
	    public double[] d15r = new double[N51];
	    public double[] d15i = new double[N51];
	    
	    public double[] d33r = new double[N51];
	    public double[] d33i = new double[N51];
	    
	    public double[] f15r = new double[N51];
	    public double[] f15i = new double[N51];
	    
	    public double[] f37r = new double[N51];
	    public double[] f37i = new double[N51];
	    
	    public double[] simr = new double[N13];
	    public double[] simi = new double[N13];
	    
	    public double[] pimr = new double[N13];
	    public double[] pimi = new double[N13];
		
// PHIAMPL
	    public double[][][] phi_amp = new double[6][2][3];
	    
// REMULT
	    public float[][][][] Rint = new float[6][3][6][2];
	    public float[][][] Mlp = new float[6][3][2];
	    public float[][][] Mlm = new float[6][3][2];
	    public float[][][] Elp = new float[6][3][2];
	    public float[][][] Elm = new float[6][3][2];
	    public float[][][] Slp = new float[6][3][2];
	    public float[][][] Slm = new float[6][3][2];
	    
// REPHI 
	    public double[][][] PHisot = new double[2][6][3];
	    
// RESFUNC
	    public double[] sigmaT = new double[2];
	    public double[] sigmaL = new double[2];
	    public double[] sigmaTT = new double[2];
	    public double[] sigmaTL = new double[2];
	    public double[] sigma = new double[2];
	    public double[] sigmaTLP = new double[2];
	    public double[] sigmaTy = new double[2];
	    public double[] sigmaLy = new double[2];
	    public double[] sigmaTTx = new double[2];
	    public double[] sigmaTTy = new double[2];
	    public double[] sigmaTTz = new double[2];
	    public double[] sigmaTLx = new double[2];
	    public double[] sigmaTLy = new double[2];
	    public double[] sigmaTLz = new double[2];
	    public double[] sigmaTTPx = new double[2];
	    public double[] sigmaTTPz = new double[2];
	    public double[] sigmaTLPx = new double[2];
	    public double[] sigmaTLPy = new double[2];
	    public double[] sigmaTLPz = new double[2];
	    public double[][] robs = new double[10][2];
	    
// RESONANCES
	    public static final int MAXres = 40;
	    
	    public double[] mres = new double[MAXres];
	    public double[] gres = new double[MAXres];
	    public double[] xres = new double[MAXres];
	    public double[] eta = new double[MAXres];
	    public double[] lres = new double[MAXres];
	    
	    public int[] langul = new int[MAXres];
	    public int[] lprime = new int[MAXres];
	    public int[] jres = new int[MAXres];
	    public int[] ires = new int[MAXres];
	    	
}
