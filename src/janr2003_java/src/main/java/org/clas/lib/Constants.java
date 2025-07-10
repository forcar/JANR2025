package org.clas.lib;

import java.util.Arrays;

import javax.swing.filechooser.FileSystemView;

public class Constants {
	
        //String janrPath = FileSystemView.getFileSystemView().getHomeDirectory().toString()+"/janr/";
	    public static String janrPath = "run/";
	
        public static final int MAX_PHI_INDEX = 6;
        public static final int MAX_IS_INDEX = 3;
        public static final int MAX_J_INDEX = 2;
	
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
	    public static float[] BM3 = new float[8];
	    public static float[] BE3 = new float[8];
	    public static float[] BS3 = new float[8];
	    public static float[] BbM3 = new float[8];
	    public static float[] BbE3 = new float[8];
	    public static float[] BbS3 = new float[8];
	    public static float[] RM3 = new float[8];
	    public static float[] RE3 = new float[8];
	    public static float[] RS3 = new float[8];

	    public static float[][] BM1 = new float[2][9];
	    public static float[][] BE1 = new float[2][9];
	    public static float[][] BS1 = new float[2][9];
	    public static float[][] BbM1 = new float[2][9];
	    public static float[][] BbE1 = new float[2][9];
	    public static float[][] BbS1 = new float[2][9];
	    public static float[][] RM1 = new float[2][9];
	    public static float[][] RE1 = new float[2][9];
	    public static float[][] RS1 = new float[2][9];

	    public static float[] E0 = new float[2];
	    public static float[] M1 = new float[2];
	    public static float[] S0 = new float[2];
	    public static float[] S1 = new float[2];

	
// BORN_TERMS
		public static float [][][] br = new float[2][6][3];
		
// CROSS_SEC 
		/*
		public static float[] sigmaT  = new float[3];
		public static float[] sigmaL  = new float[3];
		public static float[] sigmaTT = new float[3];
		public static float[] sigmaTL = new float[3];
		public static float[] sigma   = new float[3];
		public static float[] sigmaP  = new float[3];
		public static float[] sigmaI  = new float[3];
*/
		
// DELTA_DAT 
		public static float[] xM  = new float[74];
		public static float[] xE  = new float[74];
		public static float[] xS  = new float[74];
		public static float[] xOM = new float[74];
		public static float[] xOE = new float[74];
		public static float[] xOS = new float[74];
		public static float[] xH  = new float[74];
		public static float[] yM  = new float[74];
		public static float[] yE  = new float[74];
		public static float[] yS  = new float[74];
		public static float[] yOM = new float[74];
		public static float[] yOE = new float[74];
		public static float[] yOS = new float[74];
		public static float[] yH  = new float[74];

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
	    
	    public static float ChiMax;
	    public static float fix_param;
	    public static float tolorance;
	    public static float maxcalls;
	    public static float errorup;
	
// GPAR 
		public static double w2,E,Egamma,Epion,m2pion,q2pion,meta,qeta;
		public static double pi,a2,mn,mp,mn22,mp22,mn2,u,t,qk,F2p,F2n;
	
// HELICITY 
		public static float [][]  h1 = new float[3][2];
		public static float [][]  h2 = new float[3][2];
		public static float [][]  h3 = new float[3][2];
		public static float [][]  h4 = new float[3][2];
		public static float [][]  h5 = new float[3][2];
		public static float [][]  h6 = new float[3][2];
		public static float [][] sig = new float[3][2];
	
//  HIGHMUL 
		public static float [][] phih = new float[6][3];
		public static float[]     E4m = new float[3];
		public static float[]     E4p = new float[3];
		public static float[]     M4m = new float[3];
		public static float[]     M4p = new float[3];
		public static float[]     S4m = new float[3];
		public static float[]     S4p = new float[3];
		public static float[]     E5m = new float[3];
		public static float[]     E5p = new float[3];
		public static float[]     M5m = new float[3];
		public static float[]     M5p = new float[3];
		public static float[]     S5m = new float[3];
		public static float[]     S5p = new float[3];
	
// IMA 
		public static float [][] ba = new float[6][3];
			
// INPAR
		public static String 	fitopt;
		public static String 	grfopt;
		public static float 	n1;
		public static float 	n2;
		public static int   	nebac;
		
// INT_FUNC
		public static float[][] h = new float[6][3];
		public static int jk,km;
		
// JANR_DATA
	    public static final int MAXpoints = 100000;
	    
	    public static int[][] ibad = new int[2][100];
	    public static int nbad;
	    public static int datform;	    
	    public static String[] fname = new String[4];
	    public static String parofile;
	    public static int[] point_stat = new int[MAXpoints];

// JANR_TAB
	    public static int maxevents;
	    public static int maxiter;
	    public static int randomseed;

	    public static final int MAX_Q2_GRID = 20;
	    public static final int MAX_W_GRID = 100;
	    public static final int MAX_COS_GRID = 20;
	    public static final int MAX_PHI_GRID = 30;

	    public static float wmin;
	    public static float wmax;
	    public static float wstep;
	    public static float q2min;
	    public static float q2max;
	    public static float q2step;
	    public static float cosmin;
	    public static float cosmax;
	    public static float cosstep;
	    public static float phimin;
	    public static float phimax;
	    public static float phistep;
	    public static float ebeam;
	    public static float IntAccur;

	    public static double[][][][] csTab = new double[MAX_Q2_GRID][MAX_W_GRID][MAX_COS_GRID][MAX_PHI_GRID];

	    public static String tabfile;
	    public static String inffile;

// MULTAMPL
	    public static final int Maxmpoints = 101;

	    // Declare arrays with the specified dimensions
	    public static double[] Wtab = new double[Maxmpoints];
	    
	    public static double[][][] AEm = new double[Maxmpoints][6][3];
	    public static double[][][] AEp = new double[Maxmpoints][6][3];
	    
	    public static double[][][] ASm = new double[Maxmpoints][6][3];
	    public static double[][][] ASp = new double[Maxmpoints][6][3];
	    
	    public static double[][][] AMm = new double[Maxmpoints][6][3];
	    public static double[][][] AMp = new double[Maxmpoints][6][3];
	    
	    public static double[][][] RAEm = new double[Maxmpoints][4][3];
	    public static double[][][] RAEp = new double[Maxmpoints][4][3];
	    
	    public static double[][][] RASm = new double[Maxmpoints][4][3];
	    public static double[][][] RASp = new double[Maxmpoints][4][3];
	    
	    public static double[][][] RAMm = new double[Maxmpoints][4][3];
	    public static double[][][] RAMp = new double[Maxmpoints][4][3];
	    
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
	    
	    public static final int N51 = 51;
	    public static final int N13 = 13;
	    	    
	    public static double[] wsa = new double[N51];
	    
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
	    public static float[][][][] Rint = new float[6][3][6][2];
	    public static float[][][] Mlp = new float[6][3][2];
	    public static float[][][] Mlm = new float[6][3][2];
	    public static float[][][] Elp = new float[6][3][2];
	    public static float[][][] Elm = new float[6][3][2];
	    public static float[][][] Slp = new float[6][3][2];
	    public static float[][][] Slm = new float[6][3][2];
	    
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
	    public static double[] eta = new double[MAXres];
	    public static double[] lres = new double[MAXres];
	    
	    public static int[] langul = new int[MAXres];
	    public static int[] lprime = new int[MAXres];
	    public static int[] jres = new int[MAXres];
	    public static int[] ires = new int[MAXres];
	    	
}
