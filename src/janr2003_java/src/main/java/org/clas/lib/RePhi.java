package org.clas.lib;

public class RePhi extends Constants {
	
	double on1,on2,op1,op2,om1,om2,wm,wp,xx;
	double br1,br2,br3,br4,br5,br6;
    double[] AK = new double[6];
    double W, qq2;
    	
	public RePhi() {
		
	}

    public void rePhi(double S, double Q2) {
        W = Math.sqrt(S);
        qq2 = Q2;
        if (Q2 == 0) qq2 = 0.00001f;
        
        on1 = W - Egamma;
        on2 = W - Epion;
        op1 = on1 + mn;
        om1 = on1 - mn;
        op2 = on2 + mn;
        om2 = on2 - mn;
        
        double Wm = W - mn;
        double Wp = W + mn;
        double XX = Egamma * qk - Epion * qq2;
        
        AK[0] = (Math.sqrt(op1 * op2) / w2);
        AK[1] = (float) (Math.sqrt(om1 * om2) / w2);
        AK[2] = AK[1] * op2;
        AK[3] = AK[0] * om2;
        AK[4] = (float) (Math.sqrt(om1 * op2) / (qq2 * w2));
        AK[5] = (float) (Math.sqrt(op1 * om2) / (qq2 * w2));
        
        // Loop by Isospin states (+), (-) & (0)
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 2; j++) {            	
                br5 = br[j][5][i] + (E + qk) * br[j][6][i] / 2.0f;
                
                PHisot[j][0][i] =  AK[0] * (Wm * br[j][0][i] - br5);
                PHisot[j][1][i] = -AK[1] * (Wp * br[j][0][i] + br5);
                
                br1 = 2.0f * br[j][2][i] - br[j][1][i];
                br2 = br[j][6][i] / 2.0f - br[j][3][i];
                br3 = qq2 * br[j][0][i];
                
                PHisot[j][2][i] = AK[2] * ( br1 + Wp * br2);
                PHisot[j][3][i] = AK[3] * (-br1 + Wm * br2);
                
                br4 = op1 * (br3 + Wm * br5 + w2 * om1 * (br[j][1][i] - Wp * br[j][6][i] / 2.0f));
                PHisot[j][4][i] = AK[4] * (br4 - XX * (br1 + Wp * br2));
                
                br6 = -om1 * (br3 - Wp * br5 + w2 * op1 * (br[j][1][i] + Wm * br[j][6][i] / 2.0f));
                PHisot[j][5][i] = AK[5] * (br6 + XX * (br1 - Wm * br2));
            }
        }
    }
    
}
