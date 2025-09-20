package org.clas.lib;

// https://github.com/apc-llc/cernlib/blob/ef83476ab09eed4972aff95685fcb27c24e6a18c/2006/src/mathlib/gen/c/gammf64.F#L11 (FORTRAN)
// https://codingfleet.com/code-converter/java/ (DeepSeek V3) (JAVA)

public class GAMMF {
	
    private static final double PI = 3.14159265358979324;
    private static final double[] C = {
        3.65738772508338244,
        1.95754345666126827,
        0.33829711382616039,
        0.04208951276557549,
        0.00428765048212909,
        0.00036521216929462,
        0.00002740064222642,
        0.00000181240233365,
        0.00000010965775866,
        0.00000000598718405,
        0.00000000030769081,
        0.00000000001431793,
        0.00000000000065109,
        0.00000000000002596,
        0.00000000000000111,
        0.00000000000000004
    };

    public static double gammaf(double x) {
        double u = x;
        double h = 0;
        
        if (u <= 0) {
            if (u == Math.floor(x)) {
                System.err.println("ARGUMENT IS NON-POSITIVE INTEGER = " + u);
                return 0;
            } else {
                u = 1 - x;
            }
        }
        
        double f = 1;
        if (u < 3) {
            for (int i = 1; i <= (int)(4 - u); i++) {
                f /= u;
                u += 1;
            }
        } else {
            for (int i = 1; i <= (int)(u - 3); i++) {
                u -= 1;
                f *= u;
            }
        }
        
        h = u + u - 7;
        double alfa = h + h;
        double b0 = 0, b1 = 0, b2 = 0;
        
        for (int i = 15; i >= 0; i--) {
            b0 = C[i] + alfa * b1 - b2;
            b2 = b1;
            b1 = b0;
        }
        
        h = f * (b0 - h * b2);
        if (x < 0) {
            h = PI / (Math.sin(PI * x) * h);
        }
        
        return h;
    }

    public static void main(String[] args) {
        // Example usage
        System.out.println("Gamma(0.5) = " + gammaf(0.5));
        System.out.println("Gamma(1.0) = " + gammaf(1.0));
        System.out.println("Gamma(5) = "   + gammaf(5.0));
    }
}
