package org.clas.lib;

/******************************************************************************
 *  Compilation:  javac AdaptiveQuadrature.java
 *  Execution:    java AdaptiveQuadrature a b
 *
 *  Numerically integrate the function in the interval [a, b] using
 *  the extrapolated Simpson's rule in an adaptive recursive algorithm.
 *
 *  See https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/moler/quad.pdf
 *
 *  % java AdaptiveQuadrature -3 3          // true answer = 0.9973002040...
 *  0.9973002039366737
 *  0.9973002036280225
 *
 *  % java AdaptiveQuadrature 0 10000000    // true answer = 1/2
 *  1.9947114020071635
 *  0.500000121928651
 *
 *
 *  Note that even with 1 million subintervals, the trapezoid rule
 *  gets the integral of f(x) = exp(-x^2) / sqrt(2 pi) from 0 to 10,000,000
 *  completely wrong.
 *  In contrast, the adaptive quadrature rule is accurate to 8 significant
 *  digits. In this example, both methods evaluate f() at roughly 1 million
 *  points.
 *
 *  Note: many function calls can be saved by computing f(a), f(b), and f(c)
 *        only once per recursive call (instead of twice) and adding extra
 *        arguments to adapative() and passing the already computed values
 *        that will get reused in the next recursive call.
 *
 ******************************************************************************/


public class AdaptiveQuadrature {
    private final static double EPSILON = 1E-6;

   /***************************************************************************
    * Standard normal distribution density function.
    * Replace with any sufficiently smooth function.
    ***************************************************************************/
    static double f(double x) {
        return Math.exp(- x * x / 2) / Math.sqrt(2 * Math.PI);
    }

    // trapezoid rule
    static double trapezoid(double a, double b, int n) {
        double h = (b - a) / n;
        double sum = 0.5 *  h * (f(a) + f(b));
        for (int k = 1; k < n; k++)
            sum = sum + h * f(a + h*k);
        return sum;
    }


    // adaptive quadrature
    static double adaptive(double a, double b) {
        double h = b - a;
        double c = (a + b) / 2.0;
        double d = (a + c) / 2.0;
        double e = (b + c) / 2.0;
        double Q1 = h/6  * (f(a) + 4*f(c) + f(b));
        double Q2 = h/12 * (f(a) + 4*f(d) + 2*f(c) + 4*f(e) + f(b));
        if (Math.abs(Q2 - Q1) <= EPSILON)
            return Q2 + (Q2 - Q1) / 15;
        else
            return adaptive(a, c) + adaptive(c, b);
    }


    // sample client program
    public static void main(String[] args) {
        double a = 0;
        double b = 20;
        System.out.println(trapezoid(a, b, 1000000));
        System.out.println(adaptive(a, b));
    }

}
