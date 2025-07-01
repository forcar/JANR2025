package org.clas.lib;

import static java.lang.Math.*;

public class JanrIniPoint extends Constants {
	
	double q2;
	
	public JanrIniPoint(double q2) {
		janr_ini_hel(q2);
		janr_ini_point();		
	}

    public void janr_ini_point() {

        double[] cc = {-0.57735, 0.0, 1.2247};

        double a12, a32, s12, gm;

        // P11(1440) M1-
        a32 = 0.0;
        a12 = cc[ires[21] - 1] * aa1[21] * cmp1;
        s12 = cc[ires[21] - 1] * sa1[21] * csp1;
        atomm(lres[21], a12, a32, s12,
                am1[0][0],
                ae1[0][0],
                as1[0][0]);

        // S11(1535) E0+
        a32 = 0.0;
        a12 = cc[ires[22] - 1] * aa1[22] * cep2;
        s12 = cc[ires[22] - 1] * sa1[22] * csp2;
        atomp(lres[22], a12, a32, s12,
                am1[0][1],
                ae1[0][1],
                as1[0][1]);

        // D13(1520) M2- E2-
        a32 = cc[ires[23] - 1] * aa3[23] * cep3;
        a12 = cc[ires[23] - 1] * aa1[23] * cmp3;
        s12 = cc[ires[23] - 1] * sa1[23] * csp3;
        atomm(lres[23], a12, a32, s12,
                am1[0][2],
                ae1[0][2],
                as1[0][2]);

        // S11(1650) E0+
        a32 = 0.0;
        a12 = cc[ires[25] - 1] * aa1[25] * cep4;
        s12 = cc[ires[25] - 1] * sa1[25] * csp4;
        atomp(lres[25], a12, a32, s12,
                am1[0][4],
                ae1[0][4],
                as1[0][4]);

        // S31(1620)  E0+
        a32 = 0.0;
        a12 = cc[ires[5] - 1] * aa1[5] * cep5;
        s12 = cc[ires[5] - 1] * sa1[5] * csp5;
        atomp(lres[5], a12, a32, s12,
                am3[4],
                ae3[4],
                as3[4]);

        // F15(1680) M3- E3-
        a32 = cc[ires[26] - 1] * aa3[26] * cep6;
        a12 = cc[ires[26] - 1] * aa1[26] * cmp6;
        s12 = cc[ires[26] - 1] * sa1[26] * csp6;
        atomm(lres[26], a12, a32, s12,
                am1[0][5],
                ae1[0][5],
                as1[0][5]);

        // D33(1700) M2- E2-
        a32 = cc[ires[2] - 1] * aa3[2] * cep7;
        a12 = cc[ires[2] - 1] * aa1[2] * cmp7;
        s12 = cc[ires[2] - 1] * sa1[2] * csp7;
        atomm(lres[2], a12, a32, s12,
                am3[1],
                ae3[1],
                as3[1]);

        // D13(1700) M2- E2-
        a32 = cc[ires[24] - 1] * aa3[24];
        a12 = cc[ires[24] - 1] * aa1[24];
        s12 = cc[ires[24] - 1] * sa1[24];
        atomm(lres[24], a12, a32, s12,
                am1[0][3],
                ae1[0][3],
                as1[0][3]);

        // P11(1710) M1-
        a32 = cc[ires[30] - 1] * aa3[30];
        a12 = cc[ires[30] - 1] * aa1[30];
        s12 = cc[ires[30] - 1] * sa1[30];
        atomm(lres[30], a12, a32, s12,
                am1[0][9],
                ae1[0][9],
                as1[0][9]);

        // P13(1720) M1+ E1+
        a32 = cc[ires[27] - 1] * aa3[27] * cep8;
        a12 = cc[ires[27] - 1] * aa1[27] * cmp8;
        s12 = cc[ires[27] - 1] * sa1[27] * csp8;
        atomp(lres[27], a12, a32, s12,
                am1[0][6],
                ae1[0][6],
                as1[0][6]);

        // D15(1675) M2+ E2+
        a32 = cc[ires[29] - 1] * aa3[29] * cep9;
        a12 = cc[ires[29] - 1] * aa1[29] * cmp9;
        s12 = cc[ires[29] - 1] * sa1[29];
        atomp(lres[29], a12, a32, s12,
                am1[0][8],
                ae1[0][8],
                as1[0][8]);

        // P33(1600) M1+ E1+
        a32 = cc[ires[4] - 1] * aa3[4];
        a12 = cc[ires[4] - 1] * aa1[4];
        s12 = cc[ires[4] - 1] * sa1[4];
        atomp(lres[4], a12, a32, s12,
                am3[3],
                ae3[3],
                as3[3]);

        // P33(1232) M1+ E1+ (Tiator parametrization for G_M)
        gm = exp(-0.21 * q2) / (pow((1.0 + q2 / 0.71), 2.0 / 3.0));
        resM = cm1 * gm * (sqrt(pow((2.3933 + q2) / 2.46, 2.0) - 0.88) * 6.568 * 1.013);
        am3[0] = resM;
        ae3[0] = -ce1 * am3[0] * 0.025;
        as3[0] = -cs1 * am3[0] * 0.05;

        // F35(1905) M3- E3-
        a32 = cc[ires[6] - 1] * aa3[6] * cep10;
        a12 = cc[ires[6] - 1] * aa1[6] * cmp10;
        s12 = cc[ires[6] - 1] * sa1[6];
        atomm(lres[6], a12, a32, s12,
                am3[5],
                ae3[5],
                as3[5]);

        // F37(1950) M3+ E3+
        a32 = cc[ires[3] - 1] * aa3[3] * cep11;
        a12 = cc[ires[3] - 1] * aa1[3] * cmp11;
        s12 = cc[ires[3] - 1] * sa1[3];
        atomp(lres[3], a12, a32, s12,
                am3[2],
                ae3[2],
                as3[2]);

        // F17(1990)
        am1[0][7] = 0.0;
        ae1[0][7] = 0.0;
        as1[0][7] = 0.0;

        for (int i = 6; i <= 12; i++) {
            am3[i] = 0.0;
            ae3[i] = 0.0;
            as3[i] = 0.0;
        }
        for (int i = 10; i <= 12; i++) {
            am1[0][i] = 0.0;
            ae1[0][i] = 0.0;
            as1[0][i] = 0.0;
        }
        for (int i = 0; i <= 12; i++) {
            am1[1][i] = 0.0;
            ae1[1][i] = 0.0;
            as1[1][i] = 0.0;
        }
    }

    public void janr_ini_hel(double q2) {

        for (int i = 0; i < 50; i++) {
            aa3[i] = 0.0;
            aa1[i] = 0.0;
            sa1[i] = 0.0;
        }

        double dip = 1.0 / pow(1.0 + q2 / 0.71, 2);

        // P11(1440)
        aa3[21] = 0.0;
        aa1[21] = -58 * (1 - 0.78 * q2 - 3.08 * q2 * q2 - 0.0021 * pow(q2, 3)) * exp(-1.24 * q2);
        sa1[21] = 27 * (1 + 4.55 * q2 - 2.25 * q2 * q2 + 0.5100 * pow(q2, 3)) * exp(-1.01 * q2);

        // S11(1535)
        aa3[22] = 0.0;
        aa1[22] = 91 * (1 + 1.55 * q2 - 0.47 * q2 * q2 + 0.61 * pow(q2, 3)) * exp(-1.16 * q2);
        sa1[22] = -30 * (1 - 1.34 * q2 + 3.17 * q2 * q2 - 0.50 * pow(q2, 3)) * exp(-0.99 * q2);

        // D13(1520)
        aa1[23] = -28 * (1 + 8 * q2 - 3.89 * q2 * q2 + 1.15 * pow(q2, 3)) * exp(-0.96 * q2);
        aa3[23] = 146 * (1 - 0.42 * q2 + 3.5 * q2 * q2 - 0.0089 * pow(q2, 3)) * exp(-2.5 * q2);
        sa1[23] = -5 * (1 + 84 * q2 - 97 * q2 * q2 + 32 * pow(q2, 3)) * exp(-2.1 * q2);

        // S11(1650)
        aa3[25] = 0.0;
        aa1[25] = 22.0 * dip;
        sa1[25] = 10.0 * dip;

        // S31(1620)
        aa3[5] = 0.0;
        aa1[5] = 49.6 * dip;
        sa1[5] = -28.0 * dip;

        // F15(1680)
        aa1[26] = -17 * (1 + 15.7 * q2 - 19.7 * q2 * q2 + 14.3 * pow(q2, 3)) * exp(-1.93 * q2);
        aa3[26] = 133 * (1 + 1.06 * q2 - 1.42 * q2 * q2 + 1.37 * pow(q2, 3)) * exp(-2.0 * q2);
        sa1[26] = -4 * (1 + 28 * q2 - 18.8 * q2 * q2 + 4.7 * pow(q2, 3)) * exp(-1.1 * q2);

        // D33(1700)
        aa1[2] = 125.4 * dip;
        aa3[2] = 105.0 * dip;
        sa1[2] = 10.0 * dip;

        // D13(1700)
        aa1[24] = 0.0;
        aa3[24] = 0.0;
        sa1[24] = 0.0;

        // P11(1710)
        aa1[30] = 0.0;
        aa3[30] = 0.0;
        sa1[30] = 0.0;

        // P13(1720)
        aa1[27] = 96.6 * dip;
        aa3[27] = -39.0 * dip;
        sa1[27] = 10.0 * dip;

        // D15(1675)
        aa1[29] = 18.0 * dip;
        aa3[29] = 21.0 * dip;
        sa1[29] = 0.0;

        // P33(1600)
        aa1[4] = 0.0;
        aa3[4] = 0.0;
        sa1[4] = 0.0;

        // F35(1905)
        aa1[6] = 21.3 * dip;
        aa3[6] = -45.6 * dip;
        sa1[6] = 0.0;

        // F37(1955)
        aa1[3] = 10.0 * dip;
        aa3[3] = -10.0 * dip;
        sa1[3] = 0.0;
    }

    // Tiator nucl-th/0610113
    // mptoa: Fortran passes by reference; Java returns via array
    public void mptoa(double l, double m, double e, double s, double[] out) {
        // out[0] = a12, out[1] = a32, out[2] = s12
        out[0] = -0.5 * ((l + 2) * e + l * m);
        out[1] = 0.5 * sqrt(l * (l + 2)) * (e - m);
        out[2] = -(l + 1) * s / sqrt(2.0);
    }

    public void mmtoa(double ll, double m, double e, double s, double[] out) {
        double l = ll - 1;
        out[0] = 0.5 * ((l + 2) * m - l * e);
        out[1] = -0.5 * sqrt(l * (l + 2)) * (e + m);
        out[2] = -(l + 1) * s / sqrt(2.0);
    }

    public void atomp(double l, double a12, double a32, double s12, double m, double e, double s) {
        double aa32;
        if (l > 0) {
            aa32 = a32 / sqrt(l * (l + 2));
            m = -2.0 / (2 * l + 2) * (aa32 * (l + 2) + a12);
        } else {
            aa32 = 0.0;
            m = 0.0;
        }
        e = 2.0 / (2 * l + 2) * (aa32 * l - a12);
        s = -sqrt(2.0) * s12 / (l + 1);

    }

    public void atomm(double ll, double a12, double a32, double s12, double m, double e, double s) {
        double l = ll - 1;
        double aa32;
        if (l > 0) {
            aa32 = a32 / sqrt(l * (l + 2));
            e = -2.0 / (2 * l + 2) * (aa32 * (l + 2) + a12);
        } else {
            aa32 = 0.0;
            e = 0.0;
        }
        m = -2.0 / (2 * l + 2) * (aa32 * l - a12);
        s = -sqrt(2.0) * s12 / (l + 1);
    }

    // Overloads for convenience (to match Fortran call by reference)
    /*
    private void atomp(double l, double a12, double a32, double s12, DoubleSetter m, DoubleSetter e, DoubleSetter s) {
        double[] out = new double[3];
        atomp(l, a12, a32, s12, v -> m.set(v[0]), v -> e.set(v[0]), v -> s.set(v[0]));
    }

    private void atomm(double ll, double a12, double a32, double s12, DoubleSetter m, DoubleSetter e, DoubleSetter s) {
        double[] out = new double[3];
        atomm(ll, a12, a32, s12, v -> m.set(v[0]), v -> e.set(v[0]), v -> s.set(v[0]));
    }

    @FunctionalInterface
    interface DoubleSetter {
        void set(double v);
    }
    */

    public static void main(String[] args) {
    	double  q2 = 04;
    	JanrIniPoint jip = new JanrIniPoint(q2);
    }
}

