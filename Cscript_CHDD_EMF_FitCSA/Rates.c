#include "math.h"
#include "Rates.h"


double GammaC = 6.72828e7;
double GammaH = 2.67522e8;
double GammaD = 4.1065e7;


double dCH = -136935.20560730906;
double dCD = -21019.744986446523;
double dHD = -19192.023064957146;


double JCH(double w, double *Xarr, double Tauc)
{
    	double Sf2 = Xarr[0];
    	double Ss2 = Xarr[1];
    	double tm = Xarr[2];
    	double tf = Xarr[3];
    	double ts = Xarr[4];
    	
    	double tc = Tauc;
    	
    	double tfp = tc*tf/(tc + tf);
    	double tsp = tc*ts/(tc + ts);
    	
    	double tcm = tc*tm/(tc + tm);
    	double tfm = tc*tm*tf/(tc*tm + tm*tf + tc*tf);
    	double tsm = tc*tm*ts/(tc*tm + tm*ts + tc*ts);
    	
    	double Sm2 = pow((3.*pow(cos(109.47*M_PI/180.0), 2) - 1.)/2., 2);
    	
    	
    	double Jch = 2./5. * Sm2 * (Sf2*Ss2*tc/(1. + pow(w*tc, 2)) + (1. - Sf2) * tfp/(1. + pow(w*tfp, 2)) + Sf2*(1. - Ss2)*tsp/(1. + pow(w*tsp, 2))) + \
                    	2./5. * (1. - Sm2) * (Sf2*Ss2*tcm/(1. + pow(w*tcm, 2)) + (1. - Sf2) * tfm/(1. + pow(w*tfm, 2)) + Sf2*(1. - Ss2)*tsm/(1. + pow(w*tsm, 2)));

        return Jch;

}

double JHD(double w, double *Xarr, double Tauc)
{
    	double Sf2 = Xarr[0];
    	double Ss2 = Xarr[1];
    	double tm = Xarr[2];
    	double tf = Xarr[3];
    	double ts = Xarr[4];
    	
    	double tc = Tauc;
    	
    	double tfp = tc*tf/(tc + tf);
    	double tsp = tc*ts/(tc + ts);
    	
    	double tcm = tc*tm/(tc + tm);
    	double tfm = tc*tm*tf/(tc*tm + tm*tf + tc*tf);
    	double tsm = tc*tm*ts/(tc*tm + tm*ts + tc*tm);
    	
    	double Jhd = 2./5. * 1./4. * (Sf2*Ss2*tc/(1. + pow(w*tc, 2)) + (1. - Sf2) * tfp/(1. + pow(w*tfp, 2)) + Sf2*(1. - Ss2)*tsp/(1. + pow(w*tsp, 2))) + \
                    	2./5. * 3./4. * (Sf2*Ss2*tcm/(1. + pow(w*tcm, 2)) + (1. - Sf2) * tfm/(1. + pow(w*tfm, 2)) + Sf2*(1. - Ss2)*tsm/(1. + pow(w*tsm, 2)));

        return Jhd;

}


double JCC(double w, double *Xarr, double Tauc)
{
    	double Sf2 = Xarr[0];
    	double Ss2 = Xarr[1];
    	double tf = Xarr[3];
    	double ts = Xarr[4];
    	
    	double tc = Tauc;
    	
    	double tfp = tc*tf/(tc + tf);
    	double tsp = tc*ts/(tc + ts);
    	
    	double Jcc = 2./5. * (Sf2*Ss2*tc/(1. + pow(w*tc, 2)) + (1. - Sf2) * tfp/(1. + pow(w*tfp, 2)) + Sf2*(1. - Ss2)*tsp/(1. + pow(w*tsp, 2)));

        return Jcc;

}


double JCCH(double w, double *Xarr, double Tauc)
{
    	double Sf2 = Xarr[0];
    	double Ss2 = Xarr[1];
    	double tf = Xarr[3];
    	double ts = Xarr[4];
    	
    	double tc = Tauc;
    	
    	double tfp = tc*tf/(tc + tf);
    	double tsp = tc*ts/(tc + ts);
    	
    	double Sm2 = (3.*pow(cos(109.47*M_PI/180.0), 2) - 1.)/2.;
    	
    	double Jcc = 2./5. * Sm2 * (Sf2*Ss2*tc/(1. + pow(w*tc, 2)) + (1. - Sf2) * tfp/(1. + pow(w*tfp, 2)) + Sf2*(1. - Ss2)*tsp/(1. + pow(w*tsp, 2)));

        return Jcc;

}





struct mat2d RelaxMatrix( double B, double *Xarr, double Tauc, double *OtherInputs )
{
    	double CSAValueInPPMofC = Xarr[5]*pow(10, -6);
    	
    	double rxyCDvic = OtherInputs[0];
    	double rzCDvic = OtherInputs[1];
    	
    	double rCDvic = pow(pow(rxyCDvic,2) + pow(rzCDvic, 2), 0.5);
    	double dCDvic = dCD/pow(rCDvic, 3) * pow(1.115e-10, 3);
    	
    	double rHDvic = pow(pow(rxyCDvic - 1.115e-10*sin(109.47*M_PI/180.),2) + pow(rzCDvic - 1.115e-10*cos(109.47*M_PI/180.), 2), 0.5);
    	double dHDvic = dHD/pow(rHDvic, 3) * pow(1.8208e-10, 3);
    	
    	double tc = Tauc;
    	
    	double JCHm = JCH(B*(GammaC - GammaH), Xarr, tc);
    	double JCHwc = JCH(B*GammaC, Xarr, tc);
    	double JCHwh = JCH(B*GammaH, Xarr, tc);
    	double JCHp = JCH(B*(GammaC + GammaH), Xarr, tc);
    	
    	double JHDm = JHD(B*(GammaD - GammaH), Xarr, tc);
    	double JHDwh = JHD(B*GammaH, Xarr, tc);
    	double JHDp = JHD(B*(GammaD + GammaH), Xarr, tc);
    	
    	double JCDm = JCH(B*(GammaC - GammaD), Xarr, tc);
    	double JCDwc = JCH(B*GammaC, Xarr, tc);
    	double JCDp = JCH(B*(GammaC + GammaD), Xarr, tc);
    	
    	double JCCwc = JCC(B*GammaC, Xarr, tc);
    	
    	double JCCHwc = JCCH(B*GammaC, Xarr, tc);
    	
    	
    	
    	
    	
        struct mat2d Mat;
    	
    	double R00 = 1./3. * pow(CSAValueInPPMofC*B*GammaC, 2) * JCCwc + \
                    	1./4. * pow(dCH, 2) * (JCHm + 3.*JCHwc + 6.*JCHp) + \
                    4./3. * pow(dCD, 2) * (JCDm + 3.*JCHwc + 6.*JCDp) + \
                    2./3. * pow(dCDvic, 2) * (JCDm + 3.*JCHwc + 6.*JCDp);
                    
    	
    	double R01 = 1./4. * pow(dCH, 2) * (-JCHm + 6. * JCHp);
    	
    	double R02 = CSAValueInPPMofC * B * GammaC * dCH * JCCHwc;
    	
    	
    	double R11 = 1./4. * pow(dCH, 2) * (JCHm + 3.*JCHwh + 6.*JCHp) + \
                    	4./3. * pow(dHD, 2) * (JHDm + 3.*JHDwh + 6.*JHDp) + \
                    	2./3. * pow(dHDvic, 2) * (JHDm + 3.*JHDwh + 6.*JHDp);
    	
    	double R12 = 0.0;
    	
    	
    	double R22 = 1./3. * pow(CSAValueInPPMofC*B*GammaC, 2) * JCCwc + 3./4. * pow(dCH, 2) * (JCHwc + JCHwh) + \
                    	4./3. * pow(dCD, 2) * (JCDm + 3. * JCDwc + 6. * JCDp) + \
                        4./3. * pow(dHD, 2) * (JHDm + 3. * JHDwh + 6. * JHDp) + \
                        2./3. * pow(dCDvic, 2) * (JCDm + 3. * JCHwc + 6. * JCDp) + \
                        2./3. * pow(dHDvic, 2) * (JHDm + 3. * JHDwh + 6. * JHDp);
        	
    	
    	Mat.m[0][0] = R00;
    	Mat.m[0][1] = R01;
    	Mat.m[0][2] = 0.5049*R02;
    	
    	Mat.m[1][0] = R01;
    	Mat.m[1][1] = R11;
    	Mat.m[1][2] = R12;
    	
    	Mat.m[2][0] = 0.5049*R02;
    	Mat.m[2][1] = R12;
    	Mat.m[2][2] = R22;
    	
    	return Mat;

        
}


double R2calculation(double B, double *Xarr, double Tauc, double *OtherInputs )
{
    	double CSAValueInPPMofC = Xarr[5]*pow(10, -6);
    	
    	double rxyCDvic = OtherInputs[0];
    	double rzCDvic = OtherInputs[1];
    	
    	double rCDvic = pow(pow(rxyCDvic,2) + pow(rzCDvic, 2), 0.5);
    	double dCDvic = dCD/pow(rCDvic, 3) * pow(1.115e-10, 3);
    	
    	double tc = Tauc;
    	
    	double JCH0 = JCH(0., Xarr, tc);
    	double JCHm = JCH(B*(GammaC - GammaH), Xarr, tc);
    	double JCHwc = JCH(B*GammaC, Xarr, tc);
    	double JCHwh = JCH(B*GammaH, Xarr, tc);
    	double JCHp = JCH(B*(GammaC + GammaH), Xarr, tc);
    	
    	double JCDm = JCH(B*(GammaC - GammaD), Xarr, tc);
    	double JCDwd = JCH(B*GammaD, Xarr, tc);
    	double JCDp = JCH(B*(GammaC + GammaD), Xarr, tc);
    	
    	double JCC0 = JCC(0., Xarr, tc);
    	double JCCwc = JCC(B*GammaC, Xarr, tc);
    	

    	double RateCalc = 1./18. * pow(CSAValueInPPMofC*B*GammaC, 2) * (4. * JCC0 + 3. * JCCwc) + \
                        	1./8. * pow(dCH, 2) * (4. * JCH0 + JCHm + 3. * JCHwc + 6. * JCHwh + 6. * JCHp) + \
                        2./3. * pow(dCD, 2) * (4. * JCH0 + JCDm + 3. * JCHwc + 6. * JCDwd + 6. * JCDp) + \
                        1./3. * pow(dCDvic, 2) * (4. * JCH0 + JCDm + 3. * JCHwc + 6. * JCDwd + 6. * JCDp);

    	return RateCalc;
	
}

double R1calculation(double B, double *Xarr, double Tauc, double *OtherInputs )
{
    	double CSAValueInPPMofC = Xarr[5]*pow(10, -6);
    	
    	double rxyCDvic = OtherInputs[0];
    	double rzCDvic = OtherInputs[1];
    	
    	double rCDvic = pow(pow(rxyCDvic,2) + pow(rzCDvic, 2), 0.5);
    	double dCDvic = dCD/pow(rCDvic, 3) * pow(1.115e-10, 3);
    	
    	double tc = Tauc;
    	
    	double JCHm = JCH(B*(GammaC - GammaH), Xarr, tc);
    	double JCHwc = JCH(B*GammaC, Xarr, tc);
    	double JCHp = JCH(B*(GammaC + GammaH), Xarr, tc);
    	
    	double JCDm = JCH(B*(GammaC - GammaD), Xarr, tc);
    	double JCDp = JCH(B*(GammaC + GammaD), Xarr, tc);
    	
    	double JCCwc = JCC(B*GammaC, Xarr, tc);
    	

    	double RateCalc = 1./3. * pow(CSAValueInPPMofC*B*GammaC, 2) * JCCwc + \
                    	1./4. * pow(dCH, 2) * (JCHm + 3.*JCHwc + 6.*JCHp) + \
                    4./3. * pow(dCD, 2) * (JCDm + 3.*JCHwc + 6.*JCDp) + \
                    2./3. * pow(dCDvic, 2) * (JCDm + 3.*JCHwc + 6.*JCDp);

    	return RateCalc;

}

double Sigmacalculation(double B, double *Xarr, double Tauc, double *OtherInputs )
{
    	double tc = Tauc;
    	
    	double JCHm = JCH(B*(GammaC - GammaH), Xarr, tc);
    	double JCHp = JCH(B*(GammaC + GammaH), Xarr, tc);
    	
    	double RateCalc = 1./4. * pow(dCH, 2) * (-JCHm + 6. * JCHp);

    	return RateCalc;

}

