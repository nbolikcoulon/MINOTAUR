#include <complex.h>
#include "math.h"
#include "Rates.h"

double DiffCoeff[] = {33897446.28099173, 39850016.52892562};

double WignerAniso[] = {0.8635283767131133, 0.4765001025632064, 0.16101451089182553, 0.03627235696036371, 0.005003825663960973,
                       -0.4765001025632064, 0.6663266802787249, 0.5391666736713626, 0.19219787077042755, 0.03627235696036371,
                       0.16101451089182553, -0.5391666736713626, 0.6055966071312231, 0.5391666736713626, 0.16101451089182553,
                       -0.03627235696036371, 0.19219787077042755, -0.5391666736713626, 0.6663266802787249, 0.4765001025632064,
                       0.005003825663960973, -0.03627235696036371, 0.16101451089182553, -0.4765001025632064, 0.8635283767131133};

double WignerJump[] = {0.9962491183800163, 0.08640951145062285, 0.0045895571608420204, 0.0001625132432685409, 3.5238901836935666e-06,
                      -0.08640951145062285, 0.9906280817853166, 0.10563056872737676, 0.005617512704516113, 0.0001625132432685409,
                      0.0045895571608420204, -0.10563056872737676, 0.9887579268106004, 0.10563056872737676, 0.0045895571608420204,
                      -0.0001625132432685409, 0.005617512704516113, -0.10563056872737676, 0.9906280817853166, 0.08640951145062285,
                      3.5238901836935666e-06, -0.0001625132432685409, 0.0045895571608420204, -0.08640951145062285, 0.9962491183800163,
                      0.1676726221373039, 0.4027119145518626, 0.5923009582731981, 0.5807632577939187, 0.34871567055241387,
                      -0.4027119145518626,-0.5577449388281192, -0.2180674694800916, 0.3767018904130091, 0.5807632577939187,
                      0.5923009582731981, 0.2180674694800916, -0.4508351219308461, -0.2180674694800916, 0.5923009582731981,
                      -0.5807632577939187, 0.3767018904130091, 0.2180674694800916, -0.5577449388281192, 0.4027119145518626,
                      0.34871567055241387, -0.5807632577939187, 0.5923009582731981, -0.4027119145518626, 0.1676726221373039};
                      
double WignerCH[] = {0.5443392536840261, -0.3848798852734992, -0.33335341849327244, 0.3848798852734992, 0.5443392536840261};
double WignerHD[] = {0.0, 0.0, 1.0, 0.0, 0.0};

double WignerCSAp[] = {0.1047079207020818, 0.46111385049344367, 0.7435190222520962, -0.46111385049344367, 0.1047079207020818,
                      0.46054092252378687, -0.528865295304651, -0.12809026585391814, 0.528865295304651, 0.46054092252378687};
double WignerCSAo[] = {0.5378782808229333, 0.4003437920306589, -0.3175273317416253, -0.4003437920306589, 0.5378782808229333,
                      0.25441895248574575, -0.6035566261318802, 0.3768033855165249, 0.6035566261318802, 0.25441895248574575};
                  
double AlphaJump[] = {3.6887203750191313, 4.211993147641056};
double AlphaCSAp[] = {3.012466101010791, 2.419476536899442};
double AlphaCSAo[] = {4.003202905473472, 3.847238202648743};
double bCH = 1.9106119321581925;

double DipolarCst[] = {-136935.20560730906, -21019.744986446523, -19192.023064957146};
double CSAp[] = {25.578157181996232, 14.16505408206902};
double CSAo[] = {8.251283445795252, 2.55845109902765};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double AutoDD(double w, double DifCoeff[2], double SmallWignAniso[25], double Drot, double pop[2], double KineticEigVal[2], double KineticEigVec[4], double AlphaJump[2], double SmallWignJump[50], double AmplitudeI, double SmallWignI[5])
{	
	int a;
	int b;
	int bp;
	int c;
	int alpha;
	int beta;
	int n;
	
    double Spec = 0.0;
    
	for (a = 0; a < 5; ++a) {
    	for (b = 0; b < 5; ++b) {
            	for (bp = 0; bp < 5; ++ bp) {
                	for (c = 0; c < 5; ++c) {
                    	for (alpha = 0; alpha < 2; ++ alpha) {
                        	for (beta = 0; beta < 2; ++ beta) {
                            	double WignerProd = SmallWignAniso[a*5 + b]*SmallWignAniso[a*5 + bp] * cos((b-2)*AlphaJump[alpha] - (bp-2)*AlphaJump[beta]) * SmallWignJump[alpha*25 + b*5 + c] * SmallWignJump[beta*25 + bp*5 + c] * pow(SmallWignI[c], 2);
                            	
                            	for (n = 0; n < 2; ++ n) {
                                    double tauEffInv = 6.*DifCoeff[0] + pow(a-2, 2) * (DifCoeff[1] - DifCoeff[0]) - KineticEigVal[n] + pow(c-2,2) * Drot;
                                    double tauEff = 1./tauEffInv;
                                    double Lorentz = tauEff / (1. + pow(w*tauEff,2));
                                    
                                    Spec += pow(AmplitudeI, 2)/5. * WignerProd * Lorentz *
                                            pow(pop[alpha]*pop[beta], 0.5) * KineticEigVec[n*2 + alpha] * KineticEigVec[n*2 + beta];
    }}}}}}}
                        
    return Spec;
}


double AutoCSA(double w, double DifCoeff[2], double SmallWignAniso[25], double pop[2], double KineticEigVal[2], double KineticEigVec[4], double AlphaJump[2], double SmallWignJump[50], double AmplitudeI[2], double AlphaI[2], double SmallWignI[10])
{	
	int a;
	int b;
	int bp;
	int c;
	int cp;
	int alpha;
	int beta;
	int n;
	
    double Spec = 0.0;
    
	for (a = 0; a < 5; ++a) {
    	for (b = 0; b < 5; ++b) {
            	for (bp = 0; bp < 5; ++ bp) {
                	for (c = 0; c < 5; ++c) {
                    	for (cp = 0; cp < 5; ++cp) {
                        	for (alpha = 0; alpha < 2; ++ alpha) {
                            	for (beta = 0; beta < 2; ++ beta) {
                                	double WignerProd = SmallWignAniso[a*5 + b]*SmallWignAniso[a*5 + bp] * cos((b-2)*AlphaJump[alpha] - (bp-2)*AlphaJump[beta] + (c-2)*AlphaI[alpha] - (cp-2)*AlphaI[beta]) * SmallWignJump[alpha*25 + b*5 + c] * SmallWignJump[beta*25 + bp*5 + cp] * SmallWignI[alpha*5 + c] * SmallWignI[beta*5 + cp];
                                	
                                	for (n = 0; n < 2; ++ n) {
                                        double tauEffInv = 6.*DifCoeff[0] + pow(a-2, 2) * (DifCoeff[1] - DifCoeff[0]) - KineticEigVal[n];
                                        double tauEff = 1./tauEffInv;
                                        double Lorentz = tauEff / (1. + pow(w*tauEff,2));
                                            
                                        Spec += AmplitudeI[alpha]*AmplitudeI[beta]/5. * WignerProd * Lorentz *
                                                pow(pop[alpha]*pop[beta], 0.5) * KineticEigVec[n*2 + alpha] * KineticEigVec[n*2 + beta];
    }}}}}}}}
                        
    return Spec;
}

double CrossCSA(double w, double DifCoeff[2], double SmallWignAniso[25], double pop[2], double KineticEigVal[2], double KineticEigVec[4], double AlphaJump[2], double SmallWignJump[50], double AmplitudeCSAp[2], double AmplitudeCSAo[2], double AlphaCSAp[2], double AlphaCSAo[2], double SmallWignCSAp[10], double SmallWignCSAo[10])
{	
	int a;
	int b;
	int bp;
	int c;
	int cp;
	int alpha;
	int beta;
	int n;
	
    double Spec = 0.0;
    
	for (a = 0; a < 5; ++a) {
    	for (b = 0; b < 5; ++b) {
            	for (bp = 0; bp < 5; ++ bp) {
                	for (c = 0; c < 5; ++c) {
                    	for (cp = 0; cp < 5; ++cp) {
                        	for (alpha = 0; alpha < 2; ++ alpha) {
                            	for (beta = 0; beta < 2; ++ beta) {
                                	double WignerProd = SmallWignAniso[a*5 + b]*SmallWignAniso[a*5 + bp] * cos((b-2)*AlphaJump[alpha] - (bp-2)*AlphaJump[beta] + (c-2)*AlphaCSAp[alpha] - (cp-2)*AlphaCSAo[beta]) * SmallWignJump[alpha*25 + b*5 + c] * SmallWignJump[beta*25 + bp*5 + cp] * SmallWignCSAp[alpha*5 + c] * SmallWignCSAo[beta*5 + cp];
                                	
                                	for (n = 0; n < 2; ++ n) {
                                        double tauEffInv = 6.*DifCoeff[0] + pow(a-2, 2) * (DifCoeff[1] - DifCoeff[0]) - KineticEigVal[n];
                                        double tauEff = 1./tauEffInv;
                                        double Lorentz = tauEff / (1. + pow(w*tauEff,2));
                                            
                                        Spec += AmplitudeCSAp[alpha]*AmplitudeCSAo[beta]/5. * WignerProd * Lorentz *
                                                pow(pop[alpha]*pop[beta], 0.5) * KineticEigVec[n*2 + alpha] * KineticEigVec[n*2 + beta];
    }}}}}}}}
                        
    return Spec;
}

double CrossCSADD(double w, double DifCoeff[2], double SmallWignAniso[25], double pop[2], double KineticEigVal[2], double KineticEigVec[4], double AlphaJump[2], double SmallWignJump[50], double AmplitudeI[2], double AlphaI[2], double SmallWignI[10], double AmplitudeJ, double BetaJ)
{	
	int a;
	int b;
	int bp;
	int c;
	int alpha;
	int beta;
	int n;
	
    double Spec = 0.0;
    
    double LegendreDD = (3.*pow(cos(BetaJ), 2) - 1.)/2.;
    
	for (a = 0; a < 5; ++a) {
    	for (b = 0; b < 5; ++b) {
            	for (bp = 0; bp < 5; ++ bp) {
                	for (c = 0; c < 5; ++c) {
                    	for (alpha = 0; alpha < 2; ++ alpha) {
                        	for (beta = 0; beta < 2; ++ beta) {
                            	double WignerProd = SmallWignAniso[a*5 + b]*SmallWignAniso[a*5 + bp] * cos((b-2)*AlphaJump[alpha] - (bp-2)*AlphaJump[beta] + (c-2)*AlphaI[alpha]) * SmallWignJump[alpha*25 + b*5 + c] * SmallWignJump[beta*25 + bp*5 + 2] * SmallWignI[alpha*5 + c];
                            	
                            	for (n = 0; n < 2; ++ n) {
                                    double tauEffInv = 6.*DifCoeff[0] + pow(a-2, 2) * (DifCoeff[1] - DifCoeff[0]) - KineticEigVal[n];
                                    double tauEff = 1./tauEffInv;
                                    double Lorentz = tauEff / (1. + pow(w*tauEff,2));
                                        
                                    Spec += AmplitudeI[alpha]*AmplitudeJ/5. * WignerProd * LegendreDD * Lorentz *
                                            pow(pop[alpha]*pop[beta], 0.5) * KineticEigVec[n*2 + alpha] * KineticEigVec[n*2 + beta];
    }}}}}}}
                        
    return Spec;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct mat2d RelaxMatrix(double B, double *Xarr, double *TaucArr, double *OtherInputs) {
    struct mat2d Mat;
    
    double w0_c = 6.72828e7 * B;
    double w0_h = 2.67522e8 * B;
    double w0_d = 4.1065e7 * B;
    
    double p1 = Xarr[0];
    double Drot = Xarr[1];
    double k21 = pow(10, Xarr[2]);
    double CSAfactor = Xarr[3] * w0_c * pow(2./3., 0.5)/1000000;
    
    double pop[] = {p1, 1. - p1};
    double EigVal[] = {0.0, -k21/pop[1]};
    double EigVec[] = {pow(p1, 0.5), pow(pop[1], 0.5),
                      pow(pop[1], 0.5), -pow(p1, 0.5)};
    
    double JDD_CH_wc = AutoDD(w0_c, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_wh = AutoDD(w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_m = AutoDD(w0_c - w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_p = AutoDD(w0_c + w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    
    double JDD_CD_wc = AutoDD(w0_c, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    double JDD_CD_m = AutoDD(w0_c - w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    double JDD_CD_p = AutoDD(w0_c + w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    
    double JDD_HD_wh = AutoDD(w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[2], WignerHD);
    double JDD_HD_m = AutoDD(w0_h - w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[2], WignerHD);
    double JDD_HD_p = AutoDD(w0_h + w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[2], WignerHD);
    
    double JCSAp_wc = AutoCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, AlphaCSAp, WignerCSAp);
    double JCSAo_wc = AutoCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAo, AlphaCSAo, WignerCSAo);
    double JCSApo_wc = CrossCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, CSAo, AlphaCSAp, AlphaCSAo, WignerCSAp, WignerCSAo);
    
    double JCSApDD_wc = CrossCSADD(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, AlphaCSAp, WignerCSAp, DipolarCst[0], bCH);
    double JCSAoDD_wc = CrossCSADD(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAo, AlphaCSAo, WignerCSAo, DipolarCst[0], bCH);
    
    
    double R00 = 1./2. * (JDD_CH_m + 3.*JDD_CH_wc + 6.*JDD_CH_p) +
                  8./3. * (JDD_CD_m + 3.*JDD_CD_wc + 6.*JDD_CD_p) +
                  pow(CSAfactor, 2) * (2.*JCSApo_wc + JCSAp_wc + JCSAo_wc);
                 
    double R01 = 1./2. * (-JDD_CH_m + 6.*JDD_CH_p);
    
    double R02 = 2.*CSAfactor*pow(3./2., 0.5) * (JCSApDD_wc + JCSAoDD_wc);
    
    double R11 = 1./2. * (JDD_CH_m + 3.*JDD_CH_wh + 6.*JDD_CH_p) +
                 8./3. * (JDD_HD_m + 3.*JDD_HD_wh + 6.*JDD_HD_p);
    
    double R22 = 3./2. * (JDD_CH_wc + JDD_CH_wh) +
                 8./3. * (JDD_HD_m + 3.*JDD_HD_wh + 6.*JDD_HD_p) +
                 8./3. * (JDD_CD_m + 3.*JDD_CD_wc + 6.*JDD_CD_p) +
                 pow(CSAfactor, 2) * (2.*JCSApo_wc + JCSAp_wc + JCSAo_wc);
    
	Mat.m[0][0] = R00;
	Mat.m[0][1] = R01;
	Mat.m[0][2] = R02;
	
	Mat.m[1][0] = R01;
	Mat.m[1][1] = R11;
	Mat.m[1][2] = 0.0;
	
	Mat.m[2][0] = R02;
	Mat.m[2][1] = 0.0;
	Mat.m[2][2] = R22;
	
	return Mat;
}

double R1calculation(double B, double *Xarr, double *TaucArr, double *OtherInputs) {
    double w0_c = 6.72828e7 * B;
    double w0_h = 2.67522e8 * B;
    double w0_d = 4.1065e7 * B;
    
    double p1 = Xarr[0];
    double Drot = Xarr[1];
    double k21 = pow(10, Xarr[2]);
    double CSAfactor = Xarr[3] * w0_c * pow(2./3., 0.5)/1000000;
    
    double pop[] = {p1, 1. - p1};
    double EigVal[] = {0.0, -k21/pop[1]};
    double EigVec[] = {pow(p1, 0.5), pow(pop[1], 0.5),
                      pow(pop[1], 0.5), -pow(p1, 0.5)};
    
    double JDD_CH_wc = AutoDD(w0_c, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_m = AutoDD(w0_c - w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_p = AutoDD(w0_c + w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    
    double JDD_CD_wc = AutoDD(w0_c, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    double JDD_CD_m = AutoDD(w0_c - w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    double JDD_CD_p = AutoDD(w0_c + w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    
    double JCSAp_wc = AutoCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, AlphaCSAp, WignerCSAp);
    double JCSAo_wc = AutoCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAo, AlphaCSAo, WignerCSAo);
    double JCSApo_wc = CrossCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, CSAo, AlphaCSAp, AlphaCSAo, WignerCSAp, WignerCSAo);
    
    double RateCalc = 1./2. * (JDD_CH_m + 3.*JDD_CH_wc + 6.*JDD_CH_p) +
                 8./3. * (JDD_CD_m + 3.*JDD_CD_wc + 6.*JDD_CD_p) +
                 pow(CSAfactor, 2) * (2.*JCSApo_wc + JCSAp_wc + JCSAo_wc);

	return RateCalc;
}

double R2calculation(double B, double *Xarr, double *TaucArr, double *OtherInputs ) {
    double w0_c = 6.72828e7 * B;
    double w0_h = 2.67522e8 * B;
    double w0_d = 4.1065e7 * B;
    
    double p1 = Xarr[0];
    double Drot = Xarr[1];
    double k21 = pow(10, Xarr[2]);
    double CSAfactor = Xarr[3] * w0_c * pow(2./3., 0.5)/1000000;
    
    double pop[] = {p1, 1. - p1};
    double EigVal[] = {0.0, -k21/pop[1]};
    double EigVec[] = {pow(p1, 0.5), pow(pop[1], 0.5),
                      pow(pop[1], 0.5), -pow(p1, 0.5)};
    
    double JDD_CH_0 = AutoDD(0.0, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_wc = AutoDD(w0_c, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_wh = AutoDD(w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_m = AutoDD(w0_c - w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_p = AutoDD(w0_c + w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    
    double JDD_CD_0 = AutoDD(0.0, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    double JDD_CD_wd = AutoDD(w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    double JDD_CD_wc = AutoDD(w0_c, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    double JDD_CD_m = AutoDD(w0_c - w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    double JDD_CD_p = AutoDD(w0_c + w0_d, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[1], WignerCH);
    
    double JCSAp_0 = AutoCSA(0.0, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, AlphaCSAp, WignerCSAp);
    double JCSAo_0 = AutoCSA(0.0, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAo, AlphaCSAo, WignerCSAo);
    double JCSApo_0 = CrossCSA(0.0, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, CSAo, AlphaCSAp, AlphaCSAo, WignerCSAp, WignerCSAo);
    
    double JCSAp_wc = AutoCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, AlphaCSAp, WignerCSAp);
    double JCSAo_wc = AutoCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAo, AlphaCSAo, WignerCSAo);
    double JCSApo_wc = CrossCSA(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, CSAo, AlphaCSAp, AlphaCSAo, WignerCSAp, WignerCSAo);

    double RateCalc = 1./4. * (4.*JDD_CH_0 + JDD_CH_m + 3.*JDD_CH_wc + 6.*JDD_CH_wh + 6.*JDD_CH_p) +
                      4./3. * (4.*JDD_CD_0 + JDD_CD_m + 3.*JDD_CD_wc + 6.*JDD_CD_wd + 6.*JDD_CD_p) +
                      pow(CSAfactor, 2) * 1./6. * (8. * JCSApo_0 + 6. * JCSApo_wc + 4. * JCSAp_0 + 3. * JCSAp_wc + 4. * JCSAo_0 + 3. * JCSAo_wc);
            
	return RateCalc;
}

double Sigmacalculation(double B, double *Xarr, double *TaucArr, double *OtherInputs ) {
    double w0_c = 6.72828e7 * B;
    double w0_h = 2.67522e8 * B;
    
    double p1 = Xarr[0];
    double Drot = Xarr[1];
    double k21 = pow(10, Xarr[2]);
    
    double pop[] = {p1, 1. - p1};
    double EigVal[] = {0.0, -k21/pop[1]};
    double EigVec[] = {pow(p1, 0.5), pow(pop[1], 0.5),
                      pow(pop[1], 0.5), -pow(p1, 0.5)};
    
    double JDD_CH_m = AutoDD(w0_c - w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
    double JDD_CH_p = AutoDD(w0_c + w0_h, DiffCoeff, WignerAniso, Drot, pop, EigVal, EigVec, AlphaJump, WignerJump, DipolarCst[0], WignerCH);
	
	double RateCalc = 1./2. * (-JDD_CH_m + 6.*JDD_CH_p);

	return RateCalc;
}

double EtaXYcalculation(double B, double *Xarr, double *TaucArr, double *OtherInputs ) {
    double w0_c = 6.72828e7 * B;
    
    double p1 = Xarr[0];
    double k21 = pow(10, Xarr[2]);
    double CSAfactor = Xarr[3] * w0_c * pow(2./3., 0.5)/1000000;
    
    double pop[] = {p1, 1. - p1};
    double EigVal[] = {0.0, -k21/pop[1]};
    double EigVec[] = {pow(p1, 0.5), pow(pop[1], 0.5),
                      pow(pop[1], 0.5), -pow(p1, 0.5)};

    double JCSApDD_0 = CrossCSADD(0.0, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, AlphaCSAp, WignerCSAp, DipolarCst[0], bCH);
    double JCSAoDD_0 = CrossCSADD(0.0, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAo, AlphaCSAo, WignerCSAo, DipolarCst[0], bCH);
   
    double JCSApDD_wc = CrossCSADD(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, AlphaCSAp, WignerCSAp, DipolarCst[0], bCH);
    double JCSAoDD_wc = CrossCSADD(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAo, AlphaCSAo, WignerCSAo, DipolarCst[0], bCH);
    
    double RateCalc = pow(1./6., 0.5) * CSAfactor * (4. * JCSApDD_0 + 3. * JCSApDD_wc + 4. * JCSAoDD_0 + 3. * JCSAoDD_wc);
    
	return RateCalc;
}

double EtaZcalculation(double B, double *Xarr, double *TaucArr, double *OtherInputs ) {
    double w0_c = 6.72828e7 * B;
    
    double p1 = Xarr[0];
    double k21 = pow(10, Xarr[2]);
    double CSAfactor = Xarr[3] * w0_c * pow(2./3., 0.5)/1000000;
    
    double pop[] = {p1, 1. - p1};
    double EigVal[] = {0.0, -k21/pop[1]};
    double EigVec[] = {pow(p1, 0.5), pow(pop[1], 0.5),
                      pow(pop[1], 0.5), -pow(p1, 0.5)};
    
    double JCSApDD_wc = CrossCSADD(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAp, AlphaCSAp, WignerCSAp, DipolarCst[0], bCH);
    double JCSAoDD_wc = CrossCSADD(w0_c, DiffCoeff, WignerAniso, pop, EigVal, EigVec, AlphaJump, WignerJump, CSAo, AlphaCSAo, WignerCSAo, DipolarCst[0], bCH);
    
    double RateCalc = 2.*CSAfactor*pow(3./2., 0.5) * (JCSApDD_wc + JCSAoDD_wc);
    
	return RateCalc;
}




