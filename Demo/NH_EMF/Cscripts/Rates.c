#include "math.h"
#include "Rates.h"


///////////////////////////////// CONSTANTS /////////////////////////////////
//Constants
const double hbar = 1.05457173e-34;
const double mu = 4.0e-7 * M_PI;

//Gyromagnetic ratio
const double GammaH = 42.576e6 * 2.0 * M_PI;
const double GammaN = -4.316e6 * 2.0 * M_PI;

//Distances
const double rhn = 1.02e-10;
const double rhh = 2.1e-10;

//Dipolar coefficients
const double delta_NH = (mu / (4.0 * M_PI)) * hbar * GammaH * GammaN / (rhn*rhn*rhn);
const double delta_HH = (mu / (4.0 * M_PI)) * hbar * GammaH * GammaH/ (rhh*rhh*rhh);

//CSA
const double DELTA_N = -164.0e-6;///pow(3.0, 0.5);
const double sigma_Hxx = 14.6e-6;///pow(3.0, 0.5);
const double sigma_Hyy = 8.2e-6;///pow(3.0, 0.5);
const double sigma_Hzz = 2.1e-6;///pow(3.0, 0.5);
const double DELTA_H = 1.0/2.0*(sigma_Hxx + sigma_Hyy - 2.0*sigma_Hzz); //CSA of neighboring protons

//CSA Angles
const double theta_HxxHN = M_PI / 2.0 ;        //angle betwwen the x axis of CSA and NH bond
const double theta_HyyHN = 99.0 * M_PI / 180.0;         //angle betwwen the y axis of CSA and NH bond
const double theta_NNH = 18.0 * M_PI / 180.0;



///////////////////////////////// SPEC DENSITY ///////////////////////////////// #Discuss with Fabien
double AutoDD(double w, double *Xarr, double *tau_list)
{	
    double tauc = tau_list[0];
    double Sf2 = Xarr[0];
    double Ss2 = Xarr[1];
    double tau_f = tauc*Xarr[2] / (tauc + Xarr[2]);
    double tau_s = tauc*Xarr[3] / (tauc + Xarr[3]);
    
    double Spec = 2./5. * Sf2 * Ss2 * tauc / (1. + pow(w * tauc, 2)) +\
                2./5. * (1. - Sf2) * tau_f / (1. + pow(w * tau_f, 2)) +\
                2./5. * Sf2 * (1. - Ss2) * tau_s / (1. + pow(w * tau_s, 2));
                        
    return Spec;
}

double AutoDD_HH(double w, double *Xarr, double *tau_list)
{	
    return AutoDD(w, Xarr, tau_list);
}

double AutoCSA(double w, double *Xarr, double *tau_list)
{	
    return AutoDD(w, Xarr, tau_list);
}

double CrossCSADD(double w, double *Xarr, double *tau_list, double angle)
{	
    return AutoDD(w, Xarr, tau_list);
}


///////////////////////////////// RATES /////////////////////////////////
double R1calculation(double B, double *Xarr, double *TaucArr, double *OtherInputs)
{
    double w0_h = GammaH * B;
    double w0_n = GammaN * B;
    
    double J_wn = AutoDD(w0_n, Xarr, TaucArr);
    double J_m = AutoDD(w0_h - w0_n, Xarr, TaucArr);
    double J_p = AutoDD(w0_h + w0_n, Xarr, TaucArr);
    
    double RateCalc = pow(delta_NH, 2)/4. * (3.*J_wn + 6.*J_p + J_m) +\
                    pow(DELTA_N*w0_n, 2)/3.0 * J_wn;

 	return RateCalc;
}

double R2calculation(double B, double *Xarr, double *TaucArr, double *OtherInputs )
{
    double w0_h = GammaH * B;
    double w0_n = GammaN * B;
    
    double J_0 = AutoDD(0.0, Xarr, TaucArr);
    double J_wn = AutoDD(w0_n, Xarr, TaucArr);
    double J_wh = AutoDD(w0_h, Xarr, TaucArr);
    double J_m = AutoDD(w0_h - w0_n, Xarr, TaucArr);
    double J_p = AutoDD(w0_h + w0_n, Xarr, TaucArr);
    

    double RateCalc = pow(delta_NH, 2)/8. * (4.*J_0 + J_m + 3.*J_wn + 6.*J_wh + 6.*J_p) +\
                        pow(DELTA_N*w0_n, 2)/18. * (4.*J_0 + 3.*J_wn);
            
 	return RateCalc;
}

double Sigmacalculation(double B, double *Xarr, double *TaucArr, double *OtherInputs )
{
    double w0_h = GammaH * B;
    double w0_n = GammaN * B;
    
    double J_m = AutoDD(w0_h - w0_n, Xarr, TaucArr);
    double J_p = AutoDD(w0_h + w0_n, Xarr, TaucArr);
 	
 	double RateCalc = pow(delta_NH, 2)/4. * (-J_m + 6.*J_p);

 	return RateCalc;
}

///////////////////////////////// RELAXATION MATRIX /////////////////////////////////
struct mat2d RelaxMatrix(double B, double *Xarr, double *TaucArr, double *OtherInputs) {
    struct mat2d Mat;
    
    double Rex = OtherInputs[0];
    
    double w0_h = GammaH * B;
    double w0_n = GammaN * B;
    
    double J_wn = AutoDD(w0_n, Xarr, TaucArr);
    double J_wh = AutoDD(w0_h, Xarr, TaucArr);
    double J_m = AutoDD(w0_h - w0_n, Xarr, TaucArr);
    double J_p = AutoDD(w0_h + w0_n, Xarr, TaucArr);
    
    double J_HH_0 = AutoDD_HH(0.0, Xarr, TaucArr);
    double J_HH_wh = AutoDD_HH(w0_h, Xarr, TaucArr);
    double J_HH_2wh = AutoDD_HH(2.0*w0_h, Xarr, TaucArr);
    
    double J_CSADD_wh_xx = CrossCSADD(w0_h, Xarr, TaucArr, theta_HxxHN);
    double J_CSADD_wh_yy = CrossCSADD(w0_h, Xarr, TaucArr, theta_HyyHN);
    double J_CSADD_wn = CrossCSADD(w0_n, Xarr, TaucArr, theta_NNH);
    
    ////
    double R1N = pow(delta_NH, 2)/4. * (3.*J_wn + 6.*J_p + J_m) +\
                    pow(DELTA_N*w0_n, 2)/3. * J_wn;
                    
    double R1H = Rex +\
                  pow(delta_NH, 2)/4. * (3.*J_wh + 6.*J_p + J_m) +\
                  1./2.*pow(w0_h, 2)/3. * (pow((sigma_Hxx-sigma_Hzz)/2.0, 2) + pow((sigma_Hyy - sigma_Hzz)/2.0, 2)) * J_wh +\
                pow(delta_HH, 2)/2.0 * (J_HH_0 + 3.0 * J_HH_wh + 6.0 * J_HH_2wh) +\
                pow(DELTA_H*w0_h, 2)/3. * J_HH_0;
                
    double R1NH = Rex +\
                  3.0/4.0 * pow(delta_NH, 2) * (J_wn + J_wh) +\
                  1./2.*pow(w0_h, 2)/3. * (pow((sigma_Hxx-sigma_Hzz)/2.0, 2) + pow((sigma_Hyy - sigma_Hzz)/2.0, 2)) * J_wh +\
                  pow(DELTA_N*w0_n, 2)/3. * J_wn +\
                pow(delta_HH, 2)/2.0 * (J_HH_0 + 3.0 * J_HH_wh + 6.0 * J_HH_2wh) +\
                pow(DELTA_H*w0_h, 2)/3. * J_HH_0;
                
    double R1H_p = R1H;
                    
    double R1NH_p = R1NH;
                  
    ////
    double sigmaNH = pow(delta_NH, 2)/4. * (-J_m + 6.*J_p);
    
    double etaZ_N = delta_NH * DELTA_N*w0_n * J_CSADD_wn;
    
    double etaZ_H = 1./2. * delta_NH * w0_h * ((sigma_Hxx-sigma_Hzz)/2.0 * J_CSADD_wh_xx +\
                                              (sigma_Hyy - sigma_Hzz)/2.0 * J_CSADD_wh_yy);
    
    double sigma_ext = pow(delta_HH, 2)/4. * (-J_HH_0 + 6.*J_HH_2wh);

    
    
 	Mat.m[0][0] = R1N;
 	Mat.m[0][1] = sigmaNH;
 	Mat.m[0][2] = etaZ_N;
 	Mat.m[0][3] = 0.0;
 	Mat.m[0][4] = 0.0;
 	Mat.m[0][5] = 0.0;
 	Mat.m[0][6] = 0.0;
 	
 	Mat.m[1][0] = sigmaNH;
 	Mat.m[1][1] = R1H;
 	Mat.m[1][2] = etaZ_H;
 	Mat.m[1][3] = sigma_ext;
 	Mat.m[1][4] = 0.0;
 	Mat.m[1][5] = sigma_ext;
 	Mat.m[1][6] = 0.0;
 	
 	Mat.m[2][0] = etaZ_N;
 	Mat.m[2][1] = etaZ_H;
 	Mat.m[2][2] = R1NH;
 	Mat.m[2][3] = 0.0;
 	Mat.m[2][4] = sigma_ext;
 	Mat.m[2][5] = 0.0;
 	Mat.m[2][6] = sigma_ext;
 	
 	Mat.m[3][0] = 0.0;
 	Mat.m[3][1] = sigma_ext;
 	Mat.m[3][2] = 0.0;
 	Mat.m[3][3] = R1H_p;
 	Mat.m[3][4] = 0.0;
 	Mat.m[3][5] = 0.0;
 	Mat.m[3][6] = 0.0;
 	
 	Mat.m[4][0] = 0.0;
 	Mat.m[4][1] = 0.0;
 	Mat.m[4][2] = sigma_ext;
 	Mat.m[4][3] = 0.0;
 	Mat.m[4][4] = R1NH_p;
 	Mat.m[4][5] = 0.0;
 	Mat.m[4][6] = 0.0;
 	
 	Mat.m[5][0] = 0.0;
 	Mat.m[5][1] = sigma_ext;
 	Mat.m[5][2] = 0.0;
 	Mat.m[5][3] = 0.0;
 	Mat.m[5][4] = 0.0;
 	Mat.m[5][5] = R1H_p;
 	Mat.m[5][6] = 0.0;
 	
 	Mat.m[6][0] = 0.0;
 	Mat.m[6][1] = 0.0;
 	Mat.m[6][2] = sigma_ext;
 	Mat.m[6][3] = 0.0;
 	Mat.m[6][4] = 0.0;
 	Mat.m[6][5] = 0.0;
 	Mat.m[6][6] = R1NH_p;
 	 	
 	return Mat;
}




