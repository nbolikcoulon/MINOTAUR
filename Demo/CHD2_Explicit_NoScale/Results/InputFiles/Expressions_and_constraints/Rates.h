struct mat2d {

    double m[3][3];

};


struct mat2d RelaxMatrix(double Bo, double *X, double *Tauc, double *OtherInputs);

double R1calculation(double Bo, double *X, double *Tauc, double *OtherInputs);
double R2calculation(double Bo, double *X, double *Tauc, double *OtherInputs);
double Sigmacalculation(double Bo, double *X, double *Tauc, double *OtherInputs);
double EtaXYcalculation(double Bo, double *X, double *Tauc, double *OtherInputs);
double EtaZcalculation(double Bo, double *X, double *Tauc, double *OtherInputs);