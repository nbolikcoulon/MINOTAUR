struct mat2d {

    double m[7][7];

};


struct mat2d RelaxMatrix(double Bo, double *X, double *Tauc, double *OtherInputs);

double R1calculation(double Bo, double *X, double *Tauc, double *OtherInputs);
double R2calculation(double Bo, double *X, double *Tauc, double *OtherInputs);
double Sigmacalculation(double Bo, double *X, double *Tauc, double *OtherInputs);
