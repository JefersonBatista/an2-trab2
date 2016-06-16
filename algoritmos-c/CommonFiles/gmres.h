#include "heads.h"

void subtrair(double *uc, double beta, double *u, int tamU);
double norm(double *v, int tamV);
double produto_interno(double *v1, double *v2, int k);
double* GMRES(MAT *A, double *x0, double *b, double tol, int kmax, int lmax);
