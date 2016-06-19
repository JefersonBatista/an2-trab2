#include "heads.h"
#include "protos.h"

void subtrair(double *uc, double beta, double *u, int tamU);
double norm(double *v, int tamV);
double produto_interno(double *v1, double *v2, int k);
void matvec_pc(MAT *L, MAT *U, double *x, double *b);
void GMRES(MAT *A, double *b, double *x, double tol, int kmax, int lmax);
void GMRES_pc(MAT *A, MAT *L, MAT *U, double *b, double *x, double tol, int kmax, int lmax);

