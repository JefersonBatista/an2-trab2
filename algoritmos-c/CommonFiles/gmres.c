#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "gmres.h"

void subtrair(double *uc, double beta, double *u,int tamU) {
	int i = 0;
	for(i = 0; i < tamU; i++) {
		uc[i] -= beta * u[i];
	}
}

double norm(double *v, int tamV) {
	double soma = 0;
	int i;
	for(i = 0; i < tamV; i++) {
		soma += pow(v[i],2);
	}
	return sqrt(soma);
}

double produto_interno(double *v1, double *v2, int k) {
	int i = 0;
	double prod = 0;
	for(i = 0; i < k; i++) {
		prod += v1[i] * v2[i];
	}
	return prod;
}

void MATRIX_matvec(MAT *A, double *x, double *b) {
    double *AA = A->AA;
    int *IA = A->IA;
    int *JA = A->JA;
    
    int i, j, n = A->n;
    for(i = 0; i < n; i++) {
        b[i] = 0.0;
        for(j = IA[i]; j < IA[i+1]; j++) {
            b[i] += AA[j] * x[JA[j]];
        }
    }
}

	/*--------------------------------------------------------
	 * Método GMRES com restart
	 * Entrada: Matriz A, Vetor b, Solucao inicial x0
	 * Numero de elementos na base Kmax,
	 * Numero máximo de reinicializacao lmax, Tolerancia tol
	 *--------------------------------------------------------*/

double GMRES(MAT *A, double *x, double *b, double tol, int kmax, int lmax) {

	int n = A->n;
	double epson, *p, rol;
	int i, j, k;
	int iter = 1;

	/* Inicializar ebar, u, H, c, s e y */
	
	double *ebar = (double*) calloc((kmax + 1), sizeof(double));
	double *c = (double*) calloc(kmax, sizeof(double));
	double *s = (double*) calloc(kmax, sizeof(double));
	double *y = (double*) calloc(kmax, sizeof(double));
	double **u = (double**) calloc((kmax + 1), sizeof(double*));
	double **H = (double**) calloc((kmax + 1), sizeof(double*));

	for(k = 0; k < kmax + 1; k++) {
		u[k] = (double*) calloc(n, sizeof(double));
		H[k] = (double*) calloc(kmax, sizeof(double));
	}

	epson = tol * norm(b, n);

	do {
		i = 0;
		
		MATRIX_matvec(A, x, p);
		
		for(k = 0; k < n, k++)
			u[i][k] = b[k] - p[k];
	
		free(p);
		ebar[i] = norm(u[i], n);

		for(k = 0; k < n; k++)
			u[i][k] = u[i][k] / ebar[i];
		
		rol = ebar[i];
		
		while(rol > epson && i < kmax) {
		
			MATRIX_matvec(A, u[i], p);
			
			for(k = 0; k < n; k++)
				u[i+1][k] = p[k];
		
			free(p);
			
			/*---------------------
			 * Gram-Schmidt
			 *---------------------*/

			for(j = 0; j <= i; j++) {
				H[j][i] = produto_interno(u[i+1], u[j], n);
				subtrair(u[i+1], H[j][i], u[j], n);
			}
			
			H[i+1][i] = norm(u[i+1], n);
			
			for(k = 0; k < n; k++)
				u[i+1][k] = u[i+1][k] / H[i+1][i];
			
			/*----------------------------
			 * Aplica as rotações de Givens
			 * Algoritmo QR
			 *----------------------------*/

			for(j = 0; j <= i-1; j++) {
				double hji = c[j] * H[j][i] + s[j] * H[j+1][i];
				double hj1i = -s[j] * H[j][i] + c[j] * H[j+1][i];
				H[j][i] = hji;
				H[j+1][i] = hj1i;
			}

			r = H[i][i] * H[i][i] + H[i+1][i] * H[i+1][i];
			r = sqrt(r);
			
			c[i] = H[i][i] / r;
			s[i] = H[i+1][i] / r;
			H[i][i] = r;
			H[i+1][i] = 0.0;

			ebar[i+1] = -s[i] * ebar[i];
			ebar[i] = c[i] * ebar[i];

			rol = fabs(ebar[i+1]);

			i++;
		}
			/*----------------------------
			 * Fim da iteração do GMRES
			 *----------------------------*/

		i--;

			/*---------------------------------
			 * Resolve o Sistema Liner Hy = e
			 *---------------------------------*/

		for(j = i; j >= 0; j--) {
			int l;
			double soma = 0;
			for(l = j+1; l <= i; l++)
				soma += H[j][l] * y[l];
			y[j] = (ebar[j]-soma) / H[j][j];
		}

		for(j = 0; j <= i; j++) {
			for(k = 0; k < n; k++)
				x[k] += y[j] * u[j][k];
		}

		iter++;
	} while((rol >= epson) && (iter < lmax));

	free(ebar);
	free(c);
	free(s);
	free(y);

	for(k = 0; k < kmax + 1; k++) {
		free(u[k]);
		free(H[k]);
	}

	free(u);
	free(H);

	printf("(GMRES) Ciclo %d: iteracao = %d\n", l, iter);

	return x;
}



