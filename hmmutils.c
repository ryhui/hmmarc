/*
**	File:	hmmutils.c
**	Purpose: utilities for reading, writing HMM stuff. 
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "nrutil.h"
#include "hmmarc.h"

//#define REC_RATE 1.2e-8

void ReadHMM(FILE *fp, HMM *phmm)
{
	int i, j;

	if(fscanf(fp, "N: %d\n", &(phmm->N))){};

	if(fscanf(fp, "pi:\n")){};
	phmm->pi = (double *) dvector(1, phmm->N);
	for (i = 1; i <= phmm->N; i++) {
		if(fscanf(fp, "%lf", &(phmm->pi[i]))){};
	}
	if(fscanf(fp, "\n")){};
	if(fscanf(fp, "T: %lf\n", &(phmm->T))){};
	/* if(fscanf(fp, "A:\n")){};
	phmm->A = (double **) dmatrix(1, phmm->N, 1, phmm->N);
	for (i = 1; i <= phmm->N; i++) { 
		for (j = 1; j <= phmm->N; j++) {
			if(fscanf(fp, "%lf", &(phmm->A[i][j]))){};
		}
		if(fscanf(fp, "\n")){};
	} */

	if(fscanf(fp, "n: %d\n", &(phmm->n))){};

	if(fscanf(fp, "emit:\n")){};
	phmm->emit = (double **) dmatrix(1, phmm->N, 1, phmm->n);
	for (i = 1; i <= phmm->N; i++) {
		for (j = 1; j <= phmm->n; j++) {
			if(fscanf(fp, "%lf", &(phmm->emit[i][j]))){};
		}
		if(fscanf(fp,"\n")){};
	}
}

void FreeHMM(HMM *phmm)
{
	//free_dmatrix(phmm->A, 1, phmm->N, 1, phmm->N);
	free_dmatrix(phmm->emit, 1, phmm->N, 1, phmm->n);
	free_dvector(phmm->pi, 1, phmm->N);
}

void PrintHMM(FILE *fp, HMM *phmm)
{
	int i, j;

	fprintf(fp, "N: %d\n", phmm->N); 
 
	fprintf(fp, "pi:\n");
	for (i = 1; i <= phmm->N; i++) {
		fprintf(fp, "%g ", phmm->pi[i]);
	}
	fprintf(fp, "\n");

	fprintf(fp, "T: %g\n", phmm->T);
	/*fprintf(fp, "A:\n");
		for (i = 1; i <= phmm->N; i++) {
			for (j = 1; j <= phmm->N; j++) {
				fprintf(fp, "%g ", phmm->A[i][j] );
			}
			fprintf(fp, "\n");
		}*/

	fprintf(fp, "n: %d\n", phmm->n);

	fprintf(fp, "emit:\n");
		for (i = 1; i <= phmm->N; i++) {
			for (j = 1; j <= phmm->n; j++) {
				fprintf(fp, "%g ", phmm->emit[i][j] );
			}
			fprintf(fp, "\n");
	}

}

void computeA(HMM *phmm, double dst, double **A, int logtrans)
{
	A[1][2] = phmm->pi[2] * (1 - exp(- dst * 0.01 * phmm->T));
	A[1][1] = 1 - A[1][2];
	A[2][1] = (1 - exp(- dst * 0.01 * phmm->T)) * phmm->pi[1];
	A[2][2] = 1 - A[2][1];
	if(logtrans) {
		int i, j;
		for (i = 1; i <= phmm->N; i++) {
			for (j = 1; j<= phmm->N; j++)
				A[i][j] = log(A[i][j]);
		}
	}
}

