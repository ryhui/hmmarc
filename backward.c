/*
**	File:	backward.c
**	Purpose: Backward algorithm for computing the probabilty
**		of observing a sequence given a HMM model parameter.
*/

#include <stdio.h>
#include "nrutil.h"
#include "hmmarc.h"

void Backward(HMM *phmm, seqdat *pdat, double **beta, double *scale, double *pprob)
{
	int i, j; /* state indices */
	int64_t t; /* time index */
	double sum;
	double **A;
 
	A = dmatrix(1, phmm->N, 1, phmm->N);
	/* 1. Initialization */
	for (i = 1; i <= phmm->N; i++)
		beta[pdat->T][i] = 1.0/scale[pdat->T]; 
 
	/* 2. Induction */
	for (t = pdat->T - 1; t >= 1; t--) {
		computeA(phmm, pdat->gendist[t + 1], A, 0);
		for (i = 1; i <= phmm->N; i++) {
			sum = 0.0;
			for (j = 1; j <= phmm->N; j++) {
				if(pdat->snptype[t+1] > 0)
					sum += A[i][j] * phmm->emit[j][pdat->snptype[t+1]] * beta[t+1][j];
				else
					sum += A[i][j] * beta[t+1][j];
			}

			beta[t][i] = sum/scale[t];
			if (beta[t][i] == 0) //2.3e-308)
				beta[t][i] = 1e-300;
				//fprintf(stderr, "Warning: normalised beta underflow at position %d.\n", pdat->pos[t]);
			if (beta[t][i] > 1.7e308)
				fprintf(stderr, "Warning: normalised beta overflow at position %d.\n", pdat->pos[t]);
		}
	}
	free_dmatrix(A, 1, phmm->N, 1, phmm->N);
}
