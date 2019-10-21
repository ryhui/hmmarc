/*
**	File:	forward.c
**	Purpose: Foward algorithm for computing the probabilty 
**		of observing a sequence given a HMM model parameter.
*/

#include <stdio.h>
#include "nrutil.h"
#include "hmmarc.h"

void Forward(HMM *phmm, seqdat *pdat, double **alpha, double *scale, double *pprob)
/*  pprob is the log probability */
{
	int i, j;	/* state indices */
	int64_t t;	/* time index */

	double sum;	/* partial sum */
	double **A;

	A = dmatrix(1, phmm->N, 1, phmm->N);
	/* 1. Initialization */
	scale[1] = 0.0;
	for (i = 1; i <= phmm->N; i++) {
		if(pdat->snptype[1] > 0) {
			alpha[1][i] = phmm->pi[i]* phmm->emit[i][pdat->snptype[1]];
		} else {
			alpha[1][i] = 1;
		}
		scale[1] += alpha[1][i];
	}
/*printf("%d %f %f %d %e %e \n", pdat->dep[1], phmm->cn[i], pdat->cvg, pdat->ref[1], pemit(phmm, i, pdat, 1), alpha[1][i]);*/
	if (scale[1] == 0){
		fprintf(stderr, "error: scale underflow at position %d.\n", pdat->pos[1]);
	}
	for (i = 1; i <= phmm->N; i++) 
		alpha[1][i] /= scale[1]; 
	
	/* 2. Induction */
	for (t = 1; t <= pdat->T - 1; t++) {
		scale[t+1] = 0.0;
		computeA(phmm, pdat->gendist[t + 1], A, 0);
		for (j = 1; j <= phmm->N; j++) {
			sum = 0.0;
			for (i = 1; i <= phmm->N; i++) 
				sum += alpha[t][i]* A[i][j];
			//if(pdat->snptype[t+1] > 0) {
			alpha[t + 1][j] = sum * phmm->emit[j][pdat->snptype[t + 1]];
			//} else {
			//	alpha[t + 1][j] = sum;
			//}
			scale[t + 1] += alpha[t + 1][j];
		}
		if (scale[t + 1] == 0){
			fprintf(stderr, "error: scale underflow at position %d.\n", pdat->pos[t + 1]);
		}
/*printf("%d %f\n", t+ 1,scale[t + 1]);*/
		for (j = 1; j <= phmm->N; j++) {
			alpha[t + 1][j] /= scale[t + 1];
		}
	}

	/* 3. Termination */
	*pprob = 0.0;

	for (t = 1; t <= pdat->T; t++)
		*pprob += log(scale[t]);

	free_dmatrix(A, 1, phmm->N, 1, phmm->N);
}
