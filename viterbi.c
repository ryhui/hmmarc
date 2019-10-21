/*
**      File:   viterbi.c
**      Purpose: Viterbi algorithm for computing the maximum likelihood
**		state sequence and probablity of observing a sequence given the model. 
*/

#include <math.h>
#include "hmmarc.h"
#include "nrutil.h"

#define VITHUGE  100000000000.0

void Viterbi(HMM *phmm, seqdat *pdat, uint8_t *q, double *pprob)
{
	int i, j;   /* state indices */
	int64_t t;      /* time index */
 
	int     maxvalind;
	double  maxval, val;
	uint8_t **psi;
	double *delta;
	double *deltaprev;
	double **A;

	psi = cmatrix(1, pdat->T, 1, phmm->N);
	delta = dvector(1, phmm->N);
	deltaprev = dvector(1, phmm->N);
	A = dmatrix(1, phmm->N, 1, phmm->N);

	/* 0. Preprocessing */
	/* for (i = 1; i <= phmm->N; i++)
		phmm->pi[i] = log(phmm->pi[i]);
	for (i = 1; i <= phmm->N; i++) 
		for (j = 1; j <= phmm->N; j++) {
			phmm->A[i][j] = log(phmm->A[i][j]);
		}
	*/
	for (i = 1; i <= phmm->N; i++)
		for (j = 1; j <= phmm->n; j++) {
			phmm->emit[i][j] = log(phmm->emit[i][j]);
		}

	/* 1. Initialization  */
	for (i = 1; i <= phmm->N; i++) {
		if (pdat->snptype[1] > 0)
			deltaprev[i] = log(phmm->pi[i]) + phmm->emit[i][pdat->snptype[1]];
		else
			deltaprev[i] = log(phmm->pi[i]);
		psi[1][i] = 0;
	}
 
	/* 2. Recursion */
	for (t = 2; t <= pdat->T; t++) {
		computeA(phmm, pdat->gendist[t], A, 1);
		for (j = 1; j <= phmm->N; j++) {
			maxval = -VITHUGE;
			maxvalind = 1;
			for (i = 1; i <= phmm->N; i++) {
				val = deltaprev[i] + A[i][j];
				if (val > maxval) {
					maxval = val;
					maxvalind = i;
				}
			}
			if (pdat->snptype[t] > 0)
				delta[j] = maxval + phmm->emit[j][pdat->snptype[t]];
			else
				delta[j] = maxval;
			psi[t][j] = (uint8_t)maxvalind;
		}
		// fprintf(stderr, "\n%ld %d\t%f %f\t\%f %f", t, pdat->snptype[t], deltaprev[1], deltaprev[2], delta[1], delta[2]);
		for (j = 1; j <= phmm->N; j++)
			deltaprev[j] = delta[j];
//fprintf(stderr, "\t%d", psi[t][j]);

	}
 
	/* 3. Termination */
	*pprob = -VITHUGE;
	q[pdat->T] = 1;
	for (i = 1; i <= phmm->N; i++) {
		if (delta[i] > *pprob) {
			*pprob = delta[i];
			q[pdat->T] = (uint8_t)i;
		}
	}
 
	/* 4. Path (state sequence) backtracking */
	for (t = pdat->T - 1; t >= 1; t--)
		q[t] = (uint8_t)psi[t+1][q[t+1]];

	free_dvector(delta, 1, phmm->N);
	free_dvector(deltaprev, 1, phmm->N);
	free_cmatrix(psi, 1, pdat->T, 1, phmm->N);
}
