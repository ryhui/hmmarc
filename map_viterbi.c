/*
 * map_viterbi.c
 *
 *  Created on: 22 Aug 2016
 *      Author: ruoyun
 */

#include <math.h>
#include "hmmarc.h"
#include "nrutil.h"

#define VITHUGE  100000000000.0

void map_Viterbi(HMM *phmm, seqdat *pdat, mapdat *mdat, uint8_t *q, double *pprob)
{
	int i, j;   /* state indices */
	int64_t t;      /* time index */

	int     maxvalind;
	int64_t mappos=1;
	double  maxval, val;
	uint8_t **psi;
	double *delta;
	double *deltaprev;

	psi = cmatrix(1, pdat->T, 1, phmm->N);
	delta = dvector(1, phmm->N);
	deltaprev = dvector(1, phmm->N);

	/* 0. Preprocessing */
	/*for (i = 1; i <= phmm->N; i++)
		phmm->pi[i] = log(phmm->pi[i]); */
	while (pdat->startpos >= mdat->loc[mappos + 1])
		mappos++;
	computeA(phmm, mdat, mappos);
	for (i = 1; i <= phmm->N; i++)
	 	for (j = 1; j <= phmm->N; j++) {
	 		phmm->A[i][j] = log(phmm->A[i][j]);
	 	}
	for (i = 1; i <= phmm->N; i++)
		for (j = 1; j <= phmm->n; j++) {
			phmm->emit[i][j] = log(phmm->emit[i][j]);
		}

	/* 1. Initialization  */
	for (i = 1; i <= phmm->N; i++) {
		if (pdat->snptype[1] > 0)
			deltaprev[i] = phmm->pi[i] + phmm->emit[i][pdat->snptype[1]];
		else
			deltaprev[i] = phmm->pi[i];
		psi[1][i] = 0;
	}

	/* 2. Recursion */
	for (t = 2; t <= pdat->T; t++) {
//fprintf(stderr, "\n%d %d %d %f\t", t, pdat->gc[t], pdat->dep[t], pdat->ref[t]);
		if (mappos < mdat->l && t + pdat->startpos - 1 >= mdat->loc[mappos]) {
			computeA(phmm, mdat, t);
			/* for (i = 1; i <= phmm->N; i++) */
			/* 	for (j = 1; j <= phmm->N; j++) { */
			/* 		phmm->A[i][j] = log(phmm->A[i][j]); */
			/* 	} */
			mappos++;
			for (i = 1; i <= phmm->N; i++)
				for (j = 1; j <= phmm->N; j++) {
					phmm->A[i][j] = log(phmm->A[i][j]);
				}
		}
		for (j = 1; j <= phmm->N; j++) {
			maxval = -VITHUGE;
			maxvalind = 1;
			for (i = 1; i <= phmm->N; i++) {
				val = deltaprev[i] + (phmm->A[i][j]);
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
//fprintf(stderr, "%d %d %d %f %f\n", t, pdat->gc[t], pdat->dep[t], pdat->ref[t], lnpemit(phmm, i, pdat, t));
		for (j = 1; j <= phmm->N; j++){
			deltaprev[j] = delta[j];
//fprintf(stderr, "\t%d", psi[t][j]);
		}
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

