/*
**	File:	baumwelch.c
**	Purpose: Baum-Welch algorithm for estimating the parameters
**		of a HMM model, given an observation sequence. 
*/

#include <stdio.h> 
#include "nrutil.h"
#include "hmmarc.h"
#include <math.h>

#define DELTA 0.000000001
#define MEAN_D 0.002377596

void BaumWelch(HMM *phmm, seqdat *pdat[], int seqnum, double ***alpha, double ***beta, double ***gamma, int *pniter, double *plogprobinit, double *plogprobfinal, int getprob)
{
	int i, j;
	int s;
	int k, m;
	int64_t t, maxT = 0;
	int l = 0;
	//int sizej;

	double	logprobf, logprobb, logprob = 0.0, denominator_e, sum; //, threshold;
	double	**numerator, **denominator, **numerator_e, *numerator_pi, *pil;

	double ***xi, *scale, **Q, **el;
	double delta, /* deltaprev, */logprobprev;
	double /*Aldelta = 0.0,*/ pidelta = 0.0;

	//deltaprev = 10e-70;

	xi = d3darray(1, phmm->N, 1, phmm->N, 0, seqnum - 1);
	for (i = 0; i < seqnum; i++) {
		if (pdat[i]->T > maxT)
			maxT = pdat[i]->T;
	}
	scale = dvector(1, maxT);
	Q = dmatrix(1, phmm->N, 1, phmm->N);
	el = dmatrix(1, phmm->N, 1, phmm->n);
	numerator_e = dmatrix(1, phmm->N, 1, phmm->n);
	numerator = dmatrix(1, phmm->N, 1, phmm->N);
	numerator_pi = dvector(1, phmm->N);
	denominator = dmatrix(1, phmm->N, 1, phmm->N);
	pil = dvector(1, phmm->N);

	for (i = 1; i <= phmm->N; i++)
		pil[i] = phmm->pi[i];

	for (s = 0; s < seqnum; s++){
		Forward(phmm, pdat[s], alpha[s], scale, &logprobf);
		*plogprobinit = logprobf; /* log P(dep |intial model) */
		Backward(phmm, pdat[s], beta[s], scale, &logprobb);
		ComputeGamma(phmm, pdat[s]->T, alpha[s], beta[s], gamma[s]);
		logprob += logprobf;
	}

	if (getprob == 1)
		return;

	delta = DELTA + 1;
	logprobprev = logprob;

	do {	
		fprintf(stderr, "\n====  %d\t%f\t%f  ====\n", l, logprob, delta);

		if (l > 0){
			fprintf(stderr, "pi:\n");
			for (i = 1; i <=phmm->N; i++)
				fprintf(stderr, "%g ", phmm->pi[i]);
			fprintf(stderr, "\n");
			/* fprintf(stderr, "A:\n");
			for (i = 1; i <= phmm->N; i++) {
				for (j = 1; j <= phmm->N; j++) {

					fprintf(stderr, "%g ", Al[i][j]);
				}
				fprintf(stderr, "\n");
			} */
			fprintf(stderr, "T: %.3f\n", phmm->T);
			fprintf(stderr, "emit:\n");
				for (i = 1; i <= phmm->N; i++) {
					for (j = 1; j <= phmm->n; j++) {
						fprintf(stderr, "%g ", phmm->emit[i][j]);
					}
					fprintf(stderr, "\n");
				}
		}

		/* reestimate frequency of state i in time t=1 */
		/* for (i = 1; i <= phmm->N; i++)
			phmm->pi[i] = .001 + .999*gamma[1][i]; */ // pi remains the same

		/* reestimate transition matrix and symbol prob in each state */
		for (i = 1; i <= phmm->N; i++){
			for (j = 1; j <= phmm->N; j++) {
				numerator[i][j] = 0.0;
				denominator[i][j] = 0.0;
			}
		}

		for (i = 1; i <= phmm->N; i++) {
			numerator_pi[i] = 0.0;
			for (j = 1; j <= phmm->n; j++)
				numerator_e[i][j] = 0.0;
		}
		for (s = 0; s < seqnum; s++) {
			for (t = 1; t <= pdat[s]->T-1; t++) {
				ComputeXi(phmm, pdat[s], t, alpha[s], beta[s], xi[s]);
			/* for (i = 1; i <= phmm->L; i++){
				for (j = 1; j <= phmm->L; j++) { 
					for (k = phmm->part[i]; k < phmm->part[i + 1]; k++) {
						for (m = phmm->part[j]; m < phmm->part[j + 1]; m++) {
							numerator[i][j] += xi[k][m];
						}
						denominator[i][j] += gamma[t][k];
					}
				} */
				for (k = 1; k <= phmm->N; k++) {
					for (m = 1; m <= phmm->N; m++) {
						numerator[k][m] += xi[s][k][m];
						denominator[k][m] += gamma[s][t][k];
					}
					if(pdat[s]->snptype[t] > 0) {
						numerator_pi[k] += gamma[s][t][k];
						numerator_e[k][pdat[s]->snptype[t]] += gamma[s][t][k];
					}
				}
			}
			/* fprintf(stderr, "numerator[1][1]: %g denominator[1][1]: %g numerator_pi[1]: %g numerator_pi[2]: %g "
					"numerator_e[1][3]: %g\n", numerator[1][1], denominator[1][1], numerator_pi[1], numerator_pi[2],
					numerator_e[1][3]); */
		}

		/* estimate T and pi from Q matrix */
		for (i = 1; i<= phmm->N; i++) {
			for (j = 1; j<= phmm->N; j++)
				Q[i][j] = numerator[i][j] / denominator[i][j];
		}
		fprintf(stderr, "Q estimated: %.8f %.8f %.8f %.8f\n", Q[1][1], Q[1][2], Q[2][1], Q[2][2]);

		phmm->pi[2] = 1/(Q[2][1]/Q[1][2] + 1);
		phmm->pi[1] = 1 - phmm->pi[2];
		//phmm->T = Q[1][2]/(phmm->pi[2] * MEAN_D);

		pidelta = 0.0;
		sum = 0.0;
		for (i = 1; i <= phmm->N; i++)
			sum += numerator_pi[i];
		for (i = 1; i <= phmm->N; i++) {
			pidelta += fabs(numerator_pi[i] / sum - pil[i]);
			pil[i] = numerator_pi[i] / sum;

		}
		fprintf(stderr, "pi estimated from BW = %.8f\n", pil[2]);
		/* update emission matrix */
		for (i = 1; i <= phmm->N; i++) {
			denominator_e = 0.0;
			for (j = 1; j <= phmm->n; j++)
				denominator_e += numerator_e[i][j];
			for (j = 1; j <= phmm->n; j++)
				phmm->emit[i][j] = el[i][j] = numerator_e[i][j] / denominator_e;
		}

		phmm->emit[2][1] = (phmm->emit[2][1]+phmm->emit[2][2])/2;
		phmm->emit[2][2] = phmm->emit[2][1];
		/*Aldelta /= (phmm->L * phmm->L);*/
		/* Al = numerator / denominator;
		for (k = 1; k <= phmm->N; k++)
			for (m = 1; m <= phmm->N; m++)
				phmm->A[k][m] =  Al / phmm->N; */
		logprob = 0.0;
		for (s = 0; s < seqnum; s++){
			Forward(phmm, pdat[s], alpha[s], scale, &logprobf);
			Backward(phmm, pdat[s], beta[s], scale, &logprobb);
			ComputeGamma(phmm, pdat[s]->T, alpha[s], beta[s], gamma[s]);
			logprob += logprobf;
		}

		/* compute difference between log probability of two iterations */
		delta = logprob - logprobprev;
		logprobprev = logprob;
		
		l++;
	}
	while (delta > DELTA);
	/*while (delta > DELTA); if log probability does not change much, exit */
 
	*pniter = l;
	*plogprobfinal = logprobf; /* log P(dep|estimated model) */
	free_dvector(scale, 1, maxT);
	free_dmatrix(Q, 1, phmm->N, 1, phmm->N);
	free_dmatrix(numerator_e, 1, phmm->N, 1, phmm->n);
	free_dmatrix(el, 1, phmm->N, 1, phmm->n);
	free_dmatrix(numerator, 1, phmm->N, 1, phmm->N);
	free_dmatrix(denominator, 1, phmm->N, 1, phmm->N);
	free_d3darray(xi, 1, phmm->N, 1, phmm->N, 0, seqnum - 1);
}

/* gamma[t][i]: P(X_t = i | Y, theta) */
void ComputeGamma(HMM *phmm, int64_t T, double **alpha, double **beta, double **gamma)
{

	int i, j;
	int64_t t;
	double denominator;

	for (t = 1; t <= T; t++) {
		denominator = 0.0;
		for (j = 1; j <= phmm->N; j++) {
			gamma[t][j] = alpha[t][j]*beta[t][j];
			denominator += gamma[t][j];
		}

		for (i = 1; i <= phmm->N; i++) 
			gamma[t][i] = gamma[t][i]/denominator;
	}
}

/* Xi: P(X_t = i, X_t+1 = j | Y, theta) ; here only at one site (t) */
void ComputeXi(HMM* phmm, seqdat *pdat, int64_t t, double **alpha, double **beta, double **xi)
{
	int i, j;
	double sum;
	double **A;

	A = dmatrix(1, phmm->N, 1, phmm->N);
	computeA(phmm, pdat->gendist[t+1], A, 0);
	sum = 0.0;	
	for (i = 1; i <= phmm->N; i++) 
		for (j = 1; j <= phmm->N; j++) {
			if(pdat->snptype[t+1] > 0)
				xi[i][j] = alpha[t][i] * beta[t+1][j] * A[i][j] * (phmm->emit[j][pdat->snptype[t+1]]);
			else
				xi[i][j] = alpha[t][i] * beta[t+1][j] * A[i][j];
			sum += xi[i][j];
		}

	for (i = 1; i <= phmm->N; i++) 
		for (j = 1; j <= phmm->N; j++) 
			xi[i][j]  /= sum;

	free_dmatrix(A, 1, phmm->N, 1, phmm->N);
}
