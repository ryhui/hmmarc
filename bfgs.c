#include <stdio.h>
#include <math.h>
#include "hmmarc.h"
#include "nrutil.h"
#include "lbfgsb.h"

#define MIN_PROB 0.000001
#define MAX_PROB 0.999999

void HMMtoTheta(HMM *phmm, double *theta)
// Put hmm parameters into a vector theta
{
	theta[0] = phmm->T;
	theta[1] = phmm->pi[2];
	theta[2] = phmm->emit[1][1] / phmm->emit[1][3];
	theta[3] = phmm->emit[1][2] / phmm->emit[1][3];
	theta[4] = phmm->emit[2][1] / phmm->emit[2][3];
	theta[5] = phmm->emit[2][1] / phmm->emit[2][3];
	//theta[6] = phmm->emit[2][2];
	//theta[7] = phmm->emit[2][3];
}

void ThetatoHMM(HMM *phmm, double *theta)
{
	phmm->T = theta[0];
	phmm->pi[2] = theta[1];
	phmm->emit[1][1] = theta[2] / (theta[2] + theta[3] + 1);
	phmm->emit[1][2] = theta[3] / (theta[2] + theta[3] + 1);
	phmm->emit[1][3] = 1 / (theta[2] + theta[3] + 1);
	phmm->emit[2][1] = theta[4] / (theta[4] + theta[5] + 1);
	phmm->emit[2][2] = theta[5] / (theta[4] + theta[5] + 1);
	phmm->emit[2][3] = 1 / (theta[4] + theta[5] + 1);
	phmm->pi[1] = 1 - theta[1];
	//phmm->emit[1][3] = 1 - theta[2] - theta[3];
	//phmm->emit[2][3] = 1 - theta[4] - theta[5];
}

void drvTrans(double ***drv, HMM *phmm, double d, int parnum)
{
	int i, j, p;
	for (p = 0; p < parnum; p++) {
		for (i = 1; i <= phmm->N; i++) {
			for (j = 1; j <= phmm->N; j++)
				drv[p][i][j] = 0.0;
		}
	}
	/* T */
	drv[0][1][1] = -phmm->pi[2]*d;
	drv[0][1][2] = phmm->pi[2]*d;
	drv[0][2][1] = phmm->pi[1]*d;
	drv[0][2][2] = -phmm->pi[1]*d;
	/* alpha (phmm->pi[2]) */
	drv[1][1][1] = -phmm->T*d;
	drv[1][1][2] = phmm->T*d;
	drv[1][2][1] = -phmm->T*d;
	drv[1][2][2] = phmm->T*d;
}
void drvEmit(double **drv, HMM *phmm, int state, int parnum)
{
	int i, p;
	for (p = 0; p < parnum; p++) {
		for (i = 1; i < phmm->N; i++) {
			drv[p][i] = 0.0;
		}
	}
	switch (state) {
	case 1:
		drv[2][1] = 1;
		drv[4][2] = 1;
		break;
	case 2:
		drv[3][1] = 1;
		drv[5][2] = 1;
		break;
	case 3:
		drv[2][1] = -1;
		drv[3][1] = -1;
		drv[4][2] = -1;
		drv[5][2] = -1;
		break;
	}
}

void drvPi(double **drv, HMM *phmm, int parnum)
{
	int i, p;
	for (p = 0; p < parnum; p++) {
		for (i = 1; i <= phmm->N; i++)
			drv[p][i] = 0.0;
	}
	drv[1][1] = -1;
	drv[1][2] = 1;
}

double getFG(double *theta, double *gradient, HMM *phmm, seqdat *pdat[], int seqnum, int parnum)
{
	seqdat* pseqdat;
	int i, s, j, p;
	int64_t t;
	double logp = 0.0;
	double *lambda, *lambda_new, **psi, **psi_new, scale, scale_new;
	double ***tdrv, **edrv, **pdrv;
	tdrv = d3darray(1, phmm->N, 1, phmm->N, 0, parnum);
	edrv = dmatrix(0, parnum, 1, phmm->N);
	pdrv = dmatrix(0, parnum, 1, phmm->N);

	double **A;
	A = dmatrix(1, phmm->N, 1, phmm->N);
	lambda = dvector(1, phmm->N);
	lambda_new = dvector(1, phmm->N);
	psi = dmatrix(0, parnum, 1, phmm->N);
	psi_new = dmatrix(0, parnum, 1, phmm->N);
	ThetatoHMM(phmm, theta);

	for (i = 0; i < parnum; i++)
		gradient[i] = 0.0;
	/* Calculate log likelihood and gradient at the same time using the revised
	 * forward algorithm (Eq. 3.3, 3.4, 4.1 & 4.2)
	 * */

	for (s = 0; s < seqnum; s++){
		pseqdat = pdat[s];
		scale = 0.0;
		scale_new = 0.0;
		for (i = 1; i <= phmm->N; i++) {
			for (p = 0; p < parnum; p++) {
				psi[p][i] = 0.0;
				psi_new[p][i] = 0.0;
			}
			lambda[i] = phmm->pi[i] * phmm->emit[i][pseqdat->snptype[1]];
			scale += lambda[i];
			/* Initialise gradients */
			drvEmit(edrv, phmm, pseqdat->snptype[1], parnum);
			drvPi(pdrv, phmm, parnum);
			for (p = 0; p < parnum; p++)
				psi[p][i] = edrv[p][i] * phmm->pi[i] + phmm->emit[i][pseqdat->snptype[1]]*pdrv[p][i];
		}
		logp += log(scale);
		// fprintf(stderr, "logp = %g at t = %d\n", logp, 1);
		for (t = 2; t <= pseqdat->T; t++) {
			computeA(phmm, pseqdat->gendist[t], A, 0);
			drvEmit(edrv, phmm, pseqdat->snptype[t], parnum);
			drvTrans(tdrv, phmm, pseqdat->gendist[t]*0.01, parnum);
			for (i = 1; i <= phmm->N; i++){
				lambda_new[i] = 0.0;
				for (j = 1; j <= phmm->N; j++){
					lambda_new[i] += lambda[j]*A[j][i]*phmm->emit[i][pseqdat->snptype[t]] / scale;
					for (p = 0; p < parnum; p++) {
						psi_new[p][i] += psi[p][j]*phmm->emit[i][pseqdat->snptype[t]]*A[j][i] +
								lambda[j]*edrv[p][i]*A[j][i] +
								lambda[j]*phmm->emit[i][pseqdat->snptype[t]]*tdrv[p][j][i];
					}
				}
				scale_new += lambda_new[i];
			}
			for (i = 1; i <= phmm->N; i++) {
				lambda[i] = lambda_new[i];
				lambda_new[i] = 0;
				for (p = 0; p < parnum; p++) {
					psi[p][i] = psi_new[p][i] / scale;
					psi_new[p][i] = 0.0;
				}
			}
			if (scale_new <= 0) {
				fprintf(stderr, "Error: scale below 0 at t = %ld\n", t);
			}
			/* if (isnan(scale_new))
				fprintf(stderr, "Error: scale is nan at t = %ld\n", t); */
			scale = scale_new;
			scale_new = 0.0;
			logp += log(scale);
			//fprintf(stderr, "logp = %g at t = %d\n", logp, t);
		}

		for (p = 0; p < parnum; p++) {
			for (i = 1; i <= phmm->N; i++) {
				gradient[p] += psi[p][i]/scale;
			}
		}

		/* Add a term -M * (theta[2]+theta[3]+theta[4] - 1)^2 -M * (theta[5]+theta[6]+theta[7] - 1)^2 */
		for (p = 2; p < 4; p++)
			gradient[p] = gradient[p] * (theta[2]+theta[3]+1 - theta[p]) / pow(theta[2]+theta[3]+1, 2);
		for (p = 4; p < 6; p++)
			gradient[p] = gradient[p] * (theta[4]+theta[5]+1 - theta[p]) / pow(theta[4]+theta[5]+1, 2);

		/* take care of two Lagrange terms
		for (p = 2; p <= 4; p++)
			gradient[p] -= -theta[8];
		for (p = 5; p <= 7; p++)
			gradient[p] -= -theta[9];
		gradient[8] -= theta[2] + theta[3] + theta[4] - 1;
		gradient[9] -= theta[5] + theta[6] + theta[7] - 1; */
	}
	//fprintf(stderr, "log likelihood: %g\n", logp);
	free_dvector(lambda, 1, phmm->N);
	free_dvector(lambda_new, 1, phmm->N);
	free_dmatrix(psi, 0, parnum, 1, phmm->N);
	free_dmatrix(psi_new, 0, parnum, 1, phmm->N);

	free_dmatrix(edrv, 0, parnum, 1, phmm->N);
	free_d3darray(tdrv, 1, phmm->N, 1, phmm->N, 0, parnum);
	free_dmatrix(pdrv, 0, parnum, 1, phmm->N);
	//logp -= M * pow(theta[2]+theta[3]+theta[4] - 1, 2) + M * pow(theta[5]+theta[6]+theta[7] - 1, 2);
	fprintf(stderr, "logp = %g\n", logp);

	for (p = 0; p < parnum; p++)
		gradient[p] = -gradient[p];
	return -logp;
}

void trainbfgs(HMM *phmm, seqdat *pdat[], int seqnum)
{
	//int i, j, p, s;
	int p;
	//int64_t t;
	integer iprint = 101;
	double factr = 1e7;
	double pgtol = 1e-9;
	logical lsave[4];
	double dsave[29], wa[100000];
	integer parnum = 6, m = 50, *nbd;
	integer iwa[30];
	integer taskValue, csaveValue, *isave;
	integer *task=&taskValue;
	integer *csave=&csaveValue;
	isave = calloc(44, sizeof(integer));

	double *theta, *gradient, logp;//*lambda, *lambda_new, **psi, **psi_new, scale, scale_new;
	double /* ***tdrv, **edrv, **pdrv,*/ *lbound, *ubound;
	/* tdrv = d3darray(1, phmm->N, 1, phmm->N, 0, parnum);
	edrv = dmatrix(0, parnum, 1, phmm->N);
	pdrv = dmatrix(0, parnum, 1, phmm->N); */
	theta = dvector(0, parnum);
	gradient = dvector(0, parnum);
	lbound = dvector(0, parnum);
	ubound = dvector(0, parnum);
	nbd = calloc(parnum, sizeof(integer));

	HMMtoTheta(phmm, theta);
	//theta[8] = 1;
	//theta[9] = 1;

	for (p = 0; p <= 5; p++)
		lbound[p] = MIN_PROB; // lower bound, except for lagrange multipliers

	for (p = 2; p <= 5; p++) {
		//ubound[p] = MAX_PROB; // upper bound for admix. prop and emit probs
		nbd[p] = 1;
	}
	nbd[0] = 1; // admix time has only lower bound
	nbd[1] = 2;
	ubound[1] = 1;
	//nbd[8] = 0;
	//nbd[9] = 0;

	*task = (integer)START;

EVALUATE:
	setulb(&parnum, &m, theta, lbound, ubound, nbd, &logp, gradient, &factr, &pgtol, wa, iwa, task, &
			iprint, csave, lsave, isave, dsave);
	if (IS_FG(*task)) {
		logp = getFG(theta, gradient, phmm, pdat, seqnum, parnum);
		goto EVALUATE;
	}
	if (*task == NEW_X)
		goto EVALUATE;
	ThetatoHMM(phmm, theta);
	PrintHMM(stdout, phmm);
	fprintf(stdout, "%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n",
			phmm->T, phmm->pi[2], phmm->emit[1][1], phmm->emit[1][2], phmm->emit[1][3],
			phmm->emit[2][1], phmm->emit[2][2], phmm->emit[2][3], logp);
}

