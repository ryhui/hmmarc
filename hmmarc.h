/*
**      File:   hmmarc.h
**      Purpose: data structures used for HMMARC
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdint.h>

typedef struct {
	int N;		/* number of hidden states */
	double *pi;	/* ancestral proportion; admixture proportion */
	double T;	/* admixture time (in gens) */
	int n;	/* number of emission types (snp states) */
	double **emit; /* emission probabilities */
} HMM;

typedef struct{
	int64_t T;		/* length of data sequences */
	int *pos; /* position of snps */
	uint8_t *snptype;	/* snp type array */
	double *gendist; /* genetic distance between current and previous SNP */
} seqdat;

void ReadHMM(FILE *fp, HMM *phmm);
void PrintHMM(FILE *fp, HMM *phmm);
void FreeHMM(HMM *phmm);

int64_t countlines(FILE *fdat);
//void readrange(FILE *fp, seqdat *pdat);
void readseqs(FILE *fdat, seqdat *pdat);
//void readsites(FILE *fsites, int64_t *s);
void freeseqdat(seqdat *pdat);

// void readmap(FILE *fdat, mapdat *pmap);
// void freemap(mapdat *pmap);
void computeA(HMM *phmm, double dst, double **A, int logtrans);

void Forward(HMM *phmm, seqdat *pdat, double **alpha, double *scale, double *pprob);
void Backward(HMM *phmm, seqdat *pdat, double **beta, double *scale, double *pprob);
void BaumWelch(HMM *phmm, seqdat *pdat[], int seqnum, double ***alpha, double ***beta, double ***gamma, int *niter, double *plogprobinit, double *plogprobfinal, int getprob);
/* void map_Forward(HMM *phmm, seqdat *pdat, mapdat *mdat, double **alpha, double *scale, double *pprob);
void map_Backward(HMM *phmm, seqdat *pdat, mapdat *mdat, double **beta, double *scale, double *pprob);
void map_BaumWelch(HMM *phmm, seqdat *pdat[], int seqnum, mapdat *mdat, double ***alpha, double ***beta, double ***gamma, int *niter, double *plogprobinit, double *plogprobfinal, int getprob);
*/
void ComputeGamma(HMM *phmm, int64_t T, double **alpha, double **beta, double **gamma);
void ComputeXi(HMM* phmm, seqdat *pdat, int64_t t, double **alpha, double **beta, double **xi);

void Viterbi(HMM *phmm, seqdat *pdat, uint8_t *q, double *pprob);
void trainbfgs(HMM *phmm, seqdat *pdat[], int seqnum);
// void map_Viterbi(HMM *phmm, seqdat *pdat, mapdat *mdat, uint8_t *q, double *pprob);

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define ROUND(x) (int)((x) < 0 ? (x) - 0.5 : (x) + 0.5 )
