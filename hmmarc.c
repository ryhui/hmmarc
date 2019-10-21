/*
**      File:   hmmarc.c
**      Purpose: arg handler & setup
*/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <unistd.h>
#include "nrutil.h"
#include "hmmarc.h"

void Usage(char *name);

int main(int argc, char **argv)
{
	HMM	hmm;
	// seqdat dat;
	FILE *fp;
	// mapdat mdat;
	//int64_t *sites;
	//int64_t numsites=0;
	//int usemap=0;
	int getprob=0;
	int	vflg=0;
	int	c, i;
	extern char *optarg;
	extern int optind, opterr, optopt;

	while ((c= getopt(argc, argv, "vhm:r:s:")) != EOF)
		switch (c) {
		case 'v': 
			vflg++; 
			break;
		case 'h': 
			Usage(argv[0]);
			exit(1);
			break;
		case 'm':
			/* read hmm model */
			fp = fopen(optarg, "r");
			if (fp == NULL) {
				fprintf(stderr, "Error: File %s not found\n", optarg);
				exit (1);
			}
			ReadHMM(fp, &hmm);
			fclose(fp);
			break;
		}

  	/* read data sequences*/
	int seqnum = argc - optind - 1;

	// fprintf(stderr, "Initialized testseqdat has T: %ld, startpos: %ld\n", testseqdat.T, testseqdat.startpos);
	seqdat* pseqdat[seqnum];
	// This doesn't work: seqdat = malloc(seqnum * sizeof(seqdat));
	for (i = 0; i < seqnum; i++) {
		pseqdat[i] = (seqdat*) malloc(sizeof(seqdat));
		fp = fopen(argv[optind + i + 1], "r");
		if (fp == NULL) {
			fprintf(stderr, "Error: File %s not found\n", argv[optind + 1]);
			exit (1);
		}
		fprintf(stderr, "Reading %s...\n", argv[optind + i + 1]);
		pseqdat[i]->T = countlines(fp);
		/* dat.T = last - dat.startpos + 1; */
		/*fprintf(stderr, "%d lines\n", dat.T);*/
		rewind(fp);
		readseqs(fp, pseqdat[i]);
		fclose(fp);
	}


	/* Main option: est or vit */
	if (strcmp(argv[optind], "est") == 0){
	/* Baum-Welch estimation of HMM model parameters */
		double	***alpha;
		double	***beta;
		double	***gamma;

		int	niter;
		double	logprobinit, logprobfinal;


		/* allocate memory */
		alpha = calloc(seqnum, sizeof(double **));
		beta = calloc(seqnum, sizeof(double **));
		gamma = calloc(seqnum, sizeof(double **));
		for (i = 0; i < seqnum; i++) {
			alpha[i] = dmatrix(1, pseqdat[i]->T, 1, hmm.N);
			beta[i] = dmatrix(1, pseqdat[i]->T, 1, hmm.N);
			gamma[i] = dmatrix(1, pseqdat[i]->T, 1, hmm.N);
		}
		BaumWelch(&hmm, pseqdat, seqnum, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal, getprob);

		if (vflg) {
			fprintf(stderr, "Number of iterations: %d\n", niter);
			fprintf(stderr, "ln P(obs | init model): %E\n",
				logprobinit);	
			fprintf(stderr, "ln P(obs | estimated model): %E\n",
				logprobfinal);	
		}

		/* print the answer */
		PrintHMM(stdout, &hmm);
		/* if (getprob == 0)
			PrintHMM(stdout, &hmm);
		if (getprob == 1) {
			for (i=1; i <= numsites; i++){
				fprintf(stdout, "%ld", sites[i]);
				for (j=1; j <= hmm.N; j++)
					fprintf(stdout, "\t%g", gamma[1][sites[i] - pseqdat[0]->startpos + 1][j]);
				fprintf(stdout, "\n");
			}
		}*/


		/* free memory */
		for (i = 0; i < seqnum; i++) {
			free_dmatrix(alpha[i], 1, pseqdat[i]->T, 1, hmm.N);
			free_dmatrix(beta[i], 1, pseqdat[i]->T, 1, hmm.N);
			free_dmatrix(gamma[i], 1, pseqdat[i]->T, 1, hmm.N);
			freeseqdat(pseqdat[i]);
		}
	}
	else if (strcmp(argv[optind], "vit") == 0){
	/* calculation of Viterbi sequence */
		uint8_t *q;	/* state sequence q[1..T] */
		double 	/* proba, */ logproba;
		// int j;
		int64_t t;
		//double deltal;
		//if (getprob==1)
		//	fprintf(stderr, "Use est mode to get posterior probability from forward-backward algorithm. Running Viterbi algorithm now.\n");
		q = cvector(1, pseqdat[0]->T);
		Viterbi(&hmm, pseqdat[0], q, &logproba);

		/* for (t = 1; t <= dat.T; t++)
			fprintf(stdout, "%d\n", q[t]); */

		fprintf(stdout, "#start(inc)\tend(exc)\n");
		int in_arc = 0;
		for (t = 1; t < pseqdat[0]->T; t++){
			if (in_arc == 0 && q[t] == 2){
				fprintf(stdout, "%ld\t", pseqdat[0]->pos[t]);
				in_arc = 1;
			}
			if (in_arc == 1 && q[t] == 1){
				fprintf(stdout, "%ld\n", pseqdat[0]->pos[t - 1] + 1);
				in_arc = 0;
			}
		}
		if (in_arc == 1 && q[t] == 2)
			fprintf(stdout, "%ld\n", pseqdat[0]->pos[t]);

		free_cvector(q, 1, pseqdat[0]->T);
		for (i = 0; i < seqnum; i++)
			freeseqdat(pseqdat[i]);
	}
	else if (strcmp(argv[optind], "bfgs") == 0){
		trainbfgs(&hmm, pseqdat, seqnum);
	}
	else if (strcmp(argv[optind], "post") == 0){
		getprob = 1;
		double	***alpha;
		double	***beta;
		double	***gamma;

		int	niter;
		int64_t t;
		double	logprobinit, logprobfinal;


		/* allocate memory */
		alpha = calloc(seqnum, sizeof(double **));
		beta = calloc(seqnum, sizeof(double **));
		gamma = calloc(seqnum, sizeof(double **));
		for (i = 0; i < seqnum; i++) {
			alpha[i] = dmatrix(1, pseqdat[i]->T, 1, hmm.N);
			beta[i] = dmatrix(1, pseqdat[i]->T, 1, hmm.N);
			gamma[i] = dmatrix(1, pseqdat[i]->T, 1, hmm.N);
		}

		BaumWelch(&hmm, pseqdat, 1, alpha, beta, gamma, &niter, &logprobinit, &logprobfinal, getprob);

		for (t = 1; t <= pseqdat[0]->T; t++)
			fprintf(stdout, "%ld\t%g\n", pseqdat[0]->pos[t], gamma[0][t][2]);

		/* free memory */
		for (i = 0; i < seqnum; i++) {
			free_dmatrix(alpha[i], 1, pseqdat[i]->T, 1, hmm.N);
			free_dmatrix(beta[i], 1, pseqdat[i]->T, 1, hmm.N);
			free_dmatrix(gamma[i], 1, pseqdat[i]->T, 1, hmm.N);
			freeseqdat(pseqdat[i]);
		}

	}
	FreeHMM(&hmm);

	/* if (usemap == 1)
		freemap(&mdat);
	if (getprob == 1)
		free_livector(sites, 1, numsites);
	return 0; */
}

void Usage(char *name)
{
  printf("Usage: %s -m model.hmm est|vit|post|bfgs data.txt\n", name);
  printf("  model.hmm - file with initial model parameters\n");
  printf("  data - tab-seperated data file, columns: position, type(1,2,3,4), gendist\n");
  printf("Modes:\n");
  printf("  est: BW training\n");
  printf("  vit: get viterbi sequences of archaic segments\n");
  printf("  post: get posterior probability of archaic state at each site\n");
  printf("  bfgs: L-BFGSB training\n");
}
