/*
 * map.c
 *
 *  Created on: 22 Aug 2016
 *      Author: ruoyun
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmmarc.h"

#define GFTIME 2000

void readmap(FILE *fmap, mapdat *pmap)
{
	int64_t i;
	int64_t locval;
	double rateval;

	fprintf(stderr, "map length: %ld lines\n", pmap->l);
	
	pmap->loc = (int64_t *) livector(1, pmap->l);
	pmap->rate = dvector(1, pmap->l);
	for (i=1; i <= pmap->l ; i++) {
	  if(fscanf(fmap, "%ld\t%lf\t%*f\n", &locval, &rateval)) {
		  pmap->loc[i] = locval;
		  pmap->rate[i] = rateval;
	  }
	}
}

void freemap(mapdat *pmap)
{
	free_livector(pmap->loc, 1, pmap->l);
	free_dvector(pmap->rate, 1, pmap->l);
}

void computeA(HMM *hmm, mapdat *pmap, int64_t pos) /* Here assume only 2 states in HMM; return A in linear space */
{
	int64_t i = 1;
	double recrate;
	while (pmap->loc[i] < pos)
		i++;
	recrate = pmap->rate[i] / 100000000 * GFTIME; // Time since admixture (2000) times gen dist.
	hmm->A[1][2] = recrate * hmm->pi[2];
	hmm->A[1][1] = 1 - hmm->A[1][2];
	hmm->A[2][1] = recrate * hmm->pi[1];
	hmm->A[2][2] = 1 - hmm->A[2][1];
}
