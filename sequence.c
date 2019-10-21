/*
**	File:	sequence.c
**	Purpose: Routines for reading data.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "nrutil.h"
#include "hmmarc.h"

#define MAX_UINT16 0xffff

int64_t countlines(FILE *fp){
	int64_t lines = 0;
	char ch, prev = '\n';

	while((ch = fgetc(fp)) != EOF){
		if (ch == '\n'){
			lines++;
		}
		if (ch == '#'){
			lines--;
		}
		prev = ch;
	}
	if (prev != '\n')
		lines++;
	return(lines);
}

void readseqs(FILE *fdat, seqdat *pdat)
{
/* read data records. expect: snp position, snp type (0, 1, 2, 3, 4) */
	int64_t t;
	int posval;
	float dstval;
	// int64_t endpos = pdat->startpos + pdat->T - 1;
	int snpval;

	//fprintf(stderr, "Begin reading observations. Range of positions: %ld - %ld\n", pdat->startpos, endpos);
	pdat->pos = ivector(1, pdat->T);
	pdat->snptype = (uint8_t *) cvector(1, pdat->T);
	pdat->gendist = dvector(1, pdat->T);
	/* while (fgetc(fdat) == '#')
		if(fgets(line, sizeof(line), fdat)){};
	fseek(fdat, -1, SEEK_CUR);
	for (t=1; t <= pdat->T; t++)
		pdat->snptype[t] = 4; */
	for (t=1; t <= pdat->T; t++) {
		if (fscanf(fdat,"%d\t%d\t\%f\n", &posval, &snpval, &dstval) == 3){
			pdat->pos[t] = posval;
			pdat->gendist[t] = dstval;
			pdat->snptype[t] = (uint8_t) snpval;
		}
	}
}

void freeseqdat(seqdat *pdat)
{
	free_dvector(pdat->gendist, 1, pdat->T);
	free_cvector(pdat->snptype, 1, pdat->T);
	free_ivector(pdat->pos, 1, pdat->T);
}

/*void readsites(FILE *fsites, int64_t *sites)
{
	int i = 1;
	int64_t val;
	while (fscanf(fsites, "%*s\t%ld\n", &val) == 1) {
		sites[i] = val;
		i++;
	}
}
*/
