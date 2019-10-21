/*
**      File:   nrutil.c
**      Purpose: Memory allocation routines borrowed from the
**		book "Numerical Recipes" by Press, Flannery, Teukolsky,
**		and Vetterling. 
*/

#include <malloc.h>
#include <stdio.h>
#include <stdint.h>

void nrerror(error_text)
char error_text[];
{
	void exit();

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(nl,nh)
int64_t nl,nh;
{
	float *v;

	v=(float *)calloc((unsigned) (nh-nl+1),sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

uint8_t *cvector(nl,nh)
int64_t nl,nh;
{
	uint8_t *v;

	v=(uint8_t *)calloc((unsigned) (nh-nl+1),sizeof(uint8_t));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl;
}

uint16_t *sivector(nl,nh)
int64_t nl,nh;
{
	uint16_t *v;

	v=(uint16_t *)calloc((unsigned) (nh-nl+1),sizeof(uint16_t));
	if (!v) nrerror("allocation failure in sivector()");
	return v-nl;
}

int64_t *livector(nl,nh)
int64_t nl,nh;
{
	int64_t *v;

	v=(int64_t *)calloc((unsigned) (nh-nl+1),sizeof(int64_t));
	if (!v) nrerror("allocation failure in livector()");
	return v-nl;
}

int *ivector(nl,nh)
int64_t nl,nh;
{
	int *v;

	v=(int *)calloc((unsigned) (nh-nl+1),sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(nl,nh)
int64_t nl,nh;
{
	double *v;

	v=(double *)calloc((unsigned) (nh-nl+1),sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}

float **matrix(nrl,nrh,ncl,nch)
int64_t nrl,nrh,ncl,nch;
{
	int i;
	float **m;

	m=(float **) calloc((unsigned) (nrh-nrl+1),sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) calloc((unsigned) (nch-ncl+1),sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int64_t nrl,nrh,ncl,nch;
{
	int64_t i;
	double **m;

	m=(double **) calloc((unsigned) (nrh-nrl+1),sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
		m[i] -= ncl;
	}
	return m;
}

uint8_t **cmatrix(nrl,nrh,ncl,nch)
int64_t nrl,nrh,ncl,nch;
{
	int64_t i;
	uint8_t **m;

	/*fprintf(stderr, "float: %u bytes, float*: %u bytes\n", sizeof(float), sizeof(float*));
	fprintf(stderr, "uint8: %u bytes, uint8*: %u bytes\n", sizeof(uint8_t), sizeof(uint8_t*));
	fprintf(stderr, "uint16: %u bytes, uint16*: %u bytes\n", sizeof(uint16_t), sizeof(uint16_t*));
	fprintf(stderr, "reserving %u bytes\n", (nrh-nrl+1)* (sizeof(uint8_t*)+ (nch-ncl+1)*sizeof(uint8_t)));*/

	m=(uint8_t **)calloc((unsigned) (nrh-nrl+1),sizeof(uint8_t*));
	if (!m) nrerror("allocation failure 1 in cmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(uint8_t *)calloc((unsigned) (nch-ncl+1),sizeof(uint8_t));
		if (!m[i]) nrerror("allocation failure 2 in cmatrix()");
		m[i] -= ncl;
	}
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
int64_t nrl,nrh,ncl,nch;
{
	int64_t i;
	int **m;

	m=(int **)calloc((unsigned) (nrh-nrl+1),sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)calloc((unsigned) (nch-ncl+1),sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
		m[i] -= ncl;
	}
	return m;
}

double ***d3darray(nrl,nrh,ncl,nch,nll,nlh)
int64_t nrl,nrh,ncl,nch,nll,nlh;
{
	int64_t i, j;
	double ***m;

	m=(double ***) calloc((unsigned) (nlh-nll+1),sizeof(double**));
	if (!m) nrerror("allocation failure 1 in threeDdmatrix()");
	m -= nll;

	for(i=nll;i<=nlh;i++) {
		m[i]=(double **) calloc((unsigned) (nrh-nrl+1),sizeof(double*));
		if (!m[i]) nrerror("allocation failure 2 in threeDdmatrix()");
		m[i] -= nrl;

		for(j=nrl;j<=nrh;j++) {
			m[i][j]=(double *) calloc((unsigned) (nch-ncl+1),sizeof(double));
			if (!m[i][j]) nrerror("allocation failure 3 in threeDdmatrix()");
			m[i][j] -= ncl;
		}
	}
	return m;
}

void free_vector(v,nl,nh)
float *v;
int64_t nl,nh;
{
	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v;
int64_t nl,nh;
{
	free((char *) (v+nl));
}

void free_cvector(v,nl,nh)
uint8_t *v;
int64_t nl,nh;
{
	free((char*) (v+nl));
}

void free_sivector(v,nl,nh)
uint16_t *v;
int64_t nl,nh;
{
	free((char*) (v+nl));
}

void free_livector(v,nl,nh)
int64_t *v;
int64_t nl,nh;
{
	free((char*) (v+nl));
}


void free_dvector(v,nl,nh)
double *v;
int64_t nl,nh;
{
	free((char*) (v+nl));
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int64_t nrl,nrh,ncl,nch;
{
	int64_t i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_cmatrix(m,nrl,nrh,ncl,nch)
uint8_t **m;
int64_t nrl,nrh,ncl,nch;
{
	int64_t i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int64_t nrl,nrh,ncl,nch;
{
	int64_t i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

/* void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int64_t nrl,nrh,ncl,nch;
{
	int64_t i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}*/

void free_d3darray(m,nrl,nrh,ncl,nch,nll,nlh)
double ***m;
int64_t nrl,nrh,ncl,nch,nll,nlh;
{
	int64_t i, j;

	for(i=nlh;i>=nll;i--){
		for(j=nrh; j>=nrl; j--) {
			free((char*) (m[i][j]+ncl));
		}
		free((char*) (m[i] + nrl));
	}
	free((char*) (m+nll));
}
