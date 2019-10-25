/*
**      File:   nrutil.h
**      Purpose: Memory allocation routines borrowed from the
**              book "Numerical Recipes" by Press, Flannery, Teukolsky,
**              and Vetterling.
*/

#include <stdint.h>

float *vector();
float **matrix();
double *dvector();
double **dmatrix();
uint8_t *cvector();
uint8_t **cmatrix();
int *ivector();
uint16_t *sivector();
int64_t *livector();
int **imatrix();
double ***d3darray();
void free_vector();
void free_dvector();
void free_cvector();
void free_ivector();
void free_sivector();
void free_livector();
void free_matrix();
void free_dmatrix();
void free_cmatrix();
//void free_imatrix();
void free_d3darray();
void nrerror();
