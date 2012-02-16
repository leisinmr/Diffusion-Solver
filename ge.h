#include "nrutil.h"


void mprint(double **matrix, int m, char *label);
void vprint(double *vector, int m, char *label);

bool gauss_elim(double **A, double *b, double *x, int m);
bool upper_triangulate(double **A, double *b, int m);
bool back_sub(double **A, double *x, double *b, int m);

