#include "nrutil.h"
#include "ge.h"
#include<stdio.h>
#include<stdlib.h>
#include<stdbool.h>
#include<math.h>

/*
int main(int argc, char **argv){
    double **A, *b, *x;
    int n = 3;
    A = dmatrix(1,n,1,n);
    b = dvector(1,n);
    x = dvector(1,n);
    A[1][1] = 1; A[1][2]= 1; A[1][3] = 2;
    A[2][1] = 2; A[2][2]= 4; A[2][3] = -3;
    A[3][1] = 3; A[3][2]= 6; A[3][3] = -5;

    b[1] = 9; b[2] = 1; b[3] = 0;
    mprint(A,n,"A original");
    vprint(b,n,"b original");

    gauss_elim(A,b,x,n);
    mprint(A,n,"A");
    vprint(b,n,"b");
    vprint(x,n,"x");
}*/

void mprint(double **matrix, int m, char *label){
    int i, j,k;
    printf("%s:\n",label);

    for (i = 1; i <= m; ++i){
        for (j = 1; j <= m; ++j){
            printf("%10.2f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n------------------------\n");
}

void vprint(double *vector, int m, char *label){
    int i, j;
    printf("%s:\n",label);

    for (i = 1; i <= m; ++i){
        printf("%10.2f ", vector[i]);
    }

    printf("\n------------------------\n");
}

bool gauss_elim(double **A, double *b, double *x, int m)
{
    if (!upper_triangulate(A,b,m)) { /* first find the upper tri */
        return false;               /* problem, maybe it was singular */
    }
    if (!back_sub(A,x,b,m)) {       /* solve using back sub */
        return false;               /* problem, somehow pivot was zero? */
    }
    return true;
}


bool upper_triangulate(double **A, double *b, int m){
    int i,j,k,l;
    double scale;
    for (j=1;j<m;++j){           /* loop over columns */
        for (i=j+1;i<=m;++i){      /* loop over rows beneath pivot */
            int max = i-1;
            for (l=i;l<m;l++) {     /* go to each row looking for max */
                if (A[l][j] > A[max][j]) {  /* found a new max */
                    max = l;                /* save the index */
                }
            }
            /*if (max != (i-1) ) {            
               double *t1 = A[max];        
               A[max] = A[i-1];
               A[i-1] = t1;
               double temp = b[i-1];      
               b[i-1] = b[max];
               b[max] = temp;
            }*/
            if (A[i][i] == 0.0)
            {
                printf("Singular matrix\n");
                return false;
            }
            if (A[i][j] != 0) {       /* if entry not zero already */
                scale = A[i][j]/A[j][j];  /* zero out based on pivot */
                for (k=1;k<=m;++k) {
                    A[i][k] = A[i][k] - A[j][k]*scale;
                }
                b[i] = b[i] - b[j]*scale; /* same for b */
            }
        }
    }
    return true;
}

bool back_sub(double **A, double *x, double *b, int m){
    int i,j;
    x[m] = b[m]/A[m][m];

    for (i=m-1;i>=1;--i){
        x[i] = b[i];
        for (j=i+1;j<=m;++j) {
            x[i] -= A[i][j]*x[j];
        }
        if (A[i][i] == 0) {
            return false;
        }
        x[i]/=A[i][i];
    }
}
