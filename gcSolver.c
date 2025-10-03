#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

void main(){
int n=10;
int k=3;
double* A = (double*) malloc(n*n*sizeof(double));
double* B = (double*) malloc(n*sizeof(double));

srandom(20252);
criaKDiagonal(n,k,A,B);

prnsis(A,B,n);
}
