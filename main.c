#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

void main(){
    int n = 10;
    int k = 5;

    double** A = (double**) malloc(n*n*sizeof(double));
    double* B = malloc (n*sizeof(double));

    criaKDiagonal(n, k, A, B);
    prnsis(A,B,n);
