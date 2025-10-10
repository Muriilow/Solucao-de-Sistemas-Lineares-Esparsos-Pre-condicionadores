#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

int main(){
    srandom(20252);
    int n=5;
    int k=3;
    uint maxit = 100;
    double eps = 1e-6;

    double a1[25] = {2,1,0,0,0, 1,2,1,0,0, 0,1,2,1,0, 0,0,1,2,1, 0,0,0,1,2};
    struct Matrix A = {a1, n, n, k};

    double b1[5] = {2,3,5,7,11};
    struct Matrix b = {b1, n, 1, 0};

    struct LinearSis SL = {&A, &b, n, k};
    double* X = (double*) calloc(n, sizeof(double));
    double* r = (double*) malloc(n * sizeof(double));
    double time = timestamp();

    
    conjGradient(&SL, X, r, maxit, eps);
    printSis(&SL);

    free(X);
    free(r);
    return 0;
    
    /* teste
     * srandom(20252);
    int n=3;
    int k=3;
    double v1[9] = {3, 4, 5, 1, 2, 3, 5, 6, 7};
    double v2[9] = {-1, 0, 3, 2, 1, 9, 9, 9, 9};
    double* v3 = (double*) malloc(n*n*sizeof(double));

    struct Matrix B = {v2, n, n, n, n};
    struct Matrix C = {v3, n, n, n, n};
    struct LinearSis SL = {v3, v3, n, k};

    double time = timestamp();
   
    multMatrix(&A, &B, &C);
    printSis(&SL);

    return 0; */
}
