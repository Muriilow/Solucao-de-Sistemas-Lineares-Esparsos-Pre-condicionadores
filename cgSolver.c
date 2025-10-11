#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

int main(){
    srandom(20252);
    int n;
    int k;
    uint maxit;
    float w;
    double eps;

    scanf("%d", &n);
    scanf("%d", &k);
    scanf("%f", &w);
    scanf("%d", &maxit);
    scanf("%lf", &eps);

    double *a1 = malloc(n*n*sizeof(double));
    struct Matrix A = {a1, n, n, k};

    double *b1 = malloc(n*sizeof(double));
    struct Matrix b = {b1, n, 1, 0};

    struct LinearSis SL = {&A, &b, n, k};
    genKDiagonal(&SL);

    printSis(&SL);
    double* X = (double*) calloc(n, sizeof(double));
    double* r = (double*) malloc(n * sizeof(double));
    double time = timestamp();

    if(w == -1) { 
        conjGradient(&SL, X, r, maxit, eps);
    }
    else if(w == 0){
        double time; 
        double *Mv = calloc(n, sizeof(double));
        struct Matrix M = {Mv, n, 1, SL.k};
        genPreCond(SL.A, w, SL.n, SL.k, &M, &time);
        conjGradientPre(&SL, X, r, &M, maxit, eps);
        free(Mv);
    }
    printf("%d\n",n);
    
    printVetor(X, n);
    printVetor(r, n);

    free(X);
    free(a1);
    free(b1);
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
