#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

int main(){
    srandom(20252);
    int n=10;
    int k=3;
    double* A = (double*) malloc(n*n*sizeof(double));
    double* b = (double*) malloc(n*sizeof(double));
    double* X = (double*) calloc(n, sizeof(double));
    struct LinearSis SL = {A, b, n, k};
    double time = timestamp();

    genKDiagonal(&SL);

    printSis(&SL);

    return 0; 
}
