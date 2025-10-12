#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>
#include <likwid.h>

#include "utils.h"
#include "sislin.h"

int main(){
    LIKWID_MARKER_INIT;
    srandom(20252);
    int status = 0;
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

    LIKWID_MARKER_START("INICIO");
    /*Criando o sistema Linear, A: Matriz, b: Vetor Indep*/
    double *av = malloc(n*n*sizeof(double));
    struct Matrix A = {av, n, n, k};

    double *bv = malloc(n*sizeof(double));
    struct Matrix b = {bv, n, 1, 0};

    struct LinearSis SL = {&A, &b, n, k};
    genKDiagonal(&SL);
    //printSis(&SL);

    /*Gerando simetrica positiva, criando Matriz A2, b2 para guardar os novos valores*/
    double *av1 = malloc(n*n*sizeof(double));
    struct Matrix A2 = {av1, n, n, k};

    double *bv1 = malloc(n*sizeof(double));
    struct Matrix b2 = {bv1, n, 1, 0};
    double timePC; //Arrumar esse time e os outros do trabalho
    genSymmetricPositive(&SL, &A2, &b2, &timePC); //TODO: TALVEZ FAZER UMA VERIFICACAAO PARA EVITAR EM MATRIZES MAL CONDICIONADAAS

    free(av);
    free(bv); //Liberando o vetor A e b pois nao preciso mais deles
    
    SL.A = &A2;
    SL.b = &b2; //Agora os valores a matriz e vetor independente sao simetricos
    printSis(&SL);

    double* X = (double*) calloc(n, sizeof(double));
    double* r = (double*) malloc(n * sizeof(double));
    double* norma = (double*) malloc(n * sizeof(double));
    double timeM; 
    double timeGrad;
    double timeRes;
    //TODO: Fazer uma checagem para w
    //(nao comparar igualdade com valores ponto flutuante, fazer uma margem)
    
    double *Mv = calloc(n, sizeof(double));
    struct Matrix M = {Mv, n, 1, SL.k};

    genPreCond(SL.A, w, SL.n, SL.k, &M, &timeM);
    status = conjGradientPre(&SL, X, r, norma, &M, maxit, eps, &timeGrad);
    if (status == -1)
        return -1;
    free(Mv);
    calcResidue(&SL, X, r, &timeRes);
    
    
    printf("output:\n%d\n",n);
    
    double normaR = calcNormaEuclidiana(r, n);
    printf("%.8g\n", *norma);
    printVetor(X,n);
    printf("%.8g\n", normaR);
    printf("%.8g\n", timePC + timeM);
    printf("%.8g\n", timeGrad);
    printf("%.8g\n", timeRes);


    free(X);
    free(av1);
    free(bv1);
    free(norma);
    free(r); //Valgrind nao esta dando erro

    LIKWID_MARKER_STOP("INICIO");
    LIKWID_MARKER_CLOSE;
    return 0;
}
