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

    /*Criando o sistema Linear, A: Matriz, b: Vetor Indep*/
    double *av = malloc(n*n*sizeof(double));
    struct Matrix A = {av, n, n, k};

    double *bv = malloc(n*sizeof(double));
    struct Matrix b = {bv, n, 1, 0};

    struct LinearSis SL = {&A, &b, n, k};
    genKDiagonal(&SL);
    printSis(&SL);

    /*Gerando simetrica positiva, criando Matriz A2, b2 para guardar os novos valores*/
    double *av1 = malloc(n*n*sizeof(double));
    struct Matrix A2 = {av1, n, n, k};

    double *bv1 = malloc(n*sizeof(double));
    struct Matrix b2 = {bv1, n, 1, 0};
    double time2; //Arrumar esse time e os outros do trabalho
    genSymmetricPositive(&SL, &A2, &b2, &time2); //TODO: TALVEZ FAZER UMA VERIFICACAAO PARA EVITAR EM MATRIZES MAL CONDICIONADAAS

    free(av);
    free(bv); //Liberando o vetor A e b pois nao preciso mais deles
    
    SL.A = &A2;
    SL.b = &b2; //Agora os valores a matriz e vetor independente sao simetricos
    printSis(&SL);


    double* X = (double*) calloc(n, sizeof(double));
    double* r = (double*) malloc(n * sizeof(double));
    double time = timestamp();

    //TODO: Fazer uma checagem para w
    //(nao comparar igualdade com valores ponto flutuante, fazer uma margem)
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
    free(av1);
    free(bv1);
    free(r); //Valgrind nao esta dando erro
    return 0;
}
