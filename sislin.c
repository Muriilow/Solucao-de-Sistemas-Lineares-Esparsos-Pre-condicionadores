#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "sislin.h"


static double sqrVector(double* x, double* y, uint n){
    double sqrVector = 0;
    for (uint i = 0; i < n; i++){
        sqrVector += x[i]*y[i];
    }
    return sqrVector;
}
/* Cria matriz 'A' k-diagonal e Termos independentes B */
void genKDiagonal(struct LinearSis *SL){
    int k = SL->k;
    int n = SL->n;

    for(int i = 0; i < n; i++){
        SL->b->v[i] = genRandomB(k);

        for(int j = 0; j < n; j++){
            if ((j > i+k/2 || j < i - k/2)) {
                SL->A->v[n*i+j] = 0.0;
                continue;
            }

            SL->A->v[n*i+j] = genRandomA(i,j,k);
        }
    }
}

/* Gera matriz simetrica positiva */
void genSymmetricPositive(struct Matrix *A, struct Matrix *b, int n, int k, struct Matrix *ASP, double *bsp, double *time)
{
    *time = timestamp();
    struct Matrix AT ={malloc(A->row*A->column*sizeof(double)), A->column, A->row, A->k};

    genTranspose(A, &AT);

    sumMatrix(A,&AT,ASP);
    *time = timestamp() - *time;
}


void genDLU(struct Matrix* A, struct Matrix* D, struct Matrix* L, struct Matrix* U, double *time)
{
    int n = A->row;
    int k = A->k;
    *time = timestamp();

    for(int i = 0; i < n; i++){
        for(int j = i-k/2; j < i+k/2; j++){
            if(i < j)
                L->v[i*n+j] = A->v[i*n+j];
            else if(i == j)
                D->v[i*n+j] = A->v[i*n+j];
            else
                U->v[i*n+j] = A->v[i*n+j];
        }
    }


    *time = timestamp() - *time;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void genPreCond(struct Matrix *A, double w, int n, int k,
        struct Matrix *M, double *time)
{   
    if (w == -1)
        return;

    double* dv = calloc(A->row*A->column, sizeof(double));
    double* lv = calloc(A->row*A->column, sizeof(double));
    double* uv = calloc(A->row*A->column, sizeof(double));
    struct Matrix D = {dv, A->row, A->column, A->k};
    struct Matrix L = {lv, A->row, A->column, A->k};
    struct Matrix U = {uv, A->row, A->column, A->k};
    double *DLUTime;

    //genDLU(A,&D,&L,&U,DLUTime); TODO: Consertar o ERRO 

    *time = timestamp();
    if (w == 0){
        for(int i = 0; i < M->row; i++){
            M->v[i] = 1.0/A->v[i*n + i];
        }
    } else if (w == 1){
        ; //placeholder
    } else if (w > 1 && w < 2){
        ; //placeholder
    }

    *time = timestamp() - *time;
}
void inverseDLU(){

};


void genTranspose(struct Matrix *A, struct Matrix *T)
{
    uint n = A->column;
    uint m = A->row;

    for (uint i = 0; i < n; i++)
        for(uint j = 0; j < m; j++)
            A->v[j*n+i] = T->v[i*n+j];
}

void conjGradient(struct LinearSis *SL, double *x, double *r, uint maxit, double eps){
    calcResidue(SL, x, r);

    uint n = SL->n;

    /*Criando a matriz d e c usados para calculos*/
    double *v1 = malloc(n * sizeof(double));
    memcpy(v1, r, n * sizeof(double)); 
    struct Matrix d = {v1, n, 1, 0};

    double *v2 = malloc(n * sizeof(double));
    struct Matrix c = {v2, n, 1, 0};

    double cAd = 0.0; //dkt * Adk
    double alpha; // ak
    double deltaOld = 0.0;
    double deltaNew = 0.0;
    double beta = 0.0;
    double tIter = timestamp();
    uint it = 1;
    deltaOld = sqrVector(r, r, n); //Como rk * rkt eh o quadrado nao precisamos multplicar matrizes
    do {
        /*Calculando ak = rk * rkt / dkt * A * dk */
        multMatrix(SL->A, &d, &c); //Precisamos multiplicar matrizes pois A e d nao sao quadrados
       
        //Fazemos a multiplicacao de matrizes manual, ja q n vou criar outro vetor
        cAd = 0.0;
        for (int i = 0; i < n; i++) 
            cAd += c.v[i]*d.v[i];

        alpha = deltaOld / cAd; //Calculando ak
         
        //printf("%f\n", alpha);
       
        deltaNew = 0.0;
        for (uint i = 0; i < n; i++) {
            /*Xk+1 = Xk + akdk*/
            x[i] += alpha * d.v[i];
            
            //printf("%f - %d\n", x[i], i);
            /*rk+1 = rk - akAdk*/
            r[i] -= alpha * c.v[i];
            deltaNew += r[i] * r[i];
        }

        beta = deltaNew / deltaOld;
        for(uint i = 0; i < n; i++)
            d.v[i] = r[i] + beta *d.v[i];

        deltaOld = deltaNew;
        it++;
        tIter -= timestamp();
    }while (it < maxit && sqrt(deltaNew) >= eps);
    
    printf("%f - %f\n", deltaOld, cAd);

    free(v1);
    free(v2);
}

void conjGradientPre(struct LinearSis *SL, double *x, double *r, struct Matrix *M, uint maxit, double eps){

    double time;
    calcResidue(SL, x, r);
    uint n = SL->n;
    
    double *Yv = malloc(n * sizeof(double)); // Y para calcular o SL com condicionador
    struct Matrix y = {Yv, 1, SL->n, 0};
    struct Matrix rMatrix = {r, 1, SL->n, 0};
    for (int i = 0; i < n; i++) {
        y.v[i] = M->v[i] * r[i];
        //printf("%lf - %lf - %lf\n", y.v[i], M->v[i], r[i]);
    }
    //multMatrix(M, &rMatrix, &y); // y = M-1 * r

    /*Criando a matriz d e c usados para calculos*/
    double *v1 = calloc(n,sizeof(double));
    struct Matrix d = {v1, n, 1, 0};
    for (int i = 0; i < n; i++) {
        d.v[i] = M->v[i] * SL->b->v[i];
        //printf("%lf - %lf - %lf\n", d.v[i], M->v[i], SL->b->v[i]);
    }
    multMatrix(M, SL->b, &d); // v = M-1 * b

    double *v2 = malloc(n * sizeof(double));
    struct Matrix c = {v2, n, 1, 0};

    double cAd = 0.0; //dkt * Adk
    double alpha; // ak
    double deltaOld = 0.0;
    double deltaNew = 0.0;
    double valueNew = 0.0;
    double beta = 0.0;
    double tIter = timestamp();
    uint it = 1;
    deltaOld = sqrVector(r, y.v, n); //Como rk * rkt eh o quadrado nao precisamos multplicar matrizes
    do {
        /*Calculando ak = rk * rkt / dkt * A * dk */
        multMatrix(SL->A, &d, &c); //Precisamos multiplicar matrizes pois A e d nao sao quadrados
       
        //Fazemos a multiplicacao de matrizes manual, ja q n vou criar outro vetor
        cAd = 0.0;
        for (int i = 0; i < n; i++) 
            cAd += c.v[i]*d.v[i];

        alpha = deltaOld / cAd; //Calculando ak
         
        //printf("%f\n", alpha);
       
        deltaNew = 0.0;
        for (uint i = 0; i < n; i++) {
            /*Xk+1 = Xk + akdk*/
            x[i] += alpha * d.v[i];
            
            //printf("%f - %d\n", x[i], i);
            /*rk+1 = rk - akAdk*/
            r[i] -= alpha * c.v[i];

            deltaNew += y.v[i] * r[i];
            valueNew += r[i] * r[i];
        }

        beta = deltaNew / deltaOld;
        for(uint i = 0; i < n; i++)
            d.v[i] = r[i] + beta *d.v[i];

        deltaOld = deltaNew;
        it++;
        tIter -= timestamp();
    }while (it < maxit && sqrt(valueNew) >= eps);
    
    free(v1);
    free(v2);
}

void calcResidue(struct LinearSis *SL, double *x, double *r)
{
    uint n = SL->n;
    double sum = 0;

    for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++)
            sum += SL->A->v[n*i + j] * x[j];
    
        r[i] = SL->b->v[i] - sum;
        //printf("%f\n", r[i]);
    }
}

void printVetor(double* vet, int n){
    for(int i = 0; i < n; i++){
        printf("%.16g   ", vet[i]);
    }
    printf("\n");
}
void printSis(struct LinearSis *SL){
    uint n = SL->n;

    for (int i = 0; i < n; i++)
    {
        printf("  [  ");
        for(int j = 0; j < n; j++){
            if (SL->A->v[i*n+j] == 0)
                printf("            ");
            else
                printf("%.4e  ", SL->A->v[i*n+j]);
        }
        printf("]  [ %.4e ]\n", SL->b->v[i]);
    }
}

void multMatrix(struct Matrix *A, struct Matrix *B, struct Matrix *C) {
    if(A->column != B->row)
        return; 

    double sum = 0.0;
    uint aSize = A->column;
    uint bSize = B->column;

    for (uint i = 0; i < A->row; i++) {
        for (uint j = 0; j < B->column ; j++) {
            sum = 0.0;
            for (uint k = 0; k < A->column; k++)
                sum += A->v[i*aSize + k] * B->v[k*bSize + j];
            C->v[i*bSize + j] = sum;
            //printf("%f - C[%d]\n", C->v[i*bSize + j], i*bSize + j);
        }
    }
}

void sumMatrix(struct Matrix *A, struct Matrix *B, struct Matrix *C) {
    if(A->column != B->row)
        return; 

    uint n = A->column;

    for (uint i = 0; i < A->row; i++) {
        for (uint j = 0; j < B->column ; j++) {
            C->v[i*n+j] = A->v[i*n+j] + B->v[i*n+j];
        }
    }
}
