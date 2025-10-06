#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "sislin.h"


static double norma(double* x, uint n){
    double norma = 0;
    for (uint i = 0; i < n; i++){
        norma += x[i]*x[i];
    }
    return sqrt(norma);
}
/* Cria matriz 'A' k-diagonal e Termos independentes B */
void genKDiagonal(struct LinearSis *SL){
    int k = SL->k;
    int n = SL->n;

    for(int i = 0; i < n; i++){
        SL->b[i] = genRandomB(k);

        for(int j = 0; j < n; j++){
            if ((j > i+k/2 || j < i - k/2)) {
                SL->A[n*i+j] = 0.0;
                continue;
            }

            SL->A[n*i+j] = genRandomA(i,j,k);
        }
    }
}

/* Gera matriz simetrica positiva */
void genSymmetricPositive(double *A, double *b, int n, int k, double **ASP, double **bsp, double *time)
{
    *time = timestamp();
    *time = timestamp() - *time;

}


void genDLU (double *A, int n, int k, double **D, double **L, double **U, double *time)
{
    *time = timestamp();


    *time = timestamp() - *time;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(double *D, double *L, double *U, double w, int n, int k,
        double **M, double *time)
{
    *time = timestamp();


    *time = timestamp() - *time;
}


void genTranspose(struct LinearSis *SL, struct LinearSis* SLT)
{
    uint n = SL->n;

    for (uint i = 0; i < n; i++)
        for(uint j = 0; j < n; j++)
            SLT->A[j*n+i] = SL->A[i*n+j];
}

void conjGradient(struct LinearSis *SL, double *x, double *r, uint maxit, double eps){
    calcResidue(SL, x, r);

    uint n = SL->n;
    double *d = r; 
    double tam_passo; 

    uint it = 0;
    while (it < maxit) {
        tam_passo = norma(r, n) / norma(d, n); 
    }
}

void calcResidue(struct LinearSis *SL, double *x, double *r)
{
    uint n = SL->n;
    double sum = 0;

    for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++)
            sum += SL->A[n*i + j] * x[j];
    
        r[i] = SL->b[i] - sum;
        printf("%f\n", r[i]);
    }
}

void printSis(struct LinearSis *SL){
    uint n = SL->n;

    for (int i = 0; i < n; i++)
    {
        printf("  [  ");
        for(int j = 0; j < n; j++){
            if (SL->A[i*n+j] == 0)
                printf("            ");
            else
                printf("%.4e  ", SL->A[i*n+j]);
        }
        printf("]  [ %.4e ]\n", SL->b[i]);
    }
}

void multMatrix(struct Matrix *A, struct Matrix *B, struct Matrix *C) {
    if(A->column != B->row)
        return; 

    double sum = 0.0;
    uint aSize = A->rowSize;
    uint bSize = B->rowSize;

    for (uint i = 0; i < A->row - 1; i++) {
        for (uint j = 0; j < B->column - 1; j++) {
            sum = 0.0;
            for (uint k = 0; k < A->column -1; k++)
                sum += A->v[i*aSize + k] * B->v[k*bSize + j];

            C->v[i*A->row + j] = sum;
        }
    }
}

