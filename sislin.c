#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "sislin.h"

/* Cria matriz 'A' k-diagonal e Termos independentes B */
void genKDiagonal(struct LinearSis *SL){
    uint k = SL->k;
    uint n = SL->n;

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

void conjGradient(struct LinearSis *SL, double *X, double *R, uint maxit, double eps){
    

}

double static norma(double* x, uint n){
    double norma = 0;
    for (uint i = 0; i < n; i++){
        norma += x[i]*x[i];
    }
    return sqrt(norma);
}

double calcResidue(struct LinearSis *SL, double *X, double *time)
{
    uint n = SL->n;
    *time = timestamp();
    double sum = 0;
    double *r = calloc(n, sizeof(double));

    for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++)
            sum += SL->A[n*i + j] * X[j];
    
        r[i] = SL->b[i] - sum;
        printf("%f\n", r[i]);
    }


    *time = timestamp() - *time;

    return 0.0;
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

