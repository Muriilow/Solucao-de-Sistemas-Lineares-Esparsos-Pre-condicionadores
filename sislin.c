#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "utils.h"
#include "sislin.h"

/* Cria matriz 'A' k-diagonal e Termos independentes B */
void criaKDiagonal(int n, int k, double *A, double *B){
    for(int i = 0; i < n; i++){
        B[i] = generateRandomB(k);

        for(int j = 0; j < n; j++){
            if ((j > i+k/2)||(j < i-k/2))
                A[n*i+j] = 0.0;
            else 
                A[n*i+j] = generateRandomA(i,j,k);
        }
    }
}

/* Gera matriz simetrica positiva */
void genSimetricaPositiva(double *A, double *b, int n, int k, 
        double **ASP, double **bsp, double *tempo)
{
    *tempo = timestamp();

    *tempo = timestamp() - *tempo;

}


void geraDLU (double *A, int n, int k,
        double **D, double **L, double **U, double *tempo)
{
    *tempo = timestamp();


    *tempo = timestamp() - *tempo;
}

/**
 * Devolve matriz M⁻¹
 *
 */
void geraPreCond(double *D, double *L, double *U, double w, int n, int k,
        double **M, double *tempo)
{
    *tempo = timestamp();


    *tempo = timestamp() - *tempo;
}


double calcResiduoSL (double *A, double *b, double *X,
        int n, int k, double *tempo)
{
    *tempo = timestamp();

    double *r = calloc(n, sizeof(double));
    /*for (int i; i< ;) {
    }*/


    *tempo = timestamp() - *tempo;
}

void prnsis(double *A, double *B, int n){
    for (int i = 0; i < n; i++)
    {
        printf("  [  ");
        for(int j = 0; j < n; j++){
            if (A[i*n+j] == 0)
                printf("            ");
            else
                printf("%.4e  ", A[i*n+j]);
        }
        printf("]  [ %.4e ]\n", B[i]);
    }
}

