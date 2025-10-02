#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"


/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline double generateRandomA( unsigned int i, unsigned int j, unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ( (i==j) ? (double)(k<<1) : 1.0 )  * (double)random() * invRandMax;
};

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline double generateRandomB( unsigned int k )
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return (double)(k<<2) * (double)random() * invRandMax;
};

void criaKDiagonal(int n, int k, double *A, double *B);

void genSimetricaPositiva(double *A, double *b, int n, int k, double **ASP, double **bsp, double *tempo);
void geraDLU (double *A, int n, int k, double **D, double **L, double **U, double *tempo);
void geraPreCond(double *D, double *L, double *U, double w, int n, int k, double **M, double *tempo);
double calcResiduoSL (double *A, double *b, double *X, int n, int k, double *tempo);
void prnsis(double *A, double *B, int k);

#endif // __SISLIN_H__

