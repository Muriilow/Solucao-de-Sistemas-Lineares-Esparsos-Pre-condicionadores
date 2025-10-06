#ifndef __SISLIN_H__
#define __SISLIN_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "utils.h"


/*ATENCAO: Funcoes inline tem sua implementacao no .h*/
/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline double genRandomA(uint i, uint j, uint k)
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return ( (i==j) ? (double)(k<<1) : 1.0 )  * (double)random() * invRandMax;
};

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline double genRandomB(uint k)
{
  static double invRandMax = 1.0 / (double)RAND_MAX;
  return (double)(k<<2) * (double)random() * invRandMax;
};

void genKDiagonal(struct LinearSis *SL);

void genSymmetricPositive(double *A, double *b, int n, int k, double **ASP, double **bsp, double *time);
void genDLU(double *A, int n, int k, double **D, double **L, double **U, double *time);
void genPreCond(double *D, double *L, double *U, double w, int n, int k, double **M, double *time);
void genTranspose(struct LinearSis *SL, struct LinearSis *SLT);
void conjGradient(struct LinearSis *SL, double *x, double *r, uint maxit, double eps);
void calcResidue(struct LinearSis *SL, double *x, double *r);
void printSis(struct LinearSis *SL);

#endif // __SISLIN_H__

