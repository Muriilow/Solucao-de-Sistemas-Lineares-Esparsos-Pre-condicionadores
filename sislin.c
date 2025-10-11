#include <stdio.h>
#include <stdlib.h>    /* for exit e random/srandom */
#include <string.h>
#include <math.h>

#include "sislin.h"

/*TODO: DEVEMOS VERIFICAR DIVISAO POR ZERO NAS FUNCOES 
 NO ENUNCIADO BASTA UMA SIMPLES PRINT DE AVISO E RETORNAR O MAIN COM != 0 */

/**
 * Função que gera os coeficientes de um sistema linear k-diagonal
 * @param i,j coordenadas do elemento a ser calculado (0<=i,j<n)
 * @param k numero de diagonais da matriz A
 */
static inline double genRandomA(uint i, uint j, uint k){
    static double invRandMax = 1.0 / (double)RAND_MAX;
    //Eh sempre diagonal dominante
    return ( (i==j) ? (double)(k<<1) : 1.0 )  * (double)random() * invRandMax;
}

/**
 * Função que gera os termos independentes de um sistema linear k-diagonal
 * @param k numero de diagonais da matriz A
 */
static inline double genRandomB(uint k){
    static double invRandMax = 1.0 / (double)RAND_MAX;
    return (double)(k<<2) * (double)random() * invRandMax;
}

static double sqrVector(double* x, double* y, uint n){
    double sqrVector = 0;
    for (uint i = 0; i < n; i++){
        sqrVector += x[i]*y[i];
    }
    return sqrVector;
}

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

/* Gera matriz simetrica positiva. Isso pode ser ruim ou bom. A maneira que estamos usando (CGNE) pode piorar uma matriz mal condicionada. O OBJETIVO DESSA FUNCAO EH TRANSFORMAR UMA MATRIZ NAO SDP EM SDP, MAS CASO A MATRIZ SEJA MAL CONDICIONADA (NESSE CASO INDEPENDE DE SER OU NAO SIMETRICA E POSITIVA) ESSA FUNCAO IRA PIORAR EM MUITO A RESOLUCAO DA MATRIZ
 * PORTANTO, CHECAR COM O PROFESSOR SE EH NECESSARIO FAZER VERIFICACAO ANTES DE USAR A CGNE, PARA EVITAR ESSES PROBLEMAS*/
void genSymmetricPositive(struct LinearSis *SL, struct Matrix *ASP, struct Matrix *bsp, double *time)
{
    struct Matrix *A = SL->A;
    struct Matrix *b = SL->b;

    *time = timestamp();

    double *Atv = malloc(A->row*A->column*sizeof(double));
    struct Matrix AT ={Atv, A->row, A->column, A->k};

    genTranspose(&AT, A);
    
    //ESSE CALCULO EH O MESMO QUE ELEVAR O QUADRADO DE TODOS OS ELEMENTOS DE A
    //CASO ACHE INTERESSANTE FAZER UMA FUNCAO QUE FAZ SO ISSO SERIA MENOS CUSTOSO ACREDITO EU.
    multMatrix(A, &AT, ASP);
    multMatrix(&AT, b, bsp);

    *time = timestamp() - *time;
    free(Atv);
}

//TODO: Deu erro ultima vez que usei, nao sei se eh a funcao ou a forma como coloquei 
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

/*Um pre condicionamento melhora um SL simetrico, positivo, definido e mal condicionado.*/
void genPreCond(struct Matrix *A, double w, int n, int k,
        struct Matrix *M, double *time)
{   
    if (w == -1)
        for(int i = 0; i < M->row; i++){
            if(A->v[i*n + i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[%d][%d] = %.6f -> A não é definida positiva\n", i, i, A->v[i*n + i]);
                return; 
            }
            else
                M->v[i] = 1.0/A->v[i*n + i];
            return;
        }
        

    /*
    PELO VISTO GERAR A DLU PODE SER IGNORADA (CASO SO FIZERMOS O METODO JACOBI)
    TEM QUE VER COM O PROFESSOR
    double* dv = calloc(A->row*A->column, sizeof(double));
    double* lv = calloc(A->row*A->column, sizeof(double));
    double* uv = calloc(A->row*A->column, sizeof(double));
    struct Matrix D = {dv, A->row, A->column, A->k};
    struct Matrix L = {lv, A->row, A->column, A->k};
    struct Matrix U = {uv, A->row, A->column, A->k};
    double *DLUTime;

    genDLU(A,&D,&L,&U,DLUTime); TODO: Consertar o ERRO */

    *time = timestamp();
    if (w == 0){
        for(int i = 0; i < M->row; i++){
            if(A->v[i*n + i] <= 0) { //Matriz nao positiva definida 
                fprintf(stderr, "ERRO: A[%d][%d] = %.6f -> A não é definida positiva\n", i, i, A->v[i*n + i]);
                return; 
            }
            else
                M->v[i] = 1.0/A->v[i*n + i];
        }
    }

    *time = timestamp() - *time;
}

void genTranspose(struct Matrix *A, struct Matrix *T)
{
    uint n = A->column;
    uint m = A->row;

    for (uint i = 0; i < n; i++)
        for(uint j = 0; j < m; j++)
            A->v[j*n+i] = T->v[i*n+j];
}
/*TODO: COMO VOCE DEVE TER VISTO, EU CRIEI DOIS METODOS SIMILARES
 * UM PARA FAZER SEM PRECOND E OUTRO PARA FAZER COM PRECOND
 * PRECISAMOS DE UMA FORMA MAIS ELEGANTE DE RESOLVER ISSO COM APENAS UMA FUNCAO
 * TEM QUE VER SE O PROFESSOR DESCONTA NOTA DISSO, PROVAVELMENTE SIM!!!
 * */
void conjGradient(struct LinearSis *SL, double *x, double *r,double *norma, uint maxit, double eps, double* time){

    calcResidue(SL, x, r, NULL);
    uint n = SL->n;

    /*Criando a matriz d e c usados para calculos*/
    double *v1 = malloc(n * sizeof(double));
    memcpy(v1, r, n * sizeof(double)); 
    struct Matrix d = {v1, n, 1, 0};

    double *v2 = malloc(n * sizeof(double));
    struct Matrix c = {v2, n, 1, 0};
    
    double* prevx = calloc(n, sizeof(double));
    double cAd = 0.0; //dkt * Adk
    double alpha; // ak
    double deltaOld = 0.0;
    double deltaNew = 0.0;
    double beta = 0.0;
    double tIter;
    uint it = 1;
    double totalTime = 0.0;
    deltaOld = sqrVector(r, r, n); //Como rk * rkt eh o quadrado nao precisamos multplicar matrizes
    do {
        tIter = timestamp();
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

            prevx[i] = x[i];
            /*Xk+1 = Xk + akdk*/
            x[i] += alpha * d.v[i];
            
            //printf("%f - %d\n", x[i], i);
            /*rk+1 = rk - akAdk*/
            r[i] -= alpha * c.v[i];
            deltaNew += r[i] * r[i];*time = timestamp();
        }

        beta = deltaNew / deltaOld;
        for(uint i = 0; i < n; i++)
            d.v[i] = r[i] + beta *d.v[i];

        deltaOld = deltaNew;
        *norma = calcNormaMax(x, prevx, n);
        it++;
        tIter = timestamp() - tIter;
        totalTime +=tIter;
    }while (it < maxit &&  *norma >= eps);
    *time = totalTime/it;

    free(v1);
    free(v2);
}

void conjGradientPre(struct LinearSis *SL, double *x, double *norma, double *r, struct Matrix *M, uint maxit, double eps, double* time){

    calcResidue(SL, x, r, NULL);
    uint n = SL->n;

    // Y para calcular o SL com condicionador
    double *Yv = malloc(n * sizeof(double)); 
    struct Matrix y = {Yv, 1, SL->n, 0};
    struct Matrix rMatrix = {r, 1, SL->n, 0};
    for (int i = 0; i < n; i++)
        y.v[i] = M->v[i] * r[i]; // y = M^-1 * r

    /*Criando a matriz d e c usados para calculos*/
    double *v1 = calloc(n,sizeof(double));
    struct Matrix d = {v1, n, 1, 0};
    for (int i = 0; i < n; i++)
        d.v[i] = M->v[i] * SL->b->v[i];

    multMatrix(M, SL->b, &d); // v = M^-1 * b

    double *v2 = malloc(n * sizeof(double));
    struct Matrix c = {v2, n, 1, 0};

    double *prevx = calloc(n,sizeof(double));
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
        printf("alpha %f = deltaold %f / cAd %f\n", alpha, deltaOld,cAd);
         
        //printf("%f\n", alpha);
       
        deltaNew = 0.0;
        for (uint i = 0; i < n; i++) {
            prevx[i] = x[i];
            printf("vet antes:\n");
            printVetor(x,n);
            /*Xk+1 = Xk + akdk*/
            x[i] += alpha * d.v[i];
            printf("vet depois:\n");
            printVetor(x,n);
            printf("\n\n");
            //printf("%f - %d\n", x[i], i);
            /*rk+1 = rk - akAdk*/
            r[i] -= alpha * c.v[i];

            
            deltaNew += y.v[i] * r[i];
            valueNew += r[i] * r[i];
            norma[i] -= x[i];
            printf("deltanew: %f\ny.v[%d]: %f\nr[%d]:%f\n",deltaNew,i,y.v[i],i,r[i]);
        }

        printf("Bdeltanew %f\n", deltaNew);
        beta = deltaNew / deltaOld;
        for(uint i = 0; i < n; i++)
            d.v[i] = r[i] + beta *d.v[i];

        printf("deltanew %f\n", deltaNew);
        deltaOld = deltaNew;
        it++;
        *norma = calcNormaMax(x, prevx, n); 
        tIter = timestamp() - tIter;
        printf("it:%d  maxit: %d\n", it, maxit);
    }while (it < maxit && *norma >= eps); //ESSA VERIFICACAO DE COVNERGENCIA ESTA ERRADA, OLHAR ENUNCIADO
    *time = tIter/it;
    free(v1);
    free(v2);
}

double calcNormaMax(double *x,double* y, int n){
    double max = 0.0;
    double aux = 0.0;
    for (int i = 0; i < n; i++){
        aux = fabs(x[i] - y[i]);
        if (max < aux)
            max = aux;
    }
    printf("max: %f\n\n", max);
    return max;
}

double calcNormaEuclidiana(double *x, int n){
    double aux = 0.0;
    for (int i = 0; i < n; i++){
        aux += x[i]*x[i];
    }
    return sqrt(aux);
}

void calcResidue(struct LinearSis *SL, double *x, double *r, double *time)
{
    if (time)
        *time = timestamp();
    uint n = SL->n;
    double sum = 0;

    for (uint i = 0; i < n; i++) {
        for (uint j = 0; j < n; j++)
            sum += SL->A->v[n*i + j] * x[j];
    
        r[i] = SL->b->v[i] - sum;
        //printf("%f\n", r[i]);
    }
    if (time)
        *time = timestamp() - *time;
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
/*ESSA FUNCAO MAIS GERAL FUNCIONA PARA MATRIZES E VETORES, MAS DEPENDE DE COMO O VETOR ESTA ORGANIZADO
 * UM VETOR LINHA PODE CAUSAR PROBLEMAS ONDE UM VETOR COLUNA NAO, FAZER UMA MULT ENTRE MATRIZ E VETOR SERIA MAIS SEGURO*/
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

/*ME PARECE DESNECESSARIO, UTILIZAR A MULTIPLICACAO ENTRE MATRIZES JA  BASTA PARA MANTER SIMETRICA POSITIVA DEFINIDA, VER CGNE */
/*void sumMatrix(struct Matrix *A, struct Matrix *B, struct Matrix *C) {
    if(A->column != B->row)
        return; 

    uint n = A->column;

    for (uint i = 0; i < A->row; i++) {
        for (uint j = 0; j < B->column ; j++) {
            C->v[i*n+j] = A->v[i*n+j] + B->v[i*n+j];
        }
    }
}*/
