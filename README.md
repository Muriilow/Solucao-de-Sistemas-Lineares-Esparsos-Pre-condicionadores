# Solução de sistemas lineares esparsos usando Pré-condicionadores

Programa que utilizado o método dos Gradientes Conjugados para achar uma solução **x̄** dentro de uma certa tolerância **ε** ou limite de iterações **maxit**.

### Entradas: 
lê da entrada padrão (stdin) 5 cinco valores:

    n: (n>10)  a dimensão do Sistema Linear.
    k: (k>1 e k ímpar)  o número de diagonais da matriz A.
    ω:  o pré-condicionador a ser utilizado:

        ω=-1: sem pré-condicionador
        ω=0.0 : pré-condicionador de Jacobi

    maxit:  o número máximo de iterações a serem executadas.
    ε: o erro aproximado absoluto máximo da norma do vetor X.


## Autoria
Guilherme Vitoriano Santana de OLiveira - GRR20245396  
Murilo de Paula Bob - GRR2024

## Estruturas de Dados:

Struct que guarda uma matriz e informações importantes par trabalhar sobre ela:

```
struct Matrix {
    double *v;
    uint row; 
    uint column;
    uint k;
};
```
   
> v: Matriz diagonal  
row: Quantidade de linhas   
column: Quantidade de coluna  
k: Quatidade de diagonais da matriz, 0 para matriz não diagonal

 Struct que define um sistema Linear:
```
struct LinearSis {
    struct Matrix *A;
    struct Matrix *b;
    uint n; 
    uint k; 
};
```
> A: Matriz k-diagonal  
b: Vetor de valores independentes do sistema  
n: Ordem da matriz  
k: Valor que define quantas diagonais não zero


# Módulos
    utils
    sislin


