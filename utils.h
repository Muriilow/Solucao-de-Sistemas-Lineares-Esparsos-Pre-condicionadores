#ifndef __UTILS_H__
#define __UTILS_H__

#include <stdlib.h>
#include <time.h>
#include <sys/time.h>

// Valor absoluto de um número. Alternativa ao uso da função 'fabs()'
#define ABS(num)  ((num) < 0.0 ? -(num) : (num))

// Número máximo de dígitos em um número
#define numDigits(n)  6  // ( (int) log10(n) + 1 )

// Macro para verificar se valor 'n' é potência de 2 ou não
#define isPot2(n) (n && !(n & (n - 1)))

// Funções
double timestamp(void);
char* markerName(char* baseName, int n);

#endif // __UTILS_H__

