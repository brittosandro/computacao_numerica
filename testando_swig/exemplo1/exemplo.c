/* File : exemplo.c */

#include "exemplo.h"
#include <stdio.h>
#include <time.h>

unsigned long fact (unsigned long n) {

    if (n < 2)
        return 1;
    else
        return n*fact(n-1);
}

void exec_fact (int n) {

    int cont = 0;
    clock_t tempo_inicial;
    clock_t tempo_final;

    while (cont <= n) {
        tempo_inicial = clock();
        printf("%d! = %ld \n", cont, fact(cont));
        tempo_final = clock();
        cont ++;
    }
    printf("%e \n",((double)(tempo_final - tempo_inicial) / CLOCKS_PER_SEC));


}
