#include <stdio.h>
#include <time.h>
#include "fatorial.h"

void fatorial(int n) {

    int a[1000]; // Vetor com capacidade de armazenar até 200 digitos
    int i, j, temp, m, x, t = 0; // variáveis de interesse

    while (t <= n) {
        a[0] = 1;  // Inicializa o vetor com um único digito cujo digito é 1.
        m = 1;     // segundo contador, Observe o primeiro é t!

        temp = 0;
        for (i = 1; i <= t; i++) {
           for (j=0;j<m;j++) {
               x = a[j]*i + temp;  // x é o produto
               a[j] = x % 10;      // armazenando na posição j
               temp = x / 10;
           }
           while (temp > 0) {
              a[m] = temp % 10;
              temp = temp / 10;
              m++;
           }
      }

    for(i=m-1;i>=0;i--) //printing answer
    printf("%d", a[i]);
    printf("\n");
    t ++;
    }

}

int execute_fatorial (int n) {

    clock_t tempo_inicial;
    clock_t tempo_final;

    tempo_inicial = clock();
    fatorial(n);
    tempo_final = clock();

    printf("\n################################################################\n");
    printf("# O tempo que o C levou para o calculo foi de:                 #\n");
    printf("# %e                                                 #\n",((double)(tempo_final - tempo_inicial) / CLOCKS_PER_SEC));
    printf("################################################################\n");

    return 0;
}
