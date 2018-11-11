/* novo_fatorial.i */
%module fatorial

%{
/* Coloque arquivos de cabeçalho aqui! ou declarações de funções como as abaixo */
#define SWING_FILE_WITH_INIT
#include "fatorial.h"
%}

void fatorial(int n);
int execute_fatorial (int n);
