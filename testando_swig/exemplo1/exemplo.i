/* exemplo.i */
%module exemplo

%{
/* Coloque arquivos de cabeçalho aqui! ou declarações de funções como as abaixo */
#define SWING_FILE_WITH_INIT
#include "exemplo.h"
%}

unsigned long fact (unsigned long n);
void exec_fact (int n);
