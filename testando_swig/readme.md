<h2> Informações Sobre uso do SWIG </h2>

<p>
    Para utilizar o compilador de interface <em>SWIG</em> deveremos </br>
    construir três módulos ou arquivos. Os quais são discriminados  </br>
    da seguinte maneira:</br></br>  
</p>

<p>
   1) O primeiro módulo que deveremos criar é um arquivo.i; </br>

      - Neste módulo utilizamos um cabeçalho que usa um conjunto de diretivas buscando</br>
      empacotar o programa, a partir desse arquivo criamos dois outro módulos: </br>

	      - arq1.py 
	      - arq1_wrap.c 

   2) O segundo módulo que deveremos criar é um arquivo de cabeçalho arquivo.h; </br>
      
   3) O terceiro módulo que criaremos é o programa em C. </br>               
</p>

<p>
   A execução do programa deve seguir os seguintes passos: </br>
   
   1) swig  -python  exemplo.i
   2) gcc  -O2  -fPIC  -c  exemplo.c
   3) gcc  -O2  -fPIC  -c  exemplo_wrap.c  -I/usr/include/python3.6
   4) gcc  -shared  exemplo.o  exemplo_wrap.o  -o  _exemplo.so
   
</p>
