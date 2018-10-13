import matplotlib.pyplot as plt

def gera_multiplos_plots(cor='r'):

    '''
        A função gera_multiplos_plots recebe um argumento opcional que irá
        descrever a cor do texto dentro da caixa e devolve um gráfico com
        múltiplas figuras, neste caso com 12 figuras.
    '''

    plt.figure(figsize=(16, 10))     #figsize=(coluna_x, coluna_y)

    for i in range(12):
        plt.subplot(4, 3, i+1)                                     # subplot(n_linhas, n_colunas, n_plots)
        plt.text(
                  0.5,                                             # Ordena posição x
                  0.5,                                             # Ordena posição y
                  'Subplot {} (2, 2, {})'.format(i+1, i+1),        # Texto dentro do quadrante
                  horizontalalignment='center',                    # Alinha Texto
                  verticalalignment='center',
                  fontsize=20,                                     # Tamanho da fonte
                  color = cor                                      # Cor da fonte
                )

    plt.savefig('Mult_fig_6.png', dpi=300, orientation='portrait', transparent=True, format='png')
    plt.show()

if __name__ == '__main__':

    #cores = ['#0b5509', '#929901', '#01f9c6', '#017a79', '#2000b1']
    #[gera_multiplos_plots(cor) for cor in cores]
    gera_multiplos_plots()
