#############################################################################
#
# Os gráficos listados não apresentam medidas em seus eixos.
#
#############################################################################

import matplotlib.pyplot as plt

fig = plt.figure(figsize=(8, 5))     # Criando um objeto da classe figure

# --------------- Primeiro Subplot ------------------------------------------
subplot1 = fig.add_subplot(2, 2, 1)     # subplot(n_linhas, n_colunas, n_plots)
subplot1.set_xticks([])
subplot1.set_yticks([])
subplot1.text(
               0.5,                              # Ordena posição x
               0.5,                              # Ordena posição y
               'Subplot 1 (2, 2, 1)',            # Texto dentro do quadrante
               horizontalalignment='center',     # Alinha Texto
               verticalalignment='center',
               fontsize=20,                      # Tamanho da fonte
               color='b'
             )
#----------------------------------------------------------------------------

# -------------- Segundo Subplot --------------------------------------------
subplot2 = fig.add_subplot(2, 2, 4)    # subplot(n_linhas, n_colunas, n_plots)
subplot2.set_xticks([])
subplot2.set_yticks([])
plt.text(
          0.5,
          0.5,
          'Subplot 2 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20, color='r'
        )
# ---------------------------------------------------------------------------

plt.savefig('Mult_fig_4.png', dpi=300, orientation='portrait', transparent=True, format='png')
plt.show()
