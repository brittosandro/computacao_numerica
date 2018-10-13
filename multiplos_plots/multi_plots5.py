import matplotlib.pyplot as plt

plt.figure(figsize=(12, 16))     #figsize=(coluna_x, coluna_y)

# --------------- Primeiro Subplot ------------------------------------------
plt.subplot(4, 3, 1)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,                              # Ordena posição x
          0.5,                              # Ordena posição y
          'Subplot 1 (2, 2, 1)',            # Texto dentro do quadrante
          horizontalalignment='center',     # Alinha Texto
          verticalalignment='center',
          fontsize=20,                      # Tamanho da fonte
          color='b'
        )
#----------------------------------------------------------------------------

# -------------- Terceiro Subplot --------------------------------------------
plt.subplot(4, 3, 12)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 2 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20, color='r'
        )
# ---------------------------------------------------------------------------

plt.savefig('Mult_fig_5.png', dpi=300, orientation='portrait', transparent=True, format='png')
plt.show()
