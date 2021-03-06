import matplotlib.pyplot as plt

plt.figure(figsize=(8, 5))     #figsize=(coluna_x, coluna_y)

# --------------- Primeiro Subplot ------------------------------------------
plt.subplot(2, 2, 1)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,                              # Ordena posição x
          0.5,                              # Ordena posição y
          'Subplot 1 (2, 2, 1)',            # Texto dentro do quadrante
          horizontalalignment='center',     # Alinha Texto
          verticalalignment='center',
          fontsize=20,                      # Tamanho da fonte
          color='b'                         # Define cor do texto
        )
#-----------------------------------------------------------------------------

# -------------- Segundo Subplot ---------------------------------------------
plt.subplot(2, 2, 2)
plt.text(
         0.5,
         0.5,
         'Subplot 2 (2, 2, 2)',
         horizontalalignment='center',
         verticalalignment='center',
         fontsize=20,
         color='k'
       )
# ----------------------------------------------------------------------------

# -------------- Terceiro Subplot --------------------------------------------
plt.subplot(2, 2, 3)                      # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 3 (2, 2, 3)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20, color='g'
        )
# ---------------------------------------------------------------------------

# -------------- Quarto Subplot --------------------------------------------
plt.subplot(2, 2, 4)                      # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 4 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20, color='r'
        )
# ---------------------------------------------------------------------------

plt.show()
