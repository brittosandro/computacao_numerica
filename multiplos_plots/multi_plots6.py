import matplotlib.pyplot as plt

plt.figure(figsize=(16, 24))     #figsize=(coluna_x, coluna_y)

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

# -------------- Segundo Subplot --------------------------------------------
plt.subplot(4, 3, 2)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 2 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20, color='r'
        )
# ---------------------------------------------------------------------------

# -------------- Terceiro Subplot --------------------------------------------
plt.subplot(4, 3, 3)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 3 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#fcb001'
        )
# ---------------------------------------------------------------------------

# -------------- Quarto Subplot --------------------------------------------
plt.subplot(4, 3, 4)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 4 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#fa5ff7'
        )
# ---------------------------------------------------------------------------

# -------------- Quinto Subplot --------------------------------------------
plt.subplot(4, 3, 5)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 5 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#01a049'
        )
# ---------------------------------------------------------------------------

# -------------- Sexto Subplot --------------------------------------------
plt.subplot(4, 3, 6)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 6 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#ff474c'
        )
# ---------------------------------------------------------------------------

# -------------- Sétimo Subplot --------------------------------------------
plt.subplot(4, 3, 7)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 7 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#acbf69'
        )
# ---------------------------------------------------------------------------

# -------------- Oitavo Subplot --------------------------------------------
plt.subplot(4, 3, 8)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 8 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#7e1e9c'
        )
# ---------------------------------------------------------------------------

# -------------- Nono Subplot --------------------------------------------
plt.subplot(4, 3, 9)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 9 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#5b7c99'
        )
# ---------------------------------------------------------------------------

# -------------- Décimo Subplot --------------------------------------------
plt.subplot(4, 3, 10)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 10 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#789b73'
        )
# ---------------------------------------------------------------------------

# -------------- Décimo Primeiro Subplot ------------------------------------
plt.subplot(4, 3, 11)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 11 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#8c000f'
        )
# ---------------------------------------------------------------------------

# -------------- Décimo segundo Subplot ------------------------------------
plt.subplot(4, 3, 12)                        # subplot(n_linhas, n_colunas, n_plots)
plt.text(
          0.5,
          0.5,
          'Subplot 12 (2, 2, 4)',
          horizontalalignment='center',
          verticalalignment='center',
          fontsize=20,
          color='#580f41'
        )
# ---------------------------------------------------------------------------

plt.savefig('Mult_fig_5.png', dpi=300, orientation='portrait', transparent=True, format='png')
plt.show()
