############################################################################
#    Autor: Sandro F. Brito
#
#    objetivo: Testar o uso do parâmetro  [-- constrained_layout --].
#
#    Descrição:
#
#    Este programa tem por finalidade fazer um teste com a variável
# constrained_layout, essa é uma variável booleana que neste programa
# será considerada com atributo True.
#    O parâmetro [-- constrained_layout --] tem como característica
# ordenar a posição dos subtítulos do gráfico para que estes não fiquem
# sobrepostos nos diferentes subplots.
#    Neste caso iremos ver se a figura torna-se ordenada.
#
#############################################################################

import matplotlib.pyplot as plt
import numpy as np

def oscilador_amortecido(t):
    a = np.cos(2*np.pi*t)
    b = np.exp(-t)
    return a*b

t1 = np.arange(0.0, 5.0, 0.1)
t2 = np.arange(0.0, 5.0, 0.02)
t3 = np.arange(0.0, 2.0, 0.01)

fig, axs = plt.subplots(2, 1, constrained_layout=True)
fig.suptitle(
              'Testando o uso do Parâmetro constrained_layout',
              fontweight = 'bold',
              fontsize=15
            )

# ------------ Subplot 1 --------------------------------------------------
axs[0].plot(t1, oscilador_amortecido(t1), 'o', color='r')
axs[0].plot(t2, oscilador_amortecido(t2), '-', lw=2, color='b')
axs[0].set_title('Subplot 1', fontsize = 'large', fontweight = 'bold')
axs[0].set_xlabel('tempo (s)', fontsize = 'large')
axs[0].set_ylabel('Oscilador Amortecido', fontsize = 'large')

# ---------- Subplot 2 -----------------------------------------------------
axs[1].plot(t3, np.cos(2*np.pi*t3), '-', lw=2, color='b')
axs[1].set_title('Subplot 2', fontsize = 'large', fontweight = 'bold')
axs[1].set_xlabel('tempo (s)', fontsize = 'large')
axs[1].set_ylabel('Oscilador', fontsize = 'large')

plt.savefig('osciladores1.png', dpi=300, orientation='portrait', transparent=True, format='png')
plt.show()
