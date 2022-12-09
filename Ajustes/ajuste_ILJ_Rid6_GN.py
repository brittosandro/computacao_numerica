################################################################################
#
# Descrição:
#    Esse script tem como característica ajustar uma curva a partir
# de um conjunto de pontos de um arquivo, os quais foram calculas previamente
# por métodos computacionais.
#    O script usa dois modelos de função:
# 1) A função Ridberg de grau 6 para obter os valores das constantes:
#    - c1, c2, c3, c4, c5, c6
#    Bem como os valores de enerdia de dissociação e R de equilíbrio:
#     - de (enerdia de dissociação)
#     - req (distância de equilíbrio).
# 2) A funlção Improved Lenard-Jones (IJL) para obter os valores da constante:
#    - beta
#    Bem como os valores de energia de dissociação e R de equilibrío:
#     - de (enerdia de dissociação)
#     - req (distância de equilíbrio).
################################################################################


from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import seaborn as sns


# A função abaixo define o potencial Lennard - Jones Improve.
def LJ_Inprove(r, de, req, beta):
    n_r = (beta + 4*((r/req)**2))
    return ((de/(n_r - 6)) * ((6 * ((req/r)**(n_r))) - (n_r * ((req/r)**6))))

def LJ(r, epsilon, sigma):
    e = epsilon
    s = sigma
    return 4 * e * ((s/r)**12 - (s/r)**6)

# Definindo Ryd 6 potential
def rydberg6(x, de, xe, c1, c2, c3, c4, c5, c6):
    potential = -de*(1 + c1*(x-xe) + c2*(x-xe)**2 + c3*(x-xe)**3
    + c4*(x-xe)**4 + c5*(x-xe)**5 + c6*(x-xe)**6)*np.exp(-c1*(x-xe))
    return potential


# Lendo os valores da função ruido afim de obter um ajuste um ajuste pelo método
# fit do módulo lmfitself.
dados = np.loadtxt('medias_Kr_amonia.dat')
x1 = dados[:, 0]
y1 = dados[:, 1]

chut_inicial = [-15.88, 4.01, 8.60,] # chute inicial do ajuste com scipy
melhor_variavel, covar = curve_fit(LJ_Inprove, x1, y1, chut_inicial, maxfev=200000)

print('---------------------------------------')
print('  Valores de Ajuste Scipy. curve_fit \n')
print('# de = {}'.format(melhor_variavel[0]))
print('# req = {}'.format(melhor_variavel[1]))
print('# BETA = {}'.format(melhor_variavel[2]))
print('---------------------------------------\n')

novo_ajuste = Model(LJ_Inprove) # definimos o modelo de ajuste com Model
resultado_ajuste = novo_ajuste.fit(y1, r=x1, de=-15.816, req=4.005, beta=8.60)

print('---------------------------------------------------')
print('      Valores de Ajuste LmFit. Fit Report \n')
print(resultado_ajuste.fit_report(modelpars=None, show_correl=True))
print('-----------------------------------------------------')


d = open('info_ajuste_CEP_Kr.txt', 'w')
d.write('--------------------------------------------------\n')
d.write('- Informações sobre o Ajuste realizado com Scipy -\n')
d.write('- A unidade de medida para De é dada em MEV     -\n')
d.write('--------------------------------------------------\n\n')
d.write('# de = {} \n'.format(melhor_variavel[0]))
d.write('# req = {} \n'.format(melhor_variavel[1]))
d.write('# BETA = {} \n'.format(melhor_variavel[2]))
d.write('-------------------------------------------------- \n\n')

d.write('--------------------------------------------------\n')
d.write('- Informações sobre o Ajuste realizado com LmFit -\n')
d.write('- A unidade de medida para De é dada em kcal/mol -\n')
d.write('--------------------------------------------------\n')
d.write(resultado_ajuste.fit_report())
d.write('--------------------------------------------------\n')
d.close()

sns.set(style="ticks")
fig, ax = plt.subplots()

ax.plot(x1, y1, 'o', color='#0000ff', label='SAPT2+/aug-cc-pvtz', linewidth=2.)
ax.plot(x1, resultado_ajuste.best_fit, '-',  color='#e80c3d' ,label='ILJ', linewidth=2.)
ax.set_xlabel('Distância ($\AA$)', labelpad = -1.0, fontsize = 'xx-large')
ax.set_ylabel('Energia (meV)', labelpad = 0.0, fontsize = 'xx-large')
#ax.set_xlim(3., 9.3)
#ax.set_ylim(-17.5, 11.5)
ax.legend(
           loc='lower right',
           ncol = 1,
           fontsize='large',
           bbox_to_anchor=(0.999, 0.75),
           fancybox = True,
           framealpha=0.1
           )
ax.set_title('Kr + NH3', fontsize = 'xx-large', fontweight = 'bold')
plt.savefig('ajuste_Kr_ILJ.png', dpi=300, orientation='portrait', transparent=True, format='png')
plt.show()
