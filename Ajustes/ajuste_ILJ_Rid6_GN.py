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
# 3) A terceira função definida é a Lennard Jones clássica. Os parâmetros que
#   estamos buscando ajustar são:
#     - epsilon.
#     - sigma.
################################################################################


from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import seaborn as sns


# A função abaixo define o potencial Lennard - Jones Improve.
def ILJ(r, de, req, beta):
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

def main():
    # Lendo os dados para realizar o ajuste.
    dados = np.loadtxt('medias_Kr_amonia.dat')
    x1 = dados[:, 0]
    y1 = dados[:, 1]

    # chute inicial do ajuste com scipy para ILJ
    chut_inicial = [-19.88, 3.80, 8.60,] # chute inicial do ajuste com scipy
    melhor_variavel, covar = curve_fit(ILJ, x1, y1, chut_inicial, maxfev=200000)

    print('---------------------------------------')
    print('  Valores de Ajuste Scipy. curve_fit \n')
    print('# de = {}'.format(melhor_variavel[0]))
    print('# req = {}'.format(melhor_variavel[1]))
    print('# BETA = {}'.format(melhor_variavel[2]))
    print('---------------------------------------\n')

    # definimos o modelo de ajuste com Model
    novo_ajuste = Model(ILJ)
    resultado_ajuste = novo_ajuste.fit(y1, r=x1, de=-19.816, req=3.80, beta=8.60)

    print('---------------------------------------------------')
    print('      Valores de Ajuste LmFit ILJ. Fit Report \n')
    print(resultado_ajuste.fit_report(modelpars=None, show_correl=True))
    print('-----------------------------------------------------')

    d = open('info_ajuste_CEP_ILJ.txt', 'w')
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
    d.write('- A unidade de medida para De é dada em mev      -\n')
    d.write('--------------------------------------------------\n')
    d.write(resultado_ajuste.fit_report())
    d.write('--------------------------------------------------\n')
    d.close()

    sns.set(style="ticks")
    fig, ax = plt.subplots()

    # Ajusta subplots.
    fig.subplots_adjust(
                         left = 0.125,
                         right = 0.90,    # {Define as distâncias entre os extremos}
                         bottom = 0.135,
                         top = 0.900,
                         hspace = 0.200,   # Organiza espaçoes entre os subplots
                         wspace = 0.200    # Organiza espaçoes entre os subplots
                        )

    ax.plot(x1, y1, 'o', color='#0000ff', label='SAPT2+/aug-cc-pvtz', linewidth=2.)
    ax.plot(x1, resultado_ajuste.best_fit, '-',  color='#e80c3d', label='ILJ', linewidth=2.)
    ax.set_xlabel('Distância ($\AA$)', labelpad = -1.0, fontsize = 'xx-large')
    ax.set_ylabel('Energia (meV)', labelpad = 0.0, fontsize = 'xx-large')
    ax.legend(
               loc='lower right',
               ncol = 1,
               fontsize='large',
               bbox_to_anchor=(0.999, 0.75),
               fancybox = True,
               framealpha=0.1
              )
    ax.set_title('$Kr + NH_{3}$ $(Sítio$ $1)$', fontsize = 'xx-large', fontweight = 'bold')
    plt.savefig('ajuste_ILJ.png', dpi=500, orientation='portrait', transparent=True, format='png')
    plt.savefig('ajuste_ILJ.pdf', dpi=500, orientation='portrait', transparent=True, format='pdf')


    # chute inicial do ajuste com scipy para IL
    chut_inicial = [-15.88, 4.01]
    melhor_variavel, covar = curve_fit(LJ, x1, y1, chut_inicial, maxfev=200000)

    print('---------------------------------------')
    print('Valores de Ajuste Scipy LJ. curve_fit \n')
    print('# epsilon = {}'.format(melhor_variavel[0]))
    print('# sigma = {}'.format(melhor_variavel[1]))
    print('---------------------------------------\n')

    # definimos o modelo de ajuste com Model
    novo_ajuste = Model(LJ)
    resultado_ajuste = novo_ajuste.fit(y1, r=x1, epsilon=-15.816, sigma=4.005)

    print('---------------------------------------------------')
    print('      Valores de Ajuste LmFit LJ. Fit Report \n')
    print(resultado_ajuste.fit_report(modelpars=None, show_correl=True))
    print('-----------------------------------------------------')

    d = open('info_ajuste_CEP_LJ.txt', 'w')
    d.write('--------------------------------------------------\n')
    d.write('- Informações sobre o Ajuste realizado com Scipy  -\n')
    d.write('- A unidade de medida para De é dada em MEV       -\n')
    d.write('--------------------------------------------------\n\n')
    d.write('# de = {} \n'.format(melhor_variavel[0]))
    d.write('# req = {} \n'.format(melhor_variavel[1]))
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
    ax.plot(x1, resultado_ajuste.best_fit, '-',  color='#e80c3d', label='LJ', linewidth=2.)
    ax.set_xlabel('Distância ($\AA$)', labelpad = -1.0, fontsize = 'xx-large')
    ax.set_ylabel('Energia (meV)', labelpad = 0.0, fontsize = 'xx-large')
    ax.legend(
               loc='lower right',
               ncol = 1,
               fontsize='large',
               bbox_to_anchor=(0.999, 0.75),
               fancybox = True,
               framealpha=0.1
              )
    ax.set_title('$Kr + NH_{3}$ $(Sítio$ $1)$', fontsize = 'xx-large', fontweight = 'bold')
    plt.savefig('ajuste_LJ.png', dpi=500, orientation='portrait', transparent=True, format='png')
    plt.savefig('ajuste_LJ.pdf', dpi=500, orientation='portrait', transparent=True, format='pdf')

    # chute inicial do ajuste com Rydberg 6
    chut_inicial = [-19.28, 3.80, 2, -1, 1, 1, 1, 1]
    melhor_variavel, covar = curve_fit(rydberg6, x1, y1, chut_inicial, maxfev=200000)

    print('---------------------------------------')
    print('  Valores de Ajuste Scipy. curve_fit \n')
    print('# de = {}'.format(melhor_variavel[0]))
    print('# req = {}'.format(melhor_variavel[1]))
    print('# c1 = {}'.format(melhor_variavel[2]))
    print('# c2 = {}'.format(melhor_variavel[3]))
    print('# c3 = {}'.format(melhor_variavel[4]))
    print('# c4 = {}'.format(melhor_variavel[5]))
    print('# c5 = {}'.format(melhor_variavel[6]))
    print('# c6 = {}'.format(melhor_variavel[7]))
    print('---------------------------------------\n')

    # definimos o modelo de ajuste com Model
    novo_ajuste = Model(rydberg6)
    resultado_ajuste = novo_ajuste.fit(y1, x=x1, de=-19.28, xe=3.80, c1=2, c2=-1,
                                   c3=1, c4=1, c5=1, c6=1)

    print('---------------------------------------------------')
    print('      Valores de Ajuste LmFit. Fit Report \n')
    print(resultado_ajuste.fit_report(modelpars=None, show_correl=True))
    print('-----------------------------------------------------')

    d = open('info_ajuste_CEP_Rydberg6.txt', 'w')
    d.write('--------------------------------------------------\n')
    d.write('- Informações sobre o Ajuste realizado com Scipy -\n')
    d.write('- A unidade de medida para De é dada em cm-1     -\n')
    d.write('--------------------------------------------------\n\n')
    d.write('# de = {} \n'.format(melhor_variavel[0]))
    d.write('# req = {} \n'.format(melhor_variavel[1]))
    d.write('# c1 = {} \n'.format(melhor_variavel[2]))
    d.write('# c2 = {} \n'.format(melhor_variavel[3]))
    d.write('# c3 = {} \n'.format(melhor_variavel[4]))
    d.write('# c4 = {} \n'.format(melhor_variavel[5]))
    d.write('# c5 = {} \n'.format(melhor_variavel[6]))
    d.write('# c6 = {} \n'.format(melhor_variavel[7]))
    d.write('-------------------------------------------------- \n\n')

    d.write('--------------------------------------------------\n')
    d.write('- Informações sobre o Ajuste realizado com LmFit -\n')
    d.write('- A unidade de medida para De é dada em mev -\n')
    d.write('--------------------------------------------------\n')
    d.write(resultado_ajuste.fit_report())
    d.write('--------------------------------------------------\n')
    d.close()

    sns.set(style="ticks")
    fig, ax = plt.subplots()

    # Ajusta subplots.
    fig.subplots_adjust(
                         left = 0.130,
                         right = 0.930,    # {Define as distâncias entre os extremos}
                         bottom = 0.140,
                          top = 0.905,
                          hspace = 0.200,   # Organiza espaçoes entre os subplots
                          wspace = 0.200    # Organiza espaçoes entre os subplots
                        )

    ax.plot(x1, y1, 'o', color='#0000ff', label='SAPT0/aug-cc-PVQZ', linewidth=2.)
    ax.plot(x1, resultado_ajuste.best_fit, '-',  color='#e80c3d' ,label='Rydberg 6', linewidth=2.)
    ax.set_xlabel('Distância ($\AA$)', labelpad = -1.5, fontsize = 'xx-large')
    ax.set_ylabel('Energia (meV)', labelpad = 0.0, fontsize = 'xx-large')
    ax.legend(
               loc='lower right',
               ncol = 1,
               fontsize='large',
               bbox_to_anchor=(0.999, 0.75),
               fancybox = True,
               frameon=False
              )
    ax.set_title('$Kr + NH_{3}$ $(Sítio$ $1)$', fontsize = 'xx-large', fontweight = 'bold')
    plt.savefig('ajuste_Ryd6.png', dpi=500, orientation='portrait', transparent=True, format='png')
    plt.savefig('ajuste_Ryd6.pdf', dpi=500, orientation='portrait', transparent=True, format='pdf')


if __name__ == '__main__':
    main()
