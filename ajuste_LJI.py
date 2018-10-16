################################################################################
#
# Autor: Sandro F. de Brito
#
# Descrição:
#    Esse script tem como característica ajustar uma curva a partir
# de um conjunto de pontos de um arquivo, os quais foram calculas previamente
# por métodos computacionais.
#    O script é um teste, pois gera os pontos a partir da função de interesse,
# que neste caso e a função Lennard Jones Improve. Buscamos ajustar os valores:
#
#     - de (enerdia de dissociação)
#     - req (distância de equilíbrio).
#
#    Os resultados são mostrados em uma curva, bem como em um arquivo de saida,
# denominado de informações_do_ajuste.txt.
#
#    No caso deste scrip, usamos valores experiemtais conhecidos para de e req.
# Esses valores são da interação entre H2O - Xe, ou seja, água xenônio.
#    Os valores experiemtais são:
#
#    de = 20,20 meV  | req = 3,93 angstron
#
# -----> O ajuste buscará encontrar valores para de e req. <------
#
################################################################################


from scipy.optimize import curve_fit
import numpy as np
import matplotlib.pyplot as plt
from lmfit import Model
import seaborn as sns

# A função abaixo define o potencial que estamos trabalhando:
# Lennard - Jones Improve.
def LJ_Inprove(r, de, req):
    n_r = (9 + 4*((r/req)**2))
    return ((de/(n_r - 6)) * ((6 * ((req/r)**(n_r))) - (n_r * ((req/r)**6))))


r = np.linspace(3.4, 7, 20) # Distância de Interesse de 3,4 até 7,0.
de = 20.20                  # Esse valor de é um valor experimental
req = 3.93                  # Esse valor de req também é um valor experiemtnal

# A função y abaixo representa os dados que deverão ser ajustados.
# Esse tipo de função é chamada de função ruido (noise)
y = LJ_Inprove(r, de, req) + np.random.normal(0, 0.3, len(r))


chut_inicial = [1, 1] # chute inicial do ajuste com scipy
melhor_variavel, covar = curve_fit(LJ_Inprove, r, y, chut_inicial)

print('---------------------------------------')
print('  Valores de Ajuste Scipy. curve_fit \n')
print('# de = {}'.format(melhor_variavel[0]))
print('# requilibrio = {}'.format(melhor_variavel[1]))
print('---------------------------------------\n')

# Armazenando os valores a partir da função ruido em no arquivo dados.dat.
arq = open('dados.dat', 'w')
for distancia in r:
    arq.write('{0:8.3f} {1:35.25f} \n'.format(distancia,
               LJ_Inprove(distancia, melhor_variavel[0], melhor_variavel[1])))
arq.close()

# Lendo os valores da função ruido afim de obter um ajuste um ajuste pelo método
# fit do módulo lmfitself.
dados = np.loadtxt('dados.dat')
x1 = dados[:, 0]
y1 = dados[:, 1]

novo_ajuste = Model(LJ_Inprove) # definimos o modelo de ajuste com Model
resultado_ajuste = novo_ajuste.fit(y1, r=x1, de=1, req=1) # ajustando

d = open('informações_do_ajuste.txt', 'w')
d.write('--------------------------------------------------\n')
d.write('- Informações sobre o Ajuste realizado com Scipy -\n')
d.write('--------------------------------------------------\n\n')
d.write('# de = {} \n'.format(melhor_variavel[0]))
d.write('# requilibrio = {} \n'.format(melhor_variavel[1]))
d.write('-------------------------------------------------- \n\n')
d.write('--------------------------------------------------\n')
d.write('- Informações sobre o Ajuste realizado com LmFit -\n')
d.write('--------------------------------------------------\n')
d.write(resultado_ajuste.fit_report())
d.write('--------------------------------------------------\n')
d.close()

print('---------------------------------------------------')
print('      Valores de Ajuste LmFit. Fit Report \n')
print(resultado_ajuste.fit_report())
print('-----------------------------------------------------')

sns.set()
sns.set_style("white")
fig, ax = plt.subplots()
ax.plot(r, y, 'o', label='Pontos do Ajuste') # y é a função que gera os pontos a serem ajustados
ax.plot(r, LJ_Inprove(r, melhor_variavel[0], melhor_variavel[1]), label='Ajuste com Scipy', linewidth=2.5, color='red')
ax.plot(x1, resultado_ajuste.best_fit, '--', label='Ajuste com LmFit', linewidth=2.5, color='dimgray')
ax.set_xlabel(r'R (Å)')
ax.set_ylabel(r'Energy ($cm^{-1}$)')
ax.legend(loc='best')
plt.title(r'Ajuste $H_{2}O - Xe$', fontsize=16)
plt.show()
