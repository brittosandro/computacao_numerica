import matplotlib.pyplot as plt
import numpy as np
from sympy import *
from sympy.stats import *

'''
#Ep = Function('Ep')
#Ep = symbols('Ep')
#kb, T = symbols('kb T', constant = True)
#pprint(series(exp(-Ep)), use_unicode=True)

x, n = symbols('x n', integer=True)
n = Symbol('n', integer=True)    # importante deixar claro que deve ser inteiro
somas = series(exp(-x))
exp_expr = ((-1)**n) * (x**n / factorial(n))   # verifique que essa é realmente a expressão
sum_exp = Sum(exp_expr, (n, 0, oo))

print(somas)
pprint(somas)
print('\n\n')

print(sum_exp)
pprint(sum_exp)

print('\n\n')
print(sum_exp.doit())
pprint(sum_exp.doit())
'''

def potencial_LJ(r, epsilon, sigma):
    '''
    Essa função recebe valores de distâncias entre dois átomos
    e retorna o valor da energia de interação.

    Entrada (r/angstron, epsilon/eV, sigma/angstron)
    -------
    A unidade de medida da distância é em angstron, epsilon é eV e
    sigma é angstrom.

    Retorno (E/milieletron-volt)
    -------
    A unidade de saida/retorno é em meV.
    '''
    e = epsilon
    s = sigma
    E = 4*e*(((s/r)**12) - ((s/r)**6))
    return E * (10**(3))

def improve_LJ(r, de, req, beta):
    '''
    Essa função recebe valores de distancias entre dois átomos,
    o valor de profundidade do poço, a distância de equilibrio
    e o valor de beta, retornando a energia de interação

    Entrada (req/angstron, de/mev, )
    -------
    A unidade da distância é angstron e do poço é em meV.

    Retorno (E/milieletron-volt)
    -------
    A unidade de saida/retorno é em meV.
    '''

    n_r = (beta + 4*((r/req)**2))
    E = ((de/(n_r - 6)) * ((6 * ((req/r)**(n_r))) - (n_r * ((req/r)**6))))
    return E

def integra_potencial(r_min, r_max):
    r = symbols('r')
    return integrate(pot, (r, r_min, r_max))

def coef_virial(T, integral):
    Na = 6.022140e23
    kb = 8.617333262e-5
    kbT = kb * T
    a = (2*np.pi)/kbT
    B = a * Na * integral
    return B * 1.0e-30

def coef_virial1(integral):
    '''
    Essa função calcula o coeficiente do virial considerando a aproximação
    em que r é pequeno. Portanto r é aproximadamente zero, enquanto o potencial
    de interaçõe é um valor positivo e grande. A função recebe o valor de
    uma integral (pm^3) e retorna o coeficiente do virial em cm^3/mol.

    Entrada (integral/pm^3)
    -------

    Retorno (coeficiente do Virial/(cm^3/mol))
    -------
    '''
    Na = 6.022140e23
    a = 2 * np.pi
    B = a * Na * integral
    return B * 1.0e-30

if __name__ == "__main__":
    distancias = np.arange(3.10, 8.20, 0.01)

    #Parâmetros para o potencial de Lennard-Jones Ar2
    #ref:
    epsilon = 0.0103332 #eV
    sigma = 3.40        #angstrom
    energias = [potencial_LJ(r, epsilon, sigma) for r in distancias]

    #Parâmetros para o potencial de Lennard-Jones Improve Ar2
    #ref: A Spectroscopic Validation of the Improved Lennard–Jones Model
    de = 12.343  #meV
    req = 3.76   #angstrom
    beta = 9.74
    energias1 = [improve_LJ(r, de, req, beta) for r in distancias]

    #Parâmetros para o potencial de Lennard-Jones Ar2.
    #-------------------------- referência -------------------------------------
    #nome: Determining Intermolecular Potentials from second virial coefficients
    #autor: Brian P. Reid
    #---------------------------------------------------------------------------

    epsilon1 = 140 #kK lembre k = 1,38064x10^23 j/K -> kK = 1,38064x10^23 j
    sigma1 = 335   #pm (picometros)
    distancias1 = np.arange(305.21, 1303.91)
    energias2 = [potencial_LJ(r, epsilon1, sigma1)*(1/(1000*epsilon1))
                                                  for r in distancias1]
    minima = 0
    for i in range(1, len(energias2)):
        if energias2[minima] > energias2[i]:
            minima = i
            j = i

    #print('*****************')
    #print(minima)
    #print(energias2[minima])
    r_min1 = distancias1[j]
    #print(r_min1)
    #print('*****************')

    new_distancias1 = [(r/376.21) for r in distancias1]

    plt.plot(new_distancias1, energias2, color='b', linewidth=2.5, label='LJ')
    plt.legend(loc='upper right', shadow=False, fontsize='large',
               bbox_to_anchor=(0.98, 0.98), frameon=False)
    plt.ylabel(r'Energia ($U=U/\epsilon$)')
    plt.xlabel(r'R ($r/r_{min}$)')
    plt.title(r'Interação $Ar_{2}$')
    plt.show()

    # Para calcular o coeficiente do virial total vamos repartir o coeficiente
    #em quatro partes de acordo com o a distância de interação.
    # B(T) = B1(T) + B2(T) + B3(T) + B4(T)
    #
    # Cálculo de B1(T)
    # Região em que r é pequeno
    # -------------------------
    # Se r varia de 0 até r1, então podemos calcular B1(T) assumindo que a
    #integral que devemos resolver é \int_{0}^{r1}r^{2}dr.

    r = symbols('r')
    r1 = distancias1[0]
    #print(r1)
    #print(f'r inicial = {r1} | rmin = {r_min1}')
    int1 = integrate(r**2, (r, 0, r1))
    #print(f'Integral = {int1}')
    B1 = coef_virial1(int1)
    print(f'B1(T) = {B1} cm³/mol')

    # Cálculo de B2(T)
    # Região em que r1 <= r <= r2
    # ---------------------------







    '''
    minima = 0
    for i in range(1, len(energias1)):
        if energias1[minima] > energias1[i]:
            minima = i
            j = i

    r_min = distancias[j]
    minimo1 = min(energias1)
    print('--------')
    print(minimo1)
    print(energias1[minima])
    print(r_min)
    print('--------')

    plt.plot(distancias, energias, color='b', linewidth=2.5, label='LJ')
    plt.plot(distancias, energias1, color='r', linewidth=2.5, label='ILJ')
    plt.legend(loc='upper right', shadow=False, fontsize='large',
               bbox_to_anchor=(0.98, 0.98), frameon=False)
    plt.ylabel(r'Energia (meV)')
    plt.xlabel(r'R ($\AA$)')
    plt.title(r'Interação $Ar_{2}$')
    plt.show()
    '''
    #r = symbols('r')
    #r1 = distancias[0]
    #print(f'r inicial = {r1} | rmin = {r_min}')
    #int = integrate(r**2 * potencial_LJ(r, epsilon, sigma), (r, r1, r_min))
    #print(int)
    #temperaturas = np.arange(100, 1200, 25)
    #B = [coef_virial(T, int) for T in temperaturas]
    #print(B)
    #plt.plot(temperaturas, B)
    #plt.show()
