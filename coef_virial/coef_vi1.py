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

def calc_pontos_pi(dist1, dist2, x):
    p = (1/2)*(dist1 + dist2) + (1/2)*(dist2 - dist1)*x
    return p

def quad_gaussian(pi, xi, wi, T):
    '''
    A quadratura gaussiana retorna o valor valor em cm³.
    '''

    kb = 1         #j/K (joule por kelvin)
    kbT = kb * T
    ep = 140       #kK (lembre k = 1,38064x10^-23 j/K -> kK = 1,38064x10^-23 j)
    si = 335       #pm (picometros)
    potencial = potencial_LJ(pi, ep, si)*(1/1000)
    # Aqui o potencial esta em kK.
    #print(potencial)
    a = -potencial/kbT
    #print(a)
    c = wi*(1 - np.exp(a))*(pi**2)
    return c * 1e-30


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

    # Parâmetros para o potencial de Lennard-Jones Ar2.
    # -------------------------- referência -------------------------------------
    # nome: Determining Intermolecular Potentials from second virial coefficients
    # autor: Brian P. Reid
    # ---------------------------------------------------------------------------

    epsilon1 = 140 #kK (lembre k = 1,38064x10^-23 j/K -> kK = 1,38064x10^-23 j)
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
    #plt.show()

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
    #  Nesse caso nos vamos resolver a integral:
    # \int_{r1}^{r2}(1-e^(U(r)/kt))r^{2}dr, pelo método de aproximação conhecido
    # como quadratura gaussiana. Nesse caso vamos considerar n = 6, ou seja,
    # seis termos para a aproximação. Teremos então:
    # \int_{r1}^{r2}(1-e^(U(r)/kt))r^{2}dr ~ \sum_{i=1}^{6}w_i(1-e^(U(p_i)/kt))p_i^{2}.
    #  Em que p_i = 1/2*(r1+r2) + 1/2*(r2-r1)xi.
    #  Os valores de x_i e w_i são tabelados e foram retirados do sítio:
    # https://pomax.github.io/bezierinfo/legendre-gauss.html

    coef_wi_xi_6 = np.loadtxt('dados_coef6.dat', comments='#')
    wi_6 = coef_wi_xi_6[:, 1]
    xi_6 = coef_wi_xi_6[:, 2]

    # Pegando o ponto associado a energia mínima do potencial.
    r2 = r_min1
    #print(r1, r2)
    p_i = [calc_pontos_pi(r1, r2, xi) for xi in xi_6]
    #print(p_i)
    T = 100  # K (Kelvin)
    valores_quadratura = [quad_gaussian(pi, xi, wi, T)
                                       for pi, xi, wi in zip(p_i, xi_6, wi_6)]

    U_para_B2 = [potencial_LJ(r, 140, 335)*(1/1000) for r in p_i]
    #print(valores_quadratura)
    Na = 6.022140e23
    B2 = [2*np.pi*Na*quadratura for quadratura in valores_quadratura]
    B2tot = sum(B2)
    #print(B2)
    #print(soma_B2)
    xii = 'xi [pm]'
    wii = 'wi [pm]'
    r = 'ri [pm]'
    U = 'U(ri) [kK]'
    B = 'B2(T) [cm³/mol]'
    print(f'{xii:^18}  {wii:^10} {r:^20} {U:^16} {B:^29}')
    for xi, wi, ri, Ui, Bi in zip(xi_6, wi_6, p_i, U_para_B2, B2):
        print(f'{xi:14.9f}  {wi:14.9f}  {ri:15.7f}   {Ui:15.7f}   {Bi:15.7f}')
    print(f'B2(T) = {B2tot} cm³/mol')

    # Cálculo de B3(T)
    # Região em que r2 <= r <= r3
    # ---------------------------
    #  Nesse caso nos vamos resolver a integral:
    # \int_{r1}^{r2}(1-e^(U(r)/kt))r^{2}dr, pelo método de aproximação conhecido
    # como quadratura gaussiana. Nesse caso vamos considerar n = 10, ou seja,
    # dez termos para a aproximação. Teremos então:
    # \int_{r1}^{r2}(1-e^(U(r)/kt))r^{2}dr ~ \sum_{i=1}^{10}w_i(1-e^(U(p_i)/kt))p_i^{2}.
    #  Em que p_i = 1/2*(r1+r2) + 1/2*(r2-r1)xi.
    #  Os valores de x_i e w_i são tabelados e foram retirados do sítio:
    # https://pomax.github.io/bezierinfo/legendre-gauss.html

    r3 = distancias1[-1]
    coef_wi_xi_10 = np.loadtxt('dados_coef10.dat', comments='#')
    wi_10 = coef_wi_xi_10[:, 1]
    xi_10 = coef_wi_xi_10[:, 2]

    #print(r2, r3)
    p_i = [calc_pontos_pi(r2, r3, xi) for xi in xi_10]
    #print(p_i)
    T = 100  # K (Kelvin)
    valores_quadratura = [quad_gaussian(pi, xi, wi, T)
                                       for pi, xi, wi in zip(p_i, xi_10, wi_10)]

    U_para_B3 = [potencial_LJ(r, 140, 335)*(1/1000) for r in p_i]
    #print(valores_quadratura)
    Na = 6.022140e23
    B3 = [2*np.pi*Na*quadratura for quadratura in valores_quadratura]
    B3tot = sum(B3)

    xii = 'xi [pm]'
    wii = 'wi [pm]'
    r = 'ri [pm]'
    U = 'U(ri) [kK]'
    B = 'B3(T) [cm³/mol]'
    print(f'{xii:^18}  {wii:^10} {r:^20} {U:^16} {B:^29}')
    for xi, wi, ri, Ui, Bi in zip(xi_10, wi_10, p_i, U_para_B3, B3):
        print(f'{xi:14.9f}  {wi:14.9f}  {ri:15.7f}   {Ui:15.7f}   {Bi:15.7f}')
    print(f'B3(T) = {B3tot} cm³/mol')


    # Cálculo de B4(T)
    # Região em que r3 <= r <= infinito
    # ---------------------------------
    # Nessa parte deveremos calcular a integral data pela relação
    # \int_{r3}^{\inf} U(r)/kT r^2 dr

    kb = 1         #j/K (joule por kelvin)
    kbT = kb * T
    r = symbols('r')
    # Aqui o potencial esta em kK.
    p = potencial_LJ(r, 140, 335)
    f = (1/kbT) * p * r**2

    infinito = 1325.0
    int_B4 = integrate(f, (r, r3, infinito))
    #print(int_B4)
    B4 = coef_virial1(int_B4)
    print(f'B4(T) = {B4} cm³/mol')

    Btot = B1 + B2tot + B3tot + B4
    print(f'Btot(T) = {Btot} cm³/mol')

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
