from sympy import *
from sympy.stats import *
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

#plt.style.use('fivethirtyeight')
sns.set(style="ticks")


def potencial_LJ(r, epsilon, sigma):
    '''
    Essa função recebe valores de distâncias entre dois átomos
    e retorna o valor da energia de interação.

    Entrada r/(angstron ou pm), epsilon/(eV, J, kK), sigma/(angstron ou pm)
    -------
    A unidade de medida da distância é em angstron ou pm, epsilon é eV, J ou KK
    e sigma é angstrom ou pm.

    Retorno
    -------
    A unidade de saida/retorno depende do valor de epsilon dado como entrada.
    '''
    e = epsilon
    s = sigma
    E = 4*e*(((s/r)**12) - ((s/r)**6))
    return E

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

def quad_gaussian(pi, xi, wi, T, r1, r2):
    '''
    A quadratura gaussiana retorna o valor valor em cm³.
    '''

    kb = 1         #j/K (joule por kelvin)
    kbT = kb * T
    ep = 140       #kK (lembre k = 1,38064x10^-23 j/K -> kK = 1,38064x10^-23 j)
    si = 335       #pm (picometros)
    potencial = potencial_LJ(pi, ep, si)
    # Aqui o potencial esta em kK.
    #print(potencial)
    a = -potencial/kbT
    #print(np.exp(a))
    c = wi * (1 - np.exp(a)) * (pi**2) * (r2 - r1)/2

    return c * 1e-30


if __name__ == "__main__":
    distancias = np.arange(3.10, 8.20, 0.01)

    dados_experimentais = np.loadtxt('dados_experimentais_Ar.dat',
                                      comments='#')

    T_exp = dados_experimentais[:, 0]
    B_exp_Ar = dados_experimentais[:, 1]

    #Parâmetros para o potencial de Lennard-Jones Ar2
    #ref:
    #epsilon = 0.0103332 #eV
    #sigma = 3.40        #angstrom
    #energias = [potencial_LJ(r, epsilon, sigma) for r in distancias]

    #Parâmetros para o potencial de Lennard-Jones Improve Ar2
    #ref: A Spectroscopic Validation of the Improved Lennard–Jones Model
    #de = 12.343  #meV
    #req = 3.76   #angstrom
    #beta = 9.74
    #energias1 = [improve_LJ(r, de, req, beta) for r in distancias]

    # Parâmetros para o potencial de Lennard-Jones Ar2.
    # -------------------------- referência -------------------------------------
    # nome: Determining Intermolecular Potentials from second virial coefficients
    # autor: Brian P. Reid
    # ---------------------------------------------------------------------------

    epsilon1 = 140 #kK (lembre k = 1,38064x10^-23 j/K -> kK = 1,38064x10^-23 j)
    sigma1 = 335   #pm (picometros)
    distancias1 = np.arange(305.21, 1303.91)
    energias2 = [potencial_LJ(r, epsilon1, sigma1)*(1/(epsilon1))
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

    new_distancias1 = [(r/298) for r in distancias1]

    plt.plot(new_distancias1, energias2, color='b', linewidth=2.5, label='LJ')
    plt.legend(loc='upper right', shadow=False, fontsize='large',
               bbox_to_anchor=(0.98, 0.98), frameon=False)
    plt.ylabel(r'Energia ($U=U/\epsilon$)')
    plt.xlabel(r'R ($r/r_{min}$)')
    plt.title(r'Interação $Ar_{2}$')
    plt.show()

    # Para calcular o coeficiente do virial total vamos repartir o coeficiente
    #em quatro partes de acordo com o a distância de interação.
    # O valor do coeficiente do virial total é dado pela soma abaixo:
    # B(T) = B1(T) + B2(T) + B3(T) + B4(T)


    # Como buscamos calcular o valor de energia para vária temperaturas, então:
    B_tot = []
    for T in T_exp:
        # Cálculo de B1(T)
        # Região em que r é pequeno
        # -------------------------
        # Se r varia de 0 até r1, então podemos calcular B1(T) assumindo que a
        #integral que devemos resolver é \int_{0}^{r1}r^{2}dr.

        r = symbols('r')
        #r1 = distancias1[0]
        r1 = 305.2
        #print(r1)
        #print(f'r inicial = {r1} | rmin = {r_min1}')
        int1 = integrate(r**2, (r, 0, r1))
        #print(f'Integral = {int1}')
        B1 = coef_virial1(int1)
        print(f'B1(T) = {B1} cm³/mol')
        print()
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
        #r2 = r_min1
        r2 = 375.78
        #print()
        #print(r1, r2)
        #print()
        p_i = [calc_pontos_pi(r1, r2, xi) for xi in xi_6]
        #print(p_i)
        valores_quadratura = [quad_gaussian(pi, xi, wi, T, r1, r2)
                                       for pi, xi, wi in zip(p_i, xi_6, wi_6)]

        U_para_B2 = [potencial_LJ(r, 119.92, 298) for r in p_i]
        #print(valores_quadratura)
        Na = 6.022140e23
        B2 = [2*np.pi*Na*quadratura for quadratura in valores_quadratura]
        B2tot = sum(B2)
        #print(B2)
        #print(soma_B2)
        print()
        xii = 'xi'
        wii = 'ci'
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

        #r3 = distancias1[-1]
        r3 = 1303.91
        coef_wi_xi_10 = np.loadtxt('dados_coef10.dat', comments='#')
        wi_10 = coef_wi_xi_10[:, 1]
        xi_10 = coef_wi_xi_10[:, 2]

        #print(r2, r3)

        p_i = [calc_pontos_pi(r2, r3, xi) for xi in xi_10]
        #print(p_i)
        valores_quadratura = [quad_gaussian(pi, xi, wi, T, r2, r3)
                                       for pi, xi, wi in zip(p_i, xi_10, wi_10)]

        U_para_B3 = [potencial_LJ(r, 119.92, 341.1) for r in p_i]
        #print(valores_quadratura)
        Na = 6.022140e23
        B3 = [2*np.pi*Na*quadratura for quadratura in valores_quadratura]
        B3tot = sum(B3)

        xii = 'xi'
        wii = 'ci'
        r = 'ri [pm]'
        U = 'U(ri) [kK]'
        B = 'B3(T) [cm³/mol]'
        print(f'{xii:^18}  {wii:^10} {r:^20} {U:^16} {B:^29}')
        for xi, wi, ri, Ui, Bi in zip(xi_10, wi_10, p_i, U_para_B3, B3):
            print(f'{xi:14.9f}  {wi:14.9f}  {ri:15.7f}   {Ui:15.7f}   {Bi:15.7f}')
        print(f'B3(T) = {B3tot} cm³/mol')
        print()

        # Cálculo de B4(T)
        # Região em que r3 <= r <= infinito
        # ---------------------------------
        # Nessa parte deveremos calcular a integral data pela relação
        # \int_{r3}^{\inf} U(r)/kT r^2 dr

        kb = 1         #j/K (joule por kelvin)
        kbT = kb * T
        r = symbols('r')
        # Aqui o potencial esta em kK.
        p = potencial_LJ(r, 119.92, 298)
        f = (1/kbT) * p * r**2

        infinito = 1625.0
        int_B4 = integrate(f, (r, r3, infinito))
        #print(int_B4)
        B4 = coef_virial1(int_B4)
        print(f'B4(T) = {B4} cm³/mol')
        print()

        Btot = B1 + B2tot + B3tot + B4
        print(f'Btot(T) = {Btot} cm³/mol')

        B_tot.append(Btot)

    # Plot do coeficiente do virial experimental com calculado
    plt.scatter(T_exp, B_exp_Ar, marker='o', color='b', linewidth=2.5,
             label='B(T) Experimental')
    plt.plot(T_exp, B_tot, color='red', linewidth=2.5,
             label='B(T) Calculado')
    plt.legend(loc='upper right', shadow=False, fontsize='large',
               bbox_to_anchor=(0.98, 0.28), frameon=False)
    plt.ylabel(r'B(T)')
    plt.xlabel(r'Temperatura ($K$)')
    plt.title(r'Coeficiente do Virial $Ar_{2}$')
    plt.show()
