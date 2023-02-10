from sympy import *
import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.constants as const
import warnings
# suppress warnings
warnings.filterwarnings('ignore')

def pega_dados():
    with open('dados_rovib.txt', 'r') as f:
        dados = f.readlines()
    return dados


def pega_dado_como_float(dado):
    '''
    Essa função recebe um dado no formato: x = y,
    trata essa igualdade pegando o elemento y e
    converte esse valor de string para um dado
    float, retornando esse valor como float.
    '''
    return float(dado.split('=')[1].split('\n')[0].strip())


def massa_reduzida(massa1, massa2):
    '''
    A função recebe dois valores de massa em unidades
    atômicas, converte esses valores para kg (quilogramas)
    e calcula a massa reduzida para esse tipo de sistema.
    '''
    # Esse fator converte de 1 au = 1.660540199E-27
    fator_conversao = 1.660540199E-27
    soma = (massa1 + massa2)
    produto = (massa1 * massa2)

    return (produto / soma) * fator_conversao


def soma_de_massas(massa1, massa2):
    # Esse fator converte de 1 au = 1.660540199E-27
    fator_conversao = 1.660540199E-27
    soma = massa1 + massa2

    return soma * fator_conversao


def const_rot_eq(massa_reduzida, re):
    '''
    Essa função calcula a constante rotacional de equilíbrio,
    ou seja, calcula o valor de Be. Para calcular Be deve
    receber o valor da massa reduzida do sistema e a
    distância de equilíbrio do sistema.
    '''
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Massa reduzida em Kg
    mu = massa_reduzida
    # Momento de inércia
    I = mu * (re**2)
    Be = (h**2) / (8 * np.pi**2 * I)

    return Be


def convert_cm_to_joule(medida):
    fator_conversao = 1.98630e-23

    return medida * fator_conversao


def energia_rovib_tot(we, wexe, weye, Be, alfa_e, gama_e, n_quant_vib, j):
    '''
    n_quant_vib = número quântico vibracional.
    j = número quântico rotacional.
    '''
    vib = n_quant_vib + 0.5
    a = vib*we - (vib**2)*wexe + (vib**3)*weye
    b = (Be - vib*alfa_e + (vib**2)*gama_e) * (j*(j+1))

    return np.abs(a + b)


def temperatura_rotacional(massa_reduzida, re):
    '''
    Essa é a função de partição rotacional \theta_{rot}.

    Parâmetros
    ----------
    - massa reduzida (mu) em Kg.
    - distância de equilíbrio (re) em metros.
    '''
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    # Massa reduzida em Kg
    mu = massa_reduzida
    # Momento de inércia
    I = mu * (re**2)
    theta_rot = (h**2) / (8 * np.pi**2 * I * k)

    return theta_rot


def funcao_part_Mcquarie(massa_elementos, Temperatura, pressao, temperatura_rotacional,
                         we, gel, de):
    '''
    Essa é a função de partição de Mcquarie.

    Parâmetros
    ----------
    - massa_elementos: A soma das massas dos elementos em Kg.
    - Temperatura: A temperatura em K (kelvin).
    - pressao: A pressão em Pa (pascal).
    - temperatura_rotacional: Temperatura rotacional em K (kelvin).
    - we: Omega_e em Joule (j).
    - gel: degenerescência sem dimensão.
    - de: Joule (j).

    Retorno
    --------
    Retorna a função de partição adimencional.
    '''
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    T = Temperatura
    p = pressao
    theta_rot = temperatura_rotacional
    M = massa_elementos
    print(f'Massa reduzida = {mu}')
    M = soma_de_massas(massa1, massa2)
    print(f'M = {M}')

    a = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    b = (k * T**2) / (p * theta_rot)
    c = 1 / (1 - np.exp(-we/(k*T)))
    d = gel * np.exp(de/(k*T))
    Q_Mcquarie = a * b * c * d

    return Q_Mcquarie


def funcao_part_Mcquarie_sympy(massa_elementos, T, pressao, temperatura_rotacional,
                         we, gel, de):
    '''
    Essa é a função de partição de Mcquarie.

    Parâmetros
    ----------
    - massa_elementos: A soma das massas dos elementos em Kg.
    - Temperatura: A temperatura em K (kelvin).
    - pressao: A pressão em Pa (pascal).
    - temperatura_rotacional: Temperatura rotacional em K (kelvin).
    - we: Omega_e em Joule (j).
    - gel: degenerescência sem dimensão.
    - de: Joule (j).

    Retorno
    --------
    Retorna a função de partição adimencional.
    '''
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    p = pressao
    theta_rot = temperatura_rotacional
    M = massa_elementos

    a = ((2 * pi * M * k * T) / h**2) ** (3/2)
    b = (k * T**2) / (p * theta_rot)
    c = 1 / (1 - exp(-we/(k*T)))
    d = gel * exp(de/(k*T))
    Q_Mcquarie = a * b * c * d

    return log(Q_Mcquarie)

def energia_interna(derivada, Temperatura):
    R = const.gas_constant
    T = Temperatura

    return R * T**2 * derivada


def funcao_part_harmonica_Allison(massa_elementos, Temperatura, pressao, we, wexe, Be,
                                 alfa_e, gel, de):
    '''
    Essa é a função de partição Harmônica de Allison.

    Parâmetros mu = massa_reduzida(massa1, massa2)
    ----------
    - massa_elementos: A soma das massas dos elementos em Kg.
    - Temperatura: A temperatura em K (kelvin).
    - pressao: A pressão em Pa (pascal).
    - temperatura_rotacional: Temperatura rotacional em K (kelvin).
    - we: Omega_e em Joule (j).
    - wexe: Omega_e xe em Joule (j).
    - Be
    - gel: degenerescência sem dimensão.
    - de: Joule (j).

    Retorno
    --------
    Retorna a função de partição adimencional.
    '''
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    T = Temperatura
    p = pressao
    M = massa_elementos

    a = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = 1 / (1 - np.exp(-(we-wexe)/(k*T)))
    d = (k * T) / (Be - (alfa_e/2))
    e = gel * np.exp(de/(k*T))
    Q_Allison = a * b * c * d * e

    return Q_Allison


def func_particao_Allison(massa_elementos, Temperatura, pressao, we, wexe, Be,
                                 alfa_e, gel, de, nu):
    '''
    Essa é a função de pa M = massa_elementosrtição Harmônica de Allison.

    Parâmetros
    ----------
    - massa_elementos: A soma das massas dos elementos em Kg.
    - Temperatura: A temperatura em K (kelvin).
    - pressao: A pressão em Pa (pascal).
    - temperatura_rotacional: Temperatura rotacional em K (kelvin).
    - we: Omega_e em Joul M = massa_elementose (j).
    - wexe: Omega_e xe em Joule (j).
    - Be
    - gel: degenerescência sem dimensão.
    - de: Joule (j).
    - nu: número quantico vibracional.

    Retorno
    --------
    Retorna a função de partição adimencional.
    '''
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzm 1.380649x10-23 J K-1
    k = 1.380649e-23

    T = Temperatura
    p = pressao
    M = massa_elementos

    a = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = 0

    for i in range(0, nu):
        c += (np.exp(-((we - wexe)*i - wexe*i**2) / (k*T))) * (1/3 + (k*T)/(Be-(alfa_e/2)-alfa_e*i) + ( (8 * Be**3 * (k*T)**2) / (we**2 * (Be - (alfa_e/2) - alfa_e*i)**3) ) )

    d = gel * np.exp(de/(k*T))
    Q_Allison = a * b * c * d

    return Q_Allison

def funcao_particao_Foglia(massa_elementos, Temperatura, pressao, we, wexe, Be,
                                 alfa_e, gel, de, nu):
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    T = Temperatura
    p = pressao
    M = massa_elementos

    a = (k*T) / Be
    b = alfa_e / Be
    c = 0
    for i in range(0, nu):
        c += (np.exp((-(we*(i + 0.5) - wexe*(i + 0.5)**2))/(k*T))) * a * (1 + b*(i + 0.5)*(1 - np.exp(-(we*(i + 0.5) - wexe*(i + 0.5)**2)/(k*T))))
    d = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    e = (k * T) / p
    f = gel * np.exp(de/(k*T))
    print(f'Valor C = {c}')
    Q_Foglia = d * e * c * f

    return Q_Foglia


def funcao_particao_Heibbe_Scalabrini(massa_elementos, Temperatura, pressao, we,
                                     wexe, weye, Be, alfa_e, gama_e, gel, de, nu):
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    T = Temperatura
    p = pressao
    M = massa_elementos

    a = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = 0
    for i in range(0, nu):
        c += (np.exp(-((we - wexe + (3/4)*weye)*i + (-wexe + (2/3)*weye)*i**2 + (weye)*i**3) / (k*T))) * \
             (1/3 + ((k*T) / (Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2)) + \
             ((Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2) / (15*k*T)) + \
             (1/720)*((12*((Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2)**2)) / ((k*T)**2)) - \
             (1/720)*(((Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2)**3) / ((k*T)**3)))

    d = gel * np.exp(de/(k*T))
    Q_HS_tot = a * b * c * d

    return Q_HS_tot


def funcao_particao_Heibbe_Scalabrini_truncada(massa_elementos, Temperatura, pressao,
                                     we, wexe, weye, Be, alfa_e, gama_e, gel, de, nu):
    '''
    A distinção ente a função de partição Heibbe Scalabrini e Heibbe Scalabrini
    truncada, consiste no fato de que a última foi truncada na soma de
    Euler-Mclaurin na parte rotacional.
    '''
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    T = Temperatura
    p = pressao
    M = massa_elementos

    a = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = 0
    for i in range(0, nu):
        c += (np.exp( -((we - wexe + (3/4)*weye)*i + (-wexe + (2/3)*weye)*i**2 + (weye)*i**3) / (k*T))) * \
             (1/3 + ((k*T) / (Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2)))
    d = gel * np.exp(de/(k*T))
    Q_HS_truncada = a * b * c * d

    return Q_HS_truncada


def funcao_particao_Scalabrini_Rotor_Rigido(massa_elementos, Temperatura, pressao,
                                            we, wexe, weye, gel, de, nu, theta_rot):
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    T = Temperatura
    p = pressao
    M = massa_elementos

    a = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = T / theta_rot
    d = np.exp(-(we/2 - wexe/4 + weye/8) / (k*T))
    e = 0
    for i in range(0, nu):
        e += np.exp(-((we - wexe + 3/4*weye)*i + (-wexe + 3/2*weye)*i**2) / (k*T))
    f = gel * np.exp(de/(k*T))
    Q_S_rr = a * b * c * d * e * f

    return Q_S_rr


def funcao_particao_tietz(massa_reduzida, Temperatura, we, de, re, alfa_e):

    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    # Constante velocidade da luz
    c = const.speed_of_light
    print(c)

    mu = massa_reduzida
    T = Temperatura

    a1 = 2 * np.pi * c * we
    a2 = ((2*mu)/de)**0.5
    a3 = (32 * np.pi**4 * c**2 * mu**2 * re**3 * alfa_e * we) / (3 * h**2)
    a = (a1 * a2) - a3 + 1/re
    print(a)

    b = ((a/(np.pi*c*we)) * ((de / (2*np.pi))**0.5) - 1) * np.exp(a*re)
    print(b)

    c1 = ((mu*de) / ((h/2*np.pi)**2 * a**2 * b**2)) * (np.exp(2*a*re) - b**2)
    print(c1)

    d = 0.5 * (1 + (1 + ((8*np.pi*de*(np.exp(a*re)+b)**2) / (((h/2*np.pi)**2 * a**2 * b**2))))**0.5)
    print(d)

    e = c1/d - d/2
    print(e)

    f = (c1/(84 + 1 + d)) - ((84 + 1 + d)/2)
    print(f)

    g = ((h/(2*np.pi))**2 * a**2) / (2*mu)
    print(g)

    l = np.exp((g * e**2) / (k*T))
    print(l)
    m = np.exp((g*f**2)/(k*T))

    #   comp é uma variável para trabalhar somente com a parte real do número complexo
    # raiz(-1).
    comp = ((-1)**0.5).real

    n = (comp * ((g) / (k*T))**0.5)
    #print(n)

    Qtietz = 0.5 * np.exp((-de) / (k*T)) * (l - m)#* (l - m + (((np.pi*k*T)/g)**0.5) *  0.1)

    # (-comp * math.erf(n * e) + comp * math.erf(n * f))
    return Qtietz


if __name__ == '__main__':

    dados = pega_dados()

    for dado in dados:
        if 'm1' in dado:
            massa1 = pega_dado_como_float(dado)
        elif 'm2' in dado:
            massa2 = pega_dado_como_float(dado)
        elif 'p' in dado:
            pressao = pega_dado_como_float(dado)
        elif 'we_' in dado:
            we = pega_dado_como_float(dado)
        elif 'alfa_e' in dado:
            alfa_e = pega_dado_como_float(dado)
        elif 'wexe' in dado:
            wexe = pega_dado_como_float(dado)
        elif 'weye' in dado:
            weye = pega_dado_como_float(dado)
        elif 'gama_e' in dado:
            gama_e = pega_dado_como_float(dado)
        elif 'De' in dado:
            de = pega_dado_como_float(dado)
        elif 'Re' in dado:
            re = pega_dado_como_float(dado)
        elif 'gel' in dado:
            gel = pega_dado_como_float(dado)

    print(f'm1 = {massa1} ua')
    print(f'm2 = {massa2} ua')
    mu = massa_reduzida(massa1, massa2)
    print(f'Massa reduzida = {mu}')
    M = soma_de_massas(massa1, massa2)
    print(f'M = {M}')
    p = pressao
    print(f'pressão = {pressao} Pa')
    print()

    print(' --- Dados Espectroscópicos ---')
    print(f'we = {we} cm-1')
    print(f'alfa_e = {alfa_e} cm-1 ')
    print(f'wexe = {wexe} cm-1')
    print(f'weye = {weye} cm-1')
    print(f'gama_e = {gama_e} cm-1')
    print()
    print('-------------------------------')

    theta_rot = temperatura_rotacional(mu, re)
    print(f'Temperatura Rotacional\ntheta_rot = {theta_rot} K')
    Be = const_rot_eq(mu, re)
    print(f'Be = {Be}')
    print('-------------------------------')
    print()

    we = convert_cm_to_joule(we)
    alfa_e = convert_cm_to_joule(alfa_e)
    wexe = convert_cm_to_joule(wexe)
    weye = convert_cm_to_joule(gama_e)
    gama_e = convert_cm_to_joule(gama_e)

    # Número quântico rotacional
    j = 0
    # Números quânticos vibracionais para nu e nu + 1
    nu = 0
    nu1 = 1
    # Energias rovibracionais para nu e nu + 1~
    en_nu_j = 0
    en_nu1_j = 0
    # Lista de energias rovibracionais
    lista_en_nu_j = []
    # Lista de numeros quanticos vibracionais
    lista_nu = []
    while (np.abs((en_nu_j - de)) >= np.abs((en_nu1_j - de))):
        en_nu_j = energia_rovib_tot(we, wexe, weye, Be, alfa_e, gama_e, nu, j)
        en_nu1_j = energia_rovib_tot(we, wexe, weye, Be, alfa_e, gama_e, nu1, j)
        lista_en_nu_j.append(en_nu_j)
        lista_nu.append(nu)
        #print(en_nu_j)
        nu += 1
        nu1 += 1
        #print(nu)

    with open('Energia_rovib.txt', 'w') as f:
        for energia in lista_en_nu_j:
            print(energia, file=f)
    #print(en_nu_j)
    #print(nu)

    plt.plot(lista_nu, lista_en_nu_j)
    plt.plot(lista_nu, [de for i in range(len(lista_nu))])
    plt.xlabel(r"$\nu$")
    plt.ylabel(r"$E(\nu, j=0)$")
    plt.show()

    Temp = 298
    T = symbols('T')

    func_part_Macquarie = funcao_part_Mcquarie(M, Temp, p, theta_rot, we, gel, de)
    print(f'Função Partição Mcquarrie {func_part_Macquarie}')
    print('\n\n')

    func_part_Macquarie_sympy = funcao_part_Mcquarie_sympy(M, T, p, theta_rot, we, gel, de)
    pprint(func_part_Macquarie_sympy)
    print('\n\n')

    df_func_part_Macquarie_sympy = diff(func_part_Macquarie_sympy, T)
    pprint(df_func_part_Macquarie_sympy)

    #df_func_part_Macquarie_sympy = diff(func_part_Macquarie_sympy, T).evalf()
    #pprint(df_func_part_Macquarie_sympy)

    df_func_part_Macquarie_sympy = diff(func_part_Macquarie_sympy, T).evalf(subs={T: Temp})
    print(f'Derivada Mcquarie = {df_func_part_Macquarie_sympy}')

    U_Mcquarie = energia_interna(df_func_part_Macquarie_sympy, Temp)
    print(f'Energia interna Mcquarie = {U_Mcquarie}')

    '''
    func_part_harm_Allison = funcao_part_harmonica_Allison(M, Temp, p, we, wexe, Be, alfa_e, gel, de)
    print(f'Função Partição Hamônica de Allison {func_part_harm_Allison}')

    func_part_Allison = func_particao_Allison(M, Temp, p, we, wexe, Be, alfa_e, gel, de, nu)
    print(f'Função de Partição de Allison {func_part_Allison}')

    func_part_Foglia = funcao_particao_Foglia(M, Temp, p, we, wexe, Be, alfa_e, gel, de, nu)
    print(f'Função de Partição Vibracional de Foglia {func_part_Foglia}')

    func_part_HS_tot = funcao_particao_Heibbe_Scalabrini(M, Temp, p, we, wexe, weye, Be, alfa_e, gama_e, gel, de, nu)
    print(f'Função de Partição Heibbe Scalabrini Total {func_part_HS_tot}')

    func_part_HS_truc = funcao_particao_Heibbe_Scalabrini_truncada(M, Temp, p, we, wexe, weye, Be, alfa_e, gama_e, gel, de, nu)
    print(f'Função de Partição Heibbe Scalabrini Truncada {func_part_HS_truc}')

    func_part_scalabrini_rot_rig = funcao_particao_Scalabrini_Rotor_Rigido(M, Temp, p, we, wexe, weye, gel, de, nu, theta_rot)
    print(f'Função de Partição Scalabrini Rotor Rígido {func_part_scalabrini_rot_rig}')



    we = 2169.8129
    de = 1.801e-18
    alfa_e = 0.01750406
    wexe = 13.2883176
    weye = 0

    func_part_tietz = funcao_particao_tietz(mu, T, we, de, re, alfa_e)
    print(f'Função de Partição Tietz {func_part_tietz}')
    '''
