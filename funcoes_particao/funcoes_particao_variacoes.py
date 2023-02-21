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


def energia_interna(derivada, Temperatura):
    R = const.gas_constant
    T = Temperatura

    return R * T**2 * derivada


def entalpia(energia_interna, Temperatura):
    '''
    A função calcula a entalpia do sistema a partir
    da energia interna e temperatura retornando
    o valor da entalpia em joule / mol.
    H = U + RT + 1084.65x10^3 J/mol
    '''
    U = energia_interna
    T = Temperatura
    R = const.gas_constant
    H = U + R * T + 1084.65e3

    return H


def energia_Gibbs(entalpia, entropia, Temperatura):
    T = Temperatura
    G = entalpia - T * entropia

    return G


def entropia(derivada, Temperatura, funcao_particao):
    R = const.gas_constant
    T = Temperatura
    S = R * T * derivada + R * funcao_particao

    return S


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


def funcao_part_monomero_sympy(massa_elemento, T, pressao, gel, de=0):
    '''
    Essa é a função de partição de Mcquarie.

    Parâmetros
    ----------
    - massa_elementos: A massas dos elementos em Kg.
    - Temperatura: A temperatura em K (kelvin).
    - pressao: A pressão em Pa (pascal).
    - gel: degenerescência eletrônica sem dimensão.
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
    M = massa_elemento

    a = ((2 * pi * M * k * T) / h**2) ** (3/2)
    b = gel * exp(de/(k*T))
    Q_Mcquarie = a * b

    return log(Q_Mcquarie)


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


def funcao_part_harmonica_Allison_sympy(massa_elementos, T, pressao, we, wexe,
                                        Be, alfa_e, gel, de):
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
    p = pressao
    M = massa_elementos

    a = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = 1 / (1 - exp(-(we-wexe)/(k*T)))
    d = (k * T) / (Be - (alfa_e/2))
    e = gel * exp(de/(k*T))
    Q_Allison = a * b * c * d * e

    return log(Q_Allison)


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
        c += (np.exp(-((we - wexe)*i - wexe*i**2) / (k*T))) * \
             (1/3 + (k*T)/(Be-(alfa_e/2)-alfa_e*i) + \
             ((8 * Be**3 * (k*T)**2) / (we**2 * (Be - (alfa_e/2) - alfa_e*i)**3)))

    d = gel * np.exp(de/(k*T))
    Q_Allison = a * b * c * d

    return Q_Allison

def func_particao_Allison_sympy(massa_elementos, T, pressao, we, wexe, Be,
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

    p = pressao
    M = massa_elementos

    a = ((2 * pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = 0

    for i in range(0, nu):
        c += (exp(-((we - wexe)*i - wexe*i**2) / (k*T))) * \
             (1/3 + (k*T)/(Be-(alfa_e/2)-alfa_e*i) + \
             ((8 * Be**3 * (k*T)**2) / (we**2 * (Be - (alfa_e/2) - alfa_e*i)**3)))

    d = gel * exp(de/(k*T))
    Q_Allison = a * b * c * d

    return log(Q_Allison)


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
        c += (np.exp((-(we*(i + 0.5) - wexe*(i + 0.5)**2))/(k*T))) * a * \
             (1 + b*(i + 0.5)*(1 - np.exp(-(we*(i + 0.5) - wexe*(i + 0.5)**2)/(k*T))))
    d = ((2 * np.pi * M * k * T) / h**2) ** (3/2)
    e = (k * T) / p
    f = gel * np.exp(de/(k*T))
    print(f'Valor C = {c}')
    Q_Foglia = d * e * c * f

    return Q_Foglia


def funcao_particao_Foglia_sympy(massa_elementos, T, pressao, we, wexe, Be,
                                 alfa_e, gel, de, nu):
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    p = pressao
    M = massa_elementos

    a = (k*T) / Be
    b = alfa_e / Be
    c = 0
    for i in range(0, nu):
        c += (exp((-(we*(i + 0.5) - wexe*(i + 0.5)**2))/(k*T))) * a * \
             (1 + b*(i + 0.5)*(1 - exp(-(we*(i + 0.5) - wexe*(i + 0.5)**2)/(k*T))))
    d = ((2 * pi * M * k * T) / h**2) ** (3/2)
    e = (k * T) / p
    f = gel * exp(de/(k*T))
    #print(f'Valor C = {c}')
    Q_Foglia = d * e * c * f

    return log(Q_Foglia)


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


def funcao_particao_Heibbe_Scalabrini_sympy(massa_elementos, T, pressao, we,
                                     wexe, weye, Be, alfa_e, gama_e, gel, de, nu):
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    p = pressao
    M = massa_elementos

    a = ((2 * pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = 0
    for i in range(0, nu):
        c += (exp(-((we - wexe + (3/4)*weye)*i + (-wexe + (2/3)*weye)*i**2 + (weye)*i**3) / (k*T))) * \
             (1/3 + ((k*T) / (Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2)) + \
             ((Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2) / (15*k*T)) + \
             (1/720)*((12*((Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2)**2)) / ((k*T)**2)) - \
             (1/720)*(((Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2)**3) / ((k*T)**3)))

    d = gel * exp(de/(k*T))
    Q_HS_tot = a * b * c * d

    return log(Q_HS_tot)


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


def funcao_particao_Heibbe_Scalabrini_truncada_sympy(massa_elementos, T, pressao,
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
    p = pressao
    M = massa_elementos

    a = ((2 * pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = 0
    for i in range(0, nu):
        c += (exp( -((we - wexe + (3/4)*weye)*i + (-wexe + (2/3)*weye)*i**2 + (weye)*i**3) / (k*T))) * \
             (1/3 + ((k*T) / (Be - alfa_e/2 + gama_e/4 - alfa_e*i + gama_e*i + gama_e*i**2)))
    d = gel * exp(de/(k*T))
    Q_HS_truncada = a * b * c * d

    return log(Q_HS_truncada)


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


def funcao_particao_Scalabrini_Rotor_Rigido_sympy(massa_elementos, T, pressao,
                                            we, wexe, weye, gel, de, nu, theta_rot):
    # Constante de Planck 6.62607015e-34 J Hz^-1
    h = 6.62607015e-34
    # Constante de Boltzmann 1.380649x10-23 J K-1
    k = 1.380649e-23
    p = pressao
    M = massa_elementos

    a = ((2 * pi * M * k * T) / h**2) ** (3/2)
    b = (k * T) / p
    c = T / theta_rot
    d = exp(-(we/2 - wexe/4 + weye/8) / (k*T))
    e = 0
    for i in range(0, nu):
        e += exp(-((we - wexe + 3/4*weye)*i + (-wexe + 3/2*weye)*i**2) / (k*T))
    f = gel * exp(de/(k*T))
    Q_S_rr = a * b * c * d * e * f

    return log(Q_S_rr)


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
            # O valor de entrada da energia de deve ser em joule.
            de = pega_dado_como_float(dado)
        elif 'Re' in dado:
            # O valor de entrada de re dever ser angstrom.
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
    print(f'Be = {Be} J')
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
    # Energias rovibracionais para nu e nu + 1
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

    with open('Energia_rovib.txt', 'w') as f1, \
         open('niveis_vibracionais.txt', 'w') as f2:
        for energia in lista_en_nu_j:
            print(energia, file=f1)
        for niveis_vib in lista_nu:
            print(niveis_vib, file=f2)
    #print(en_nu_j)
    #print(nu)

    plt.plot(lista_nu, lista_en_nu_j)
    plt.plot(lista_nu, [de for i in range(len(lista_nu))])
    plt.xlabel(r"$\nu$")
    plt.ylabel(r"$E(\nu, j=0)$")
    plt.show()

    Temp_inicial = 298
    Temp_final = 6000

    faixa_Temp = []

    dados_entalpia_Macquarie = []
    dados_entalpia_Allison_harm = []
    dados_entalpia_Allison = []
    dados_entalpia_Foglia = []
    dados_entalpia_H_S = []
    dados_entalpia_H_S_trunc = []
    dados_entalpia_H_S_rot_rig = []

    dados_entropia_Macquarie = []
    dados_entropia_Allison_harm = []
    dados_entropia_Allison = []
    dados_entropia_Foglia = []
    dados_entropia_H_S = []
    dados_entropia_H_S_trunc = []
    dados_entropia_H_S_rot_rig = []

    dados_cp_Mcquarie = []
    dados_cp_Allison_harm = []
    dados_cp_Allison = []
    dados_cp_Foglia = []
    dados_cp_H_S = []
    dados_cp_H_S_trunc = []
    dados_cp_H_S_rot_rig = []

    dados_G_Mcquarie = []
    dados_G_Allison_harm = []
    dados_G_Allison = []
    dados_G_Foglia = []
    dados_G_H_S = []
    dados_G_H_S_trunc = []
    dados_G_H_S_rot_rig = []

    with open('energia_interna.csv', 'w') as f1, \
         open('entalpia.csv', 'w') as f2, \
         open('entropia.csv', 'w') as f3, \
         open('capacidade_calorifica.csv', 'w') as f4, \
         open('energia_de_gibs.csv', 'w') as f5:

         cabecalho_U = 'Temperatura,U_Mcquarie,U_Allison_Harmonica,U_Allison,\
         U_Foglia,U_Heibbe_Scalabrini,U_Heibbe_Scalabrini_Trunc,\
         U_Heibbe_Scalabrini_ROTOR-RIG'

         cabecalho_H = 'Temperatura,H_Mcquarie,H_Allison_Harmonica,H_Allison,\
         H_Foglia,H_Heibbe_Scalabrini,H_Heibbe_Scalabrini_Trunc,\
         H_Heibbe_Scalabrini_ROTOR-RIG'

         cabecalho_S = 'Temperatura,S_Mcquarie,S_Allison_Harmonica,S_Allison,\
         S_Foglia,S_Heibbe_Scalabrini,S_Heibbe_Scalabrini_Trunc,\
         S_Heibbe_Scalabrini_ROTOR-RIG'

         cabecalho_Cp = 'Temperatura,Cp_Mcquarie,Cp_Allison_Harmonica,Cp_Allison,\
         Cp_Foglia,S_Heibbe_Scalabrini,Cp_Heibbe_Scalabrini_Trunc,\
         Cp_Heibbe_Scalabrini_ROTOR-RIG'

         cabecalho_G = 'Temperatura,G_Mcquarie,G_Allison_Harmonica,G_Allison,\
         G_Foglia,G_Heibbe_Scalabrini,G_Heibbe_Scalabrini_Trunc,\
         G_Heibbe_Scalabrini_ROTOR-RIG'


         print(cabecalho_U, end='\n', file=f1)
         print(cabecalho_H, end='\n', file=f2)
         print(cabecalho_S, end='\n', file=f3)
         print(cabecalho_Cp, end='\n', file=f4)
         print(cabecalho_G, end='\n', file=f5)

         for Temp in range(Temp_inicial, Temp_final, 600):
            faixa_Temp.append(Temp)
            '''
            func_part_Macquarie = funcao_part_Mcquarie(M, Temp, p, theta_rot, we, gel, de)
            print(f'Função Partição Mcquarrie {func_part_Macquarie}')
            print('\n\n')

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

            T = symbols('T')

            ### Mcquarie
            # Nas linhas abaixo vamos definir as funções de partição desiginadas por Mcquarie
            func_part_Mcquarie_sympy = funcao_part_Mcquarie_sympy(M, T, p, theta_rot, we, gel, de)

            # Funcoes de particao para os monomeros
            func_part_monomero1_sympy = funcao_part_monomero_sympy(massa1, T, p, gel)
            func_part_monomero2_sympy = funcao_part_monomero_sympy(massa2, T, p, gel)

            # derivadas para as funcoes de particao de Mcquarie não avaliadas
            df_func_part_Mcquarie_sympy_nao_aval = diff(func_part_Mcquarie_sympy, T)

            # Derivadas das funcoes de particao dos monomeros não avaliadas
            df_func_part_sympy_nao_aval_mono1 = diff(func_part_monomero1_sympy, T)
            df_func_part_sympy_nao_aval_mono2 = diff(func_part_monomero2_sympy, T)

            #pprint(df_func_part_Macquarie_sympy_nao_aval)

            #df_func_part_Macquarie_sympy = diff(func_part_Macquarie_sympy, T).evalf()
            #pprint(df_func_part_Macquarie_sympy)

            # Derivadas avaliadas
            df_func_part_Mcquarie_sympy = diff(func_part_Mcquarie_sympy, T).evalf(subs={T: Temp})
            df_func_part_sympy_mono1 = diff(func_part_monomero1_sympy, T).evalf(subs={T: Temp})
            df_func_part_sympy_mono2 = diff(func_part_monomero2_sympy, T).evalf(subs={T: Temp})

            print(f'Derivada Mcquarie = {df_func_part_Mcquarie_sympy}')
            print(f'Derivada Mcquarie Mono 1 = {df_func_part_sympy_mono1}')
            print(f'Derivada Mcquarie Mono 2 = {df_func_part_sympy_mono2}')
            print()

            U_Mcquarie = energia_interna(df_func_part_Mcquarie_sympy, Temp)
            U_mono1 = energia_interna(df_func_part_sympy_mono1, Temp)
            U_mono2 = energia_interna(df_func_part_sympy_mono2, Temp)

            U_Mcquarie_nao_aval =  energia_interna(df_func_part_Mcquarie_sympy_nao_aval, T)
            U_nao_aval_mono1 =  energia_interna(df_func_part_sympy_nao_aval_mono1, T)
            U_nao_aval_mono2 =  energia_interna(df_func_part_sympy_nao_aval_mono2, T)

            delta_U_Mcquarie = U_Mcquarie - (U_mono1 + U_mono2)

            delta_U_Mcquarie_nao_aval = U_Mcquarie_nao_aval - (U_nao_aval_mono1 + U_nao_aval_mono2)

            #print(f'Energia interna Mcquarie = {U_Mcquarie}')
            #print(f'Energia interna Mcquarie Mono1 = {U_Mcquarie_mono1}')
            #print(f'Energia interna Mcquarie Mono2 = {U_Mcquarie_mono2}')
            print(f'Variação Energia interna Mcquarie = {delta_U_Mcquarie}')

            H_Mcquarie = entalpia(U_Mcquarie, Temp)
            H_mono1 = entalpia(U_mono1, Temp)
            H_mono2 = entalpia(U_mono2, Temp)
            delta_H_Mcquarie = H_Mcquarie - (H_mono1 + H_mono2)

            dados_entalpia_Macquarie.append(delta_H_Mcquarie)

            print(f'Variação Entalpia Mcquarie = {delta_H_Mcquarie}')

            Cp_Mcquarie = diff(entalpia(U_Mcquarie_nao_aval, T)).evalf(subs={T: Temp})
            Cp_mono1 = diff(entalpia(U_mono1, T)).evalf(subs={T: Temp})
            Cp_mono2 = diff(entalpia(U_mono2, T)).evalf(subs={T: Temp})
            delta_Cp_Mcquarie = Cp_Mcquarie - (Cp_mono1 + Cp_mono2)
            dados_cp_Mcquarie.append(delta_Cp_Mcquarie)
            print(f'Variação Cp Mcquarie = {delta_Cp_Mcquarie}')

            S_Mcquarie = entropia(df_func_part_Mcquarie_sympy, Temp,
                                  funcao_part_Mcquarie_sympy(M, T, p, theta_rot,
                                  we, gel, de).evalf(subs={T: Temp}))

            S_mono1 = entropia(df_func_part_sympy_mono1, Temp,
                              funcao_part_monomero_sympy(massa1, T, p, gel).evalf(subs={T: Temp}))

            S_mono2 = entropia(df_func_part_sympy_mono2, Temp,
                              funcao_part_monomero_sympy(massa2, T, p, gel).evalf(subs={T: Temp}))

            delta_S_Mcquarie = S_Mcquarie - (S_mono1 + S_mono2)
            dados_entropia_Macquarie.append(delta_S_Mcquarie)
            print(f'Variação da Entropia Mcquarie = {delta_S_Mcquarie}')

            delta_G_Mcquarie = energia_Gibbs(delta_H_Mcquarie, delta_S_Mcquarie, Temp)
            dados_G_Mcquarie.append(delta_G_Mcquarie)
            print(f'Energia de Gibbs = {delta_G_Mcquarie}')

            print('-'*60)
            print('\n')


            ### Allison Harmonica

            func_part_harm_Allison_sympy = funcao_part_harmonica_Allison_sympy(M, T, p, we, wexe,
                                                                          Be, alfa_e, gel, de)
            #pprint(func_part_harm_Allison_sympy)
            #print('\n')

            df_func_part_Allison_nao_aval = diff(func_part_harm_Allison_sympy, T)
            df_func_part_Allison_harm = diff(func_part_harm_Allison_sympy, T).evalf(subs={T: Temp})
            print(f'Derivada Allison Harm = {df_func_part_Allison_harm}')

            U_Allison_harm_nao_aval = energia_interna(df_func_part_Allison_nao_aval, T)
            U_Allison_harm = energia_interna(df_func_part_Allison_harm, Temp)
            print(f'Energia interna Allison Harm = {U_Allison_harm}')

            delta_U_Allison_harm = U_Allison_harm - (U_mono1 + U_mono2)
            delta_U_Allison_harm_nao_aval = U_Allison_harm_nao_aval - (U_nao_aval_mono1 + U_nao_aval_mono2)

            H_Allison_harm = entalpia(U_Allison_harm, Temp)
            delta_H_Allison_harm = H_Allison_harm - (H_mono1 + H_mono2)
            dados_entalpia_Allison_harm.append(delta_H_Allison_harm)
            print(f'Entalpia Allison Harm = {delta_H_Allison_harm}')


            Cp_Allison_harm = diff(entalpia(U_Allison_harm_nao_aval, T)).evalf(subs={T: Temp})
            delta_Cp_Allison_harm = Cp_Allison_harm - (Cp_mono1 + Cp_mono2)
            dados_cp_Allison_harm.append(delta_Cp_Allison_harm)
            print(f'Cp Allison Harm. = {delta_Cp_Allison_harm}')

            S_Allison_harm = entropia(df_func_part_Allison_harm, Temp,
                                     funcao_part_harmonica_Allison_sympy(M, T, p, we,
                                     wexe, Be, alfa_e, gel, de).evalf(subs={T: Temp}))
            delta_S_Allison_harm = S_Allison_harm - (S_mono1 + S_mono2)
            dados_entropia_Allison_harm.append(delta_S_Allison_harm)
            print(f'Entropia Allison Harm = {delta_S_Allison_harm}')

            delta_G_Allison_harm = energia_Gibbs(delta_H_Allison_harm, delta_S_Allison_harm, Temp)
            dados_G_Allison_harm.append(delta_G_Allison_harm)
            print(f'Energia de Gibbs = {delta_G_Allison_harm}')

            print('-'*60)
            print('\n')


            '''
            ### Allison

            func_part_Allison_sympy = func_particao_Allison_sympy(M, T, pressao, we, wexe,
                                                                  Be, alfa_e, gel, de, nu)
            #pprint(func_part_Allison_sympy)
            #print('\n')

            df_func_part_Allison_nao_aval = diff(func_part_Allison_sympy, T)
            df_func_part_Allison = diff(func_part_Allison_sympy, T).evalf(subs={T: Temp})
            print(f'Derivada Allison = {df_func_part_Allison}')

            U_Allison_nao_aval = energia_interna(df_func_part_Allison_nao_aval, T)
            U_Allison = energia_interna(df_func_part_Allison, Temp)
            print(f'Energia interna Allison = {U_Allison}')

            if Temp == Temp_inicial:
                H_Allison_298 = entalpia(U_Allison, Temp_inicial)
                H_Allison = entalpia(U_Allison, Temp)
                dados_entalpia_Allison.append(H_Allison - H_Allison_298)
                print(f'Entalpia Allison Harm = {H_Allison}')
            else:
                H_Allison = entalpia(U_Allison, Temp)
                dados_entalpia_Allison.append(H_Allison - H_Allison_298)
                print(f'Entalpia Allison Harm = {H_Allison}')

            Cp_Allison = diff(entalpia(U_Allison_nao_aval, T)).evalf(subs={T: Temp})
            dados_cp_Allison.append(Cp_Allison)
            print(f'Cp Allison = {Cp_Allison}')

            S_Allison = entropia(df_func_part_Allison, Temp,
                                 func_particao_Allison_sympy(M, T, pressao, we, wexe,
                                 Be, alfa_e, gel, de, nu).evalf(subs={T: Temp}))
            dados_entropia_Allison.append(S_Allison)
            print(f'Entropia Allison  = {S_Allison}')

            G_Allison = energia_Gibbs(H_Allison, S_Allison, Temp)
            dados_G_Allison.append(G_Allison)
            print(f'Energia de Gibbs = {G_Allison}')

            print('-'*60)
            print('\n')


            ### Foglia

            func_part_Foglia_sympy = funcao_particao_Foglia_sympy(M, T, pressao, we, wexe, Be,
                                                                  alfa_e, gel, de, nu)
            #pprint(func_part_Foglia_sympy)
            #print('\n')

            df_func_part_Foglia_nao_aval = diff(func_part_Foglia_sympy, T)
            df_func_part_Foglia = diff(func_part_Foglia_sympy, T).evalf(subs={T: Temp})
            print(f'Derivada Foglia = {df_func_part_Foglia}')

            U_Foglia_nao_aval = energia_interna(df_func_part_Foglia_nao_aval, T)
            U_Foglia = energia_interna(df_func_part_Foglia, Temp)
            print(f'Energia interna Foglia = {U_Foglia}')

            if Temp == Temp_inicial:
                H_Foglia_298 = entalpia(U_Foglia, Temp_inicial)
                H_Foglia = entalpia(U_Foglia, Temp)
                dados_entalpia_Foglia.append(H_Foglia - H_Foglia_298)
                print(f'Entalpia Foglia = {H_Foglia}')
            else:
                H_Foglia = entalpia(U_Foglia, Temp)
                dados_entalpia_Foglia.append(H_Foglia - H_Foglia_298)
                print(f'Entalpia Foglia = {H_Foglia}')

            Cp_Foglia = diff(entalpia(U_Foglia_nao_aval, T)).evalf(subs={T: Temp})
            dados_cp_Foglia.append(Cp_Foglia)
            print(f'Cp Foglia = {Cp_Foglia}')

            S_Foglia = entropia(df_func_part_Foglia, Temp,
                                funcao_particao_Foglia_sympy(M, T, pressao, we, wexe, Be,
                                alfa_e, gel, de, nu).evalf(subs={T: Temp}))
            dados_entropia_Foglia.append(S_Foglia)
            print(f'Entropia Foglia  = {S_Foglia}')

            G_Foglia = energia_Gibbs(H_Foglia, S_Foglia, Temp)
            dados_G_Foglia.append(G_Foglia)
            print(f'Energia de Gibbs = {G_Foglia}')

            print('-'*60)
            print('\n')


            ### Heibbe Scalabrini

            func_part_H_S = funcao_particao_Heibbe_Scalabrini_sympy(M, T, pressao, we,
                                                 wexe, weye, Be, alfa_e, gama_e, gel,
                                                 de, nu)

            #pprint(func_part_H_S)
            #print('\n')

            df_func_part_H_S_nao_aval = diff(func_part_H_S, T)
            df_func_part_H_S = diff(func_part_H_S, T).evalf(subs={T: Temp})
            print(f'Derivada Heibbe-Scalabrini = {df_func_part_Foglia}')

            U_H_S_nao_aval = energia_interna(df_func_part_H_S_nao_aval, T)
            U_H_S = energia_interna(df_func_part_H_S, Temp)
            print(f'Energia interna Heibbe-Scalabrini = {U_H_S}')

            if Temp == Temp_inicial:
                H_H_S_298 = entalpia(U_H_S, Temp_inicial)
                H_H_S = entalpia(U_H_S, Temp)
                dados_entalpia_H_S.append(H_H_S - H_H_S_298)
                print(f'Entalpia Heibbe-Scalabrini = {H_H_S}')
            else:
                H_H_S = entalpia(U_H_S, Temp)
                dados_entalpia_H_S.append(H_H_S - H_H_S_298)
                print(f'Entalpia Heibbe-Scalabrini = {H_H_S}')

            Cp_H_S = diff(entalpia(U_H_S_nao_aval, T)).evalf(subs={T: Temp})
            dados_cp_H_S.append(Cp_H_S)
            print(f'Cp Heibbe-Scalabrini = {Cp_H_S}')

            S_H_S = entropia(df_func_part_H_S, Temp,
                             funcao_particao_Heibbe_Scalabrini_sympy(M, T, pressao, we,
                             wexe, weye, Be, alfa_e, gama_e, gel, de, nu).evalf(subs={T: Temp}))
            dados_entropia_H_S.append(S_H_S)
            print(f'Entropia Heibbe-Scalabrini  = {S_H_S}')

            G_H_S = energia_Gibbs(H_H_S, S_H_S, Temp)
            dados_G_H_S.append(G_H_S)
            print(f'Energia de Gibbs = {G_H_S}')

            print('-'*60)
            print('\n')


            ### Heibbe Scalabrini Truncada

            func_part_H_S_trunc = funcao_particao_Heibbe_Scalabrini_truncada_sympy(M, T,
                                                 pressao, we, wexe, weye, Be, alfa_e,
                                                 gama_e, gel, de, nu)

            #pprint(func_part_H_S_trunc)
            #print('\n')
            df_func_part_H_S_trunc_nao_aval = diff(func_part_H_S_trunc, T)
            df_func_part_H_S_trunc = diff(func_part_H_S_trunc, T).evalf(subs={T: Temp})
            print(f'Derivada Heibbe-Scalabrini Truncada = {df_func_part_H_S_trunc}')

            U_H_S_trunc_nao_aval = energia_interna(df_func_part_H_S_trunc_nao_aval, T)
            U_H_S_trunc = energia_interna(df_func_part_H_S_trunc, Temp)
            print(f'Energia interna Heibbe-Scalabrini Truncada = {U_H_S_trunc}')

            if Temp == Temp_inicial:
                H_H_S_trunc_298 = entalpia(U_H_S_trunc, Temp_inicial)
                H_H_S_trunc = entalpia(U_H_S_trunc, Temp)
                dados_entalpia_H_S_trunc.append(H_H_S_trunc - H_H_S_trunc_298)
                print(f'Entalpia Heibbe-Scalabrini Truncada = {H_H_S_trunc}')
            else:
                H_H_S_trunc = entalpia(U_H_S_trunc, Temp)
                dados_entalpia_H_S_trunc.append(H_H_S_trunc - H_H_S_trunc_298)
                print(f'Entalpia Heibbe-Scalabrini Truncada = {H_H_S_trunc}')

            Cp_H_S_trunc = diff(entalpia(U_H_S_trunc_nao_aval, T)).evalf(subs={T: Temp})
            dados_cp_H_S_trunc.append(Cp_H_S_trunc)
            print(f'Cp Heibbe-Scalabrini Truncada = {Cp_H_S_trunc}')

            S_H_S_trunc = entropia(df_func_part_H_S_trunc, Temp,
                                  funcao_particao_Heibbe_Scalabrini_truncada_sympy(M, T,
                                  pressao, we, wexe, weye, Be, alfa_e, gama_e, gel, de,
                                  nu).evalf(subs={T: Temp}))
            dados_entropia_H_S_trunc.append(S_H_S_trunc)
            print(f'Entropia Heibbe-Scalabrini Truncada = {S_H_S_trunc}')

            G_H_S_trunc = energia_Gibbs(H_H_S_trunc, S_H_S_trunc, Temp)
            dados_G_H_S_trunc.append(G_H_S_trunc)
            print(f'Energia de Gibbs = {G_H_S_trunc}')


            print('-'*60)
            print('\n')


            ### Heibbe - Scalabrino Rotor rígido

            func_part_H_S_rot_rig = funcao_particao_Scalabrini_Rotor_Rigido_sympy(M, T,
                                    pressao, we, wexe, weye, gel, de, nu, theta_rot)

            #pprint(func_part_H_S_rot_rig)
            #print('\n')
            df_func_part_H_S_rot_rig_nao_aval = diff(func_part_H_S_rot_rig, T)
            df_func_part_H_S_rot_rig = diff(func_part_H_S_rot_rig, T).evalf(subs={T: Temp})
            print(f'Derivada Heibbe-Scalabrini Rotor Rigido = {df_func_part_H_S_rot_rig}')

            U_H_S_rot_rig_nao_aval = energia_interna(df_func_part_H_S_rot_rig_nao_aval, T)
            U_H_S_rot_rig = energia_interna(df_func_part_H_S_rot_rig, Temp)
            print(f'Energia interna Heibbe-Scalabrini Rotor Rigido = {U_H_S_rot_rig}')

            if Temp == Temp_inicial:
               H_H_S_rot_rig_298 = entalpia(U_H_S_rot_rig, Temp_inicial)
               H_H_S_rot_rig = entalpia(U_H_S_rot_rig, Temp)
               dados_entalpia_H_S_rot_rig.append(H_H_S_rot_rig - H_H_S_rot_rig_298)
               print(f'Entalpia Heibbe-Scalabrini Rotor Rigido = {H_H_S_rot_rig}')
            else:
               H_H_S_rot_rig = entalpia(U_H_S_rot_rig, Temp)
               dados_entalpia_H_S_rot_rig.append(H_H_S_rot_rig - H_H_S_rot_rig_298)
               print(f'Entalpia Heibbe-Scalabrini Rotor Rigido = {H_H_S_rot_rig}')

            Cp_H_S_rot_rid = diff(entalpia(U_H_S_rot_rig_nao_aval, T)).evalf(subs={T: Temp})
            dados_cp_H_S_rot_rig.append(Cp_H_S_rot_rid)
            print(f'Cp Heibbe-Scalabrini Rotor Rigido = {Cp_H_S_rot_rid}')

            S_H_S_rot_rig = entropia(df_func_part_H_S_rot_rig, Temp,
                                     funcao_particao_Scalabrini_Rotor_Rigido_sympy(M, T,
                                     pressao, we, wexe, weye, gel, de, nu,
                                     theta_rot).evalf(subs={T: Temp}))
            dados_entropia_H_S_rot_rig.append(S_H_S_rot_rig)
            print(f'Entropia Heibbe-Scalabrini Rotor Rigido = {S_H_S_rot_rig}')

            G_H_S_rot_rig = energia_Gibbs(H_H_S_rot_rig, S_H_S_rot_rig, Temp)
            dados_G_H_S_rot_rig.append(G_H_S_rot_rig)
            print(f'Energia de Gibbs = {G_H_S_rot_rig}')

            print('-'*60)
            print('\n')

            print(f'{Temp},{U_Mcquarie},{U_Allison_harm},{U_Allison},{U_Foglia},{U_H_S},{U_H_S_trunc},{U_H_S_rot_rig}', end='\n', file=f1)
            print(f'{Temp},{H_Mcquarie},{H_Allison_harm},{H_Allison},{H_Foglia},{H_H_S},{H_H_S_trunc},{H_H_S_rot_rig}', end='\n', file=f2)
            print(f'{Temp},{S_Mcquarie},{S_Allison_harm},{S_Allison},{S_Foglia},{S_H_S},{S_H_S_trunc},{S_H_S_rot_rig}', end='\n', file=f3)
            print(f'{Temp},{Cp_Mcquarie},{Cp_Allison_harm},{Cp_Allison},{Cp_Foglia},{Cp_H_S},{Cp_H_S_trunc},{Cp_H_S_rot_rid}', end='\n', file=f4)
            print(f'{Temp},{G_Mcquarie},{G_Allison_harm},{G_Allison},{G_Foglia},{G_H_S},{G_H_S_trunc},{G_H_S_rot_rig}', end='\n', file=f5)
            '''


    plt.plot(faixa_Temp, dados_entalpia_Macquarie, label=r'$\Delta$Macquarie')
    plt.plot(faixa_Temp, dados_entalpia_Allison_harm, label='Allison - Harm')
    plt.xlabel(r"$Temperatura$ $[K]$")
    plt.ylabel(r"$Entalpia$ $(J/mol)$")
    plt.legend()
    plt.show()

    '''
    plt.plot(faixa_Temp, dados_entalpia_Allison, label='Allison')
    plt.plot(faixa_Temp, dados_entalpia_Foglia, label='Foglia')
    plt.plot(faixa_Temp, dados_entalpia_H_S, label='Heibbe-Scalabrini')
    plt.plot(faixa_Temp, dados_entalpia_H_S_trunc, label='Heibbe-Scalabrini-Truc')
    plt.plot(faixa_Temp, dados_entalpia_H_S_rot_rig, label='Heibbe-Scalabrini-Rot-Rig')
    '''

    plt.plot(faixa_Temp, dados_entropia_Macquarie, label='Macquarie')
    plt.plot(faixa_Temp, dados_entropia_Allison_harm, label='Allison - Harm')
    plt.xlabel(r"$Temperatura [K]$")
    plt.ylabel(r"$Entropia (J/K mol)$")
    plt.legend()
    plt.show()

    '''
    plt.plot(faixa_Temp, dados_entropia_Allison, label='Allison')
    plt.plot(faixa_Temp, dados_entropia_Foglia, label='Foglia')
    plt.plot(faixa_Temp, dados_entropia_H_S, label='Heibbe-Scalabrini')
    plt.plot(faixa_Temp, dados_entropia_H_S_trunc, label='Heibbe-Scalabrini-Truc')
    plt.plot(faixa_Temp, dados_entropia_H_S_rot_rig, label='Scalabrini-Rot-Rig')
    '''

    plt.plot(faixa_Temp, dados_cp_Mcquarie, label='Macquarie')
    plt.plot(faixa_Temp, dados_cp_Allison_harm, label='Allison - Harm')
    plt.xlabel(r"$Temperatura [K]$")
    plt.ylabel(r"$Cp$")
    plt.legend()
    plt.show()

    '''
    plt.plot(faixa_Temp, dados_cp_Allison, label='Allison')
    plt.plot(faixa_Temp, dados_cp_Foglia, label='Foglia')
    plt.plot(faixa_Temp, dados_cp_H_S, label='Heibbe-Scalabrini')
    plt.plot(faixa_Temp, dados_cp_H_S_trunc, label='Heibbe-Scalabrini-Truc')
    plt.plot(faixa_Temp, dados_cp_H_S_rot_rig, label='Scalabrini-Rot-Rig')
    '''

    plt.plot(faixa_Temp, dados_G_Mcquarie, label='Macquarie')
    plt.plot(faixa_Temp, dados_G_Allison_harm, label='Allison Harm.')
    plt.xlabel(r"$Temperatura [K]$")
    plt.ylabel(r"$Energia Gibbs$")
    plt.legend()
    plt.show()

    '''
    plt.plot(faixa_Temp, dados_G_Allison, label='Allison')
    plt.plot(faixa_Temp, dados_G_Foglia, label='Foglia')
    plt.plot(faixa_Temp, dados_G_H_S, label='Heibbe-Scalabrini')
    plt.plot(faixa_Temp, dados_G_H_S_trunc, label='Heibbe-Scalabrini-Truc')
    plt.plot(faixa_Temp, dados_G_H_S_rot_rig, label='Heibbe-Scalabrini-Rot-Rig')

    '''
