import matplotlib.pyplot as plt
import pandas as pd

dados_U = pd.read_csv('energia_interna.csv',)
dados_H = pd.read_csv('entalpia.csv',)
dados_S = pd.read_csv('entropia.csv',)

'''
print('Energia Interna U (J/mol)')
print(dados_U.head())
print()

print('Entalpia (J/mol)')
print(dados_H.head())
print()

print('Entropia (J/k mol)')
print(dados_S.head())
print()
'''

H_298_Mcquarie = dados_H['H_Mcquarie'][0]
delta_H_Mcquarie = dados_H['H_Mcquarie'] - H_298_Mcquarie

H_298_Allison_harm = dados_H['H_Allison_Harmonica'][0]
delta_H_Allison_harm = dados_H['H_Allison_Harmonica'] - H_298_Allison_harm

H_298_Allison = dados_H['H_Allison'][0]
delta_H_Allison = dados_H['H_Allison'] - H_298_Allison

H_298_Foglia = dados_H['H_Foglia'][0]
delta_H_Foglia = dados_H['H_Foglia'] - H_298_Foglia

H_298_H_H_S = dados['H_Heibbe_Scalabrini'][0]
delta_H_H_S = dados['H_Heibbe_Scalabrini'] - H_298_H_H_S

H_298_H_H_S_truc = dados['H_Heibbe_Scalabrini_Trunc'][0]
delta_H_H_S = dados['H_Heibbe_Scalabrini_Trunc'] - H_298_H_H_S_truc

H_298_H_H_S_rot_rig = dados['H_Heibbe_Scalabrini_ROTOR-RIG'][0]
delta_H_H_S = dados['H_Heibbe_Scalabrini_ROTOR-RIG'] - H_298_H_H_S_rot_rig



plt.plot(dados_H['Temperatura'], delta_H_Allison_harm)
plt.show()
