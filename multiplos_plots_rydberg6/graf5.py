import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
#import seaborn as sns

dados = np.loadtxt('He_c1.dat')
x = dados[:, 0]
y = dados[:, 9]

dados1 = np.loadtxt('He_c2.dat')
x1 = dados1[:, 0]
y1 = dados1[:, 9]

dados11 = np.loadtxt('He_c3.dat')
x11 = dados11[:, 0]
y11 = dados11[:, 9]

# ------------------niveis vibracionais He ----------------------------------
dados_vib1 = np.loadtxt('He_c1_v0.txt')
vib1_x = dados_vib1[:, 0]
vib1_y = dados_vib1[:, 1]

dados_vib_c2_0 = np.loadtxt('He_c2_v0.txt')
vib_x_c2_v0 = dados_vib_c2_0[:, 0]
vib_y_c2_v0 = dados_vib_c2_0[:, 1]

dados_vib_c2_1 = np.loadtxt('He_c2_v1.txt')
vib_x_c2_v1 = dados_vib_c2_1[:, 0]
vib_y_c2_v1 = dados_vib_c2_1[:, 1]

dados_vib_c3_0 = np.loadtxt('He_c3_v0.txt')
vib_x_c3_v0 = dados_vib_c3_0[:, 0]
vib_y_c3_v0 = dados_vib_c3_0[:, 1]

dados_vib_c3_1 = np.loadtxt('He_c3_v1.txt')
vib_x_c3_v1 = dados_vib_c3_1[:, 0]
vib_y_c3_v1 = dados_vib_c3_1[:, 1]
# ----------------------------------------------------------------------------


dados2 = np.loadtxt('Ne_c1.dat')
x2 = dados2[:, 0]
y2 = dados2[:, 9]

dados3 = np.loadtxt('Ne_c2.dat')
x3 = dados3[:, 0]
y3 = dados3[:, 9]

dados33 = np.loadtxt('Ne_c3.dat')
x33 = dados33[:, 0]
y33 = dados33[:, 9]

# ------------------niveis vibracionais Ne -----------------------------------
# ****** Caso 1
dados_Ne_vib_c1_v0 = np.loadtxt('Ne_c1_v0.txt')
vib_x_Ne_c1_v0 = dados_Ne_vib_c1_v0[:, 0]
vib_y_Ne_c1_v0 = dados_Ne_vib_c1_v0[:, 1]

dados_Ne_vib_c1_v1 = np.loadtxt('Ne_c1_v1.txt')
vib_x_Ne_c1_v1 = dados_Ne_vib_c1_v1[:, 0]
vib_y_Ne_c1_v1 = dados_Ne_vib_c1_v1[:, 1]

dados_Ne_vib_c1_v2 = np.loadtxt('Ne_c1_v2.txt')
vib_x_Ne_c1_v2 = dados_Ne_vib_c1_v2[:, 0]
vib_y_Ne_c1_v2 = dados_Ne_vib_c1_v2[:, 1]

dados_Ne_vib_c1_v3 = np.loadtxt('Ne_c1_v3.txt')
vib_x_Ne_c1_v3 = dados_Ne_vib_c1_v3[:, 0]
vib_y_Ne_c1_v3 = dados_Ne_vib_c1_v3[:, 1]
# ****** fim Caso 1

# ****** Caso 2
dados_Ne_vib_c2_v0 = np.loadtxt('Ne_c2_v0.txt')
vib_x_Ne_c2_v0 = dados_Ne_vib_c2_v0[:, 0]
vib_y_Ne_c2_v0 = dados_Ne_vib_c2_v0[:, 1]

dados_Ne_vib_c2_v1 = np.loadtxt('Ne_c2_v1.txt')
vib_x_Ne_c2_v1 = dados_Ne_vib_c2_v1[:, 0]
vib_y_Ne_c2_v1 = dados_Ne_vib_c2_v1[:, 1]

dados_Ne_vib_c2_v2 = np.loadtxt('Ne_c2_v2.txt')
vib_x_Ne_c2_v2 = dados_Ne_vib_c2_v2[:, 0]
vib_y_Ne_c2_v2 = dados_Ne_vib_c2_v2[:, 1]

dados_Ne_vib_c2_v3 = np.loadtxt('Ne_c2_v3.txt')
vib_x_Ne_c2_v3 = dados_Ne_vib_c2_v3[:, 0]
vib_y_Ne_c2_v3 = dados_Ne_vib_c2_v3[:, 1]

dados_Ne_vib_c2_v4 = np.loadtxt('Ne_c2_v4.txt')
vib_x_Ne_c2_v4 = dados_Ne_vib_c2_v4[:, 0]
vib_y_Ne_c2_v4 = dados_Ne_vib_c2_v4[:, 1]
# ****** fim Caso 2

# ****** Caso 3
dados_Ne_vib_c3_v0 = np.loadtxt('Ne_c3_v0.txt')
vib_x_Ne_c3_v0 = dados_Ne_vib_c3_v0[:, 0]
vib_y_Ne_c3_v0 = dados_Ne_vib_c3_v0[:, 1]

dados_Ne_vib_c3_v1 = np.loadtxt('Ne_c3_v1.txt')
vib_x_Ne_c3_v1 = dados_Ne_vib_c3_v1[:, 0]
vib_y_Ne_c3_v1 = dados_Ne_vib_c3_v1[:, 1]

dados_Ne_vib_c3_v2 = np.loadtxt('Ne_c3_v2.txt')
vib_x_Ne_c3_v2 = dados_Ne_vib_c3_v2[:, 0]
vib_y_Ne_c3_v2 = dados_Ne_vib_c3_v2[:, 1]

dados_Ne_vib_c3_v3 = np.loadtxt('Ne_c3_v3.txt')
vib_x_Ne_c3_v3 = dados_Ne_vib_c3_v3[:, 0]
vib_y_Ne_c3_v3 = dados_Ne_vib_c3_v3[:, 1]

dados_Ne_vib_c3_v4 = np.loadtxt('Ne_c3_v4.txt')
vib_x_Ne_c3_v4 = dados_Ne_vib_c3_v4[:, 0]
vib_y_Ne_c3_v4 = dados_Ne_vib_c3_v4[:, 1]

# ****** fim Caso 3
# ----------------------------------------------------------------------------

dados4 = np.loadtxt('Ar_c1.dat')
x4 = dados4[:, 0]
y4 = dados4[:, 9]

dados5 = np.loadtxt('Ar_c2.dat')
x5 = dados5[:, 0]
y5 = dados5[:, 9]

dados55 = np.loadtxt('Ar_c3.dat')
x55 = dados55[:, 0]
y55 = dados55[:, 9]

# ------------------niveis vibracionais Ar -----------------------------------
# ****** Caso 1
dados_Ar_vib_c1_v0 = np.loadtxt('Ar_c1_v0.txt')
vib_x_Ar_c1_v0 = dados_Ar_vib_c1_v0[:, 0]
vib_y_Ar_c1_v0 = dados_Ar_vib_c1_v0[:, 1]

dados_Ar_vib_c1_v1 = np.loadtxt('Ar_c1_v1.txt')
vib_x_Ar_c1_v1 = dados_Ar_vib_c1_v1[:, 0]
vib_y_Ar_c1_v1 = dados_Ar_vib_c1_v1[:, 1]

dados_Ar_vib_c1_v2 = np.loadtxt('Ar_c1_v2.txt')
vib_x_Ar_c1_v2 = dados_Ar_vib_c1_v2[:, 0]
vib_y_Ar_c1_v2 = dados_Ar_vib_c1_v2[:, 1]

dados_Ar_vib_c1_v3 = np.loadtxt('Ar_c1_v3.txt')
vib_x_Ar_c1_v3 = dados_Ar_vib_c1_v3[:, 0]
vib_y_Ar_c1_v3 = dados_Ar_vib_c1_v3[:, 1]

dados_Ar_vib_c1_v4 = np.loadtxt('Ar_c1_v4.txt')
vib_x_Ar_c1_v4 = dados_Ar_vib_c1_v4[:, 0]
vib_y_Ar_c1_v4 = dados_Ar_vib_c1_v4[:, 1]

dados_Ar_vib_c1_v5 = np.loadtxt('Ar_c1_v5.txt')
vib_x_Ar_c1_v5 = dados_Ar_vib_c1_v5[:, 0]
vib_y_Ar_c1_v5 = dados_Ar_vib_c1_v5[:, 1]

dados_Ar_vib_c1_v6 = np.loadtxt('Ar_c1_v6.txt')
vib_x_Ar_c1_v6 = dados_Ar_vib_c1_v6[:, 0]
vib_y_Ar_c1_v6 = dados_Ar_vib_c1_v6[:, 1]

dados_Ar_vib_c1_v7 = np.loadtxt('Ar_c1_v7.txt')
vib_x_Ar_c1_v7 = dados_Ar_vib_c1_v7[:, 0]
vib_y_Ar_c1_v7 = dados_Ar_vib_c1_v7[:, 1]

dados_Ar_vib_c1_v8 = np.loadtxt('Ar_c1_v8.txt')
vib_x_Ar_c1_v8 = dados_Ar_vib_c1_v8[:, 0]
vib_y_Ar_c1_v8 = dados_Ar_vib_c1_v8[:, 1]
# ****** fim Caso 1

# ****** Caso 2
dados_Ar_vib_c2_v0 = np.loadtxt('Ar_c2_v0.txt')
vib_x_Ar_c2_v0 = dados_Ar_vib_c2_v0[:, 0]
vib_y_Ar_c2_v0 = dados_Ar_vib_c2_v0[:, 1]

dados_Ar_vib_c2_v1 = np.loadtxt('Ar_c2_v1.txt')
vib_x_Ar_c2_v1 = dados_Ar_vib_c2_v1[:, 0]
vib_y_Ar_c2_v1 = dados_Ar_vib_c2_v1[:, 1]

dados_Ar_vib_c2_v2 = np.loadtxt('Ar_c2_v2.txt')
vib_x_Ar_c2_v2 = dados_Ar_vib_c2_v2[:, 0]
vib_y_Ar_c2_v2 = dados_Ar_vib_c2_v2[:, 1]

dados_Ar_vib_c2_v3 = np.loadtxt('Ar_c2_v3.txt')
vib_x_Ar_c2_v3 = dados_Ar_vib_c2_v3[:, 0]
vib_y_Ar_c2_v3 = dados_Ar_vib_c2_v3[:, 1]

dados_Ar_vib_c2_v4 = np.loadtxt('Ar_c2_v4.txt')
vib_x_Ar_c2_v4 = dados_Ar_vib_c2_v4[:, 0]
vib_y_Ar_c2_v4 = dados_Ar_vib_c2_v4[:, 1]

dados_Ar_vib_c2_v5 = np.loadtxt('Ar_c2_v5.txt')
vib_x_Ar_c2_v5 = dados_Ar_vib_c2_v5[:, 0]
vib_y_Ar_c2_v5 = dados_Ar_vib_c2_v5[:, 1]

dados_Ar_vib_c2_v6 = np.loadtxt('Ar_c2_v6.txt')
vib_x_Ar_c2_v6 = dados_Ar_vib_c2_v6[:, 0]
vib_y_Ar_c2_v6 = dados_Ar_vib_c2_v6[:, 1]

dados_Ar_vib_c2_v7 = np.loadtxt('Ar_c2_v7.txt')
vib_x_Ar_c2_v7 = dados_Ar_vib_c2_v7[:, 0]
vib_y_Ar_c2_v7 = dados_Ar_vib_c2_v7[:, 1]

dados_Ar_vib_c2_v8 = np.loadtxt('Ar_c2_v8.txt')
vib_x_Ar_c2_v8 = dados_Ar_vib_c2_v8[:, 0]
vib_y_Ar_c2_v8 = dados_Ar_vib_c2_v8[:, 1]

dados_Ar_vib_c2_v9 = np.loadtxt('Ar_c2_v9.txt')
vib_x_Ar_c2_v9 = dados_Ar_vib_c2_v9[:, 0]
vib_y_Ar_c2_v9 = dados_Ar_vib_c2_v9[:, 1]

dados_Ar_vib_c2_v10 = np.loadtxt('Ar_c2_v10.txt')
vib_x_Ar_c2_v10 = dados_Ar_vib_c2_v10[:, 0]
vib_y_Ar_c2_v10 = dados_Ar_vib_c2_v10[:, 1]
# ****** fim Caso 2

# ****** Caso 3
dados_Ar_vib_c3_v0 = np.loadtxt('Ar_c3_v0.txt')
vib_x_Ar_c3_v0 = dados_Ar_vib_c3_v0[:, 0]
vib_y_Ar_c3_v0 = dados_Ar_vib_c3_v0[:, 1]

dados_Ar_vib_c3_v1 = np.loadtxt('Ar_c3_v1.txt')
vib_x_Ar_c3_v1 = dados_Ar_vib_c3_v1[:, 0]
vib_y_Ar_c3_v1 = dados_Ar_vib_c3_v1[:, 1]

dados_Ar_vib_c3_v2 = np.loadtxt('Ar_c3_v2.txt')
vib_x_Ar_c3_v2 = dados_Ar_vib_c3_v2[:, 0]
vib_y_Ar_c3_v2 = dados_Ar_vib_c3_v2[:, 1]

dados_Ar_vib_c3_v3 = np.loadtxt('Ar_c3_v3.txt')
vib_x_Ar_c3_v3 = dados_Ar_vib_c3_v3[:, 0]
vib_y_Ar_c3_v3 = dados_Ar_vib_c3_v3[:, 1]

dados_Ar_vib_c3_v4 = np.loadtxt('Ar_c3_v4.txt')
vib_x_Ar_c3_v4 = dados_Ar_vib_c3_v4[:, 0]
vib_y_Ar_c3_v4 = dados_Ar_vib_c3_v4[:, 1]

dados_Ar_vib_c3_v5 = np.loadtxt('Ar_c3_v5.txt')
vib_x_Ar_c3_v5 = dados_Ar_vib_c3_v5[:, 0]
vib_y_Ar_c3_v5 = dados_Ar_vib_c3_v5[:, 1]

dados_Ar_vib_c3_v6 = np.loadtxt('Ar_c3_v6.txt')
vib_x_Ar_c3_v6 = dados_Ar_vib_c3_v6[:, 0]
vib_y_Ar_c3_v6 = dados_Ar_vib_c3_v6[:, 1]

dados_Ar_vib_c3_v7 = np.loadtxt('Ar_c3_v7.txt')
vib_x_Ar_c3_v7 = dados_Ar_vib_c3_v7[:, 0]
vib_y_Ar_c3_v7 = dados_Ar_vib_c3_v7[:, 1]

dados_Ar_vib_c3_v8 = np.loadtxt('Ar_c3_v8.txt')
vib_x_Ar_c3_v8 = dados_Ar_vib_c3_v8[:, 0]
vib_y_Ar_c3_v8 = dados_Ar_vib_c3_v8[:, 1]

dados_Ar_vib_c3_v9 = np.loadtxt('Ar_c3_v9.txt')
vib_x_Ar_c3_v9 = dados_Ar_vib_c3_v9[:, 0]
vib_y_Ar_c3_v9 = dados_Ar_vib_c3_v9[:, 1]

dados_Ar_vib_c3_v10 = np.loadtxt('Ar_c3_v10.txt')
vib_x_Ar_c3_v10 = dados_Ar_vib_c3_v10[:, 0]
vib_y_Ar_c3_v10 = dados_Ar_vib_c3_v10[:, 1]

# ****** fim Caso 3
# ----------------------------------------------------------------------------

dados6 = np.loadtxt('Kr_c1.dat')
x6 = dados6[:, 0]
y6 = dados6[:, 9]

dados7 = np.loadtxt('Kr_c2.dat')
x7 = dados7[:, 0]
y7 = dados7[:, 9]

dados77 = np.loadtxt('Kr_c3.dat')
x77 = dados77[:, 0]
y77 = dados77[:, 9]

# ------------------niveis vibracionais Kr -----------------------------------
# ****** Caso 1
dados_Kr_vib_c1_v0 = np.loadtxt('Kr_c1_v0.txt')
vib_x_Kr_c1_v0 = dados_Kr_vib_c1_v0[:, 0]
vib_y_Kr_c1_v0 = dados_Kr_vib_c1_v0[:, 1]

dados_Kr_vib_c1_v1 = np.loadtxt('Kr_c1_v1.txt')
vib_x_Kr_c1_v1 = dados_Kr_vib_c1_v1[:, 0]
vib_y_Kr_c1_v1 = dados_Kr_vib_c1_v1[:, 1]

dados_Kr_vib_c1_v2 = np.loadtxt('Kr_c1_v2.txt')
vib_x_Kr_c1_v2 = dados_Kr_vib_c1_v2[:, 0]
vib_y_Kr_c1_v2 = dados_Kr_vib_c1_v2[:, 1]

dados_Kr_vib_c1_v3 = np.loadtxt('Kr_c1_v3.txt')
vib_x_Kr_c1_v3 = dados_Kr_vib_c1_v3[:, 0]
vib_y_Kr_c1_v3 = dados_Kr_vib_c1_v3[:, 1]

dados_Kr_vib_c1_v4 = np.loadtxt('Kr_c1_v4.txt')
vib_x_Kr_c1_v4 = dados_Kr_vib_c1_v4[:, 0]
vib_y_Kr_c1_v4 = dados_Kr_vib_c1_v4[:, 1]

dados_Kr_vib_c1_v4 = np.loadtxt('Kr_c1_v4.txt')
vib_x_Kr_c1_v4 = dados_Kr_vib_c1_v4[:, 0]
vib_y_Kr_c1_v4 = dados_Kr_vib_c1_v4[:, 1]

dados_Kr_vib_c1_v5 = np.loadtxt('Kr_c1_v5.txt')
vib_x_Kr_c1_v5 = dados_Kr_vib_c1_v5[:, 0]
vib_y_Kr_c1_v5 = dados_Kr_vib_c1_v5[:, 1]

dados_Kr_vib_c1_v6 = np.loadtxt('Kr_c1_v6.txt')
vib_x_Kr_c1_v6 = dados_Kr_vib_c1_v6[:, 0]
vib_y_Kr_c1_v6 = dados_Kr_vib_c1_v6[:, 1]

dados_Kr_vib_c1_v7 = np.loadtxt('Kr_c1_v7.txt')
vib_x_Kr_c1_v7 = dados_Kr_vib_c1_v7[:, 0]
vib_y_Kr_c1_v7 = dados_Kr_vib_c1_v7[:, 1]

dados_Kr_vib_c1_v8 = np.loadtxt('Kr_c1_v8.txt')
vib_x_Kr_c1_v8 = dados_Kr_vib_c1_v8[:, 0]
vib_y_Kr_c1_v8 = dados_Kr_vib_c1_v8[:, 1]

dados_Kr_vib_c1_v9 = np.loadtxt('Kr_c1_v9.txt')
vib_x_Kr_c1_v9 = dados_Kr_vib_c1_v9[:, 0]
vib_y_Kr_c1_v9 = dados_Kr_vib_c1_v9[:, 1]

dados_Kr_vib_c1_v10 = np.loadtxt('Kr_c1_v10.txt')
vib_x_Kr_c1_v10 = dados_Kr_vib_c1_v10[:, 0]
vib_y_Kr_c1_v10 = dados_Kr_vib_c1_v10[:, 1]

dados_Kr_vib_c1_v11 = np.loadtxt('Kr_c1_v11.txt')
vib_x_Kr_c1_v11 = dados_Kr_vib_c1_v11[:, 0]
vib_y_Kr_c1_v11 = dados_Kr_vib_c1_v11[:, 1]
# ****** Fim Caso 1

# ****** Caso 2
dados_Kr_vib_c2_v0 = np.loadtxt('Kr_c2_v0.txt')
vib_x_Kr_c2_v0 = dados_Kr_vib_c2_v0[:, 0]
vib_y_Kr_c2_v0 = dados_Kr_vib_c2_v0[:, 1]

dados_Kr_vib_c2_v1 = np.loadtxt('Kr_c2_v1.txt')
vib_x_Kr_c2_v1 = dados_Kr_vib_c2_v1[:, 0]
vib_y_Kr_c2_v1 = dados_Kr_vib_c2_v1[:, 1]

dados_Kr_vib_c2_v2 = np.loadtxt('Kr_c2_v2.txt')
vib_x_Kr_c2_v2 = dados_Kr_vib_c2_v2[:, 0]
vib_y_Kr_c2_v2 = dados_Kr_vib_c2_v2[:, 1]

dados_Kr_vib_c2_v3 = np.loadtxt('Kr_c2_v3.txt')
vib_x_Kr_c2_v3 = dados_Kr_vib_c2_v3[:, 0]
vib_y_Kr_c2_v3 = dados_Kr_vib_c2_v3[:, 1]

dados_Kr_vib_c2_v4 = np.loadtxt('Kr_c2_v4.txt')
vib_x_Kr_c2_v4 = dados_Kr_vib_c2_v4[:, 0]
vib_y_Kr_c2_v4 = dados_Kr_vib_c2_v4[:, 1]

dados_Kr_vib_c2_v5 = np.loadtxt('Kr_c2_v5.txt')
vib_x_Kr_c2_v5 = dados_Kr_vib_c2_v5[:, 0]
vib_y_Kr_c2_v5 = dados_Kr_vib_c2_v5[:, 1]

dados_Kr_vib_c2_v6 = np.loadtxt('Kr_c2_v6.txt')
vib_x_Kr_c2_v6 = dados_Kr_vib_c2_v6[:, 0]
vib_y_Kr_c2_v6 = dados_Kr_vib_c2_v6[:, 1]

dados_Kr_vib_c2_v7 = np.loadtxt('Kr_c2_v7.txt')
vib_x_Kr_c2_v7 = dados_Kr_vib_c2_v7[:, 0]
vib_y_Kr_c2_v7 = dados_Kr_vib_c2_v7[:, 1]

dados_Kr_vib_c2_v8 = np.loadtxt('Kr_c2_v8.txt')
vib_x_Kr_c2_v8 = dados_Kr_vib_c2_v8[:, 0]
vib_y_Kr_c2_v8 = dados_Kr_vib_c2_v8[:, 1]

dados_Kr_vib_c2_v9 = np.loadtxt('Kr_c2_v9.txt')
vib_x_Kr_c2_v9 = dados_Kr_vib_c2_v9[:, 0]
vib_y_Kr_c2_v9 = dados_Kr_vib_c2_v9[:, 1]

dados_Kr_vib_c2_v10 = np.loadtxt('Kr_c2_v10.txt')
vib_x_Kr_c2_v10 = dados_Kr_vib_c2_v10[:, 0]
vib_y_Kr_c2_v10 = dados_Kr_vib_c2_v10[:, 1]

dados_Kr_vib_c2_v11 = np.loadtxt('Kr_c2_v11.txt')
vib_x_Kr_c2_v11 = dados_Kr_vib_c2_v11[:, 0]
vib_y_Kr_c2_v11 = dados_Kr_vib_c2_v11[:, 1]
# ****** Fim Caso 2

# ****** Caso 3
dados_Kr_vib_c3_v0 = np.loadtxt('Kr_c3_v0.txt')
vib_x_Kr_c3_v0 = dados_Kr_vib_c3_v0[:, 0]
vib_y_Kr_c3_v0 = dados_Kr_vib_c3_v0[:, 1]

dados_Kr_vib_c3_v1 = np.loadtxt('Kr_c3_v1.txt')
vib_x_Kr_c3_v1 = dados_Kr_vib_c3_v1[:, 0]
vib_y_Kr_c3_v1 = dados_Kr_vib_c3_v1[:, 1]

dados_Kr_vib_c3_v2 = np.loadtxt('Kr_c3_v2.txt')
vib_x_Kr_c3_v2 = dados_Kr_vib_c3_v2[:, 0]
vib_y_Kr_c3_v2 = dados_Kr_vib_c3_v2[:, 1]

dados_Kr_vib_c3_v3 = np.loadtxt('Kr_c3_v3.txt')
vib_x_Kr_c3_v3 = dados_Kr_vib_c3_v3[:, 0]
vib_y_Kr_c3_v3 = dados_Kr_vib_c3_v3[:, 1]

dados_Kr_vib_c3_v4 = np.loadtxt('Kr_c3_v4.txt')
vib_x_Kr_c3_v4 = dados_Kr_vib_c3_v4[:, 0]
vib_y_Kr_c3_v4 = dados_Kr_vib_c3_v4[:, 1]

dados_Kr_vib_c3_v5 = np.loadtxt('Kr_c3_v5.txt')
vib_x_Kr_c3_v5 = dados_Kr_vib_c3_v5[:, 0]
vib_y_Kr_c3_v5 = dados_Kr_vib_c3_v5[:, 1]

dados_Kr_vib_c3_v6 = np.loadtxt('Kr_c3_v6.txt')
vib_x_Kr_c3_v6 = dados_Kr_vib_c3_v6[:, 0]
vib_y_Kr_c3_v6 = dados_Kr_vib_c3_v6[:, 1]

dados_Kr_vib_c3_v7 = np.loadtxt('Kr_c3_v7.txt')
vib_x_Kr_c3_v7 = dados_Kr_vib_c3_v7[:, 0]
vib_y_Kr_c3_v7 = dados_Kr_vib_c3_v7[:, 1]

dados_Kr_vib_c3_v8 = np.loadtxt('Kr_c3_v8.txt')
vib_x_Kr_c3_v8 = dados_Kr_vib_c3_v8[:, 0]
vib_y_Kr_c3_v8 = dados_Kr_vib_c3_v8[:, 1]

dados_Kr_vib_c3_v9 = np.loadtxt('Kr_c3_v9.txt')
vib_x_Kr_c3_v9 = dados_Kr_vib_c3_v9[:, 0]
vib_y_Kr_c3_v9 = dados_Kr_vib_c3_v9[:, 1]

dados_Kr_vib_c3_v10 = np.loadtxt('Kr_c3_v10.txt')
vib_x_Kr_c3_v10 = dados_Kr_vib_c3_v10[:, 0]
vib_y_Kr_c3_v10 = dados_Kr_vib_c3_v10[:, 1]

dados_Kr_vib_c3_v11 = np.loadtxt('Kr_c3_v11.txt')
vib_x_Kr_c3_v11 = dados_Kr_vib_c3_v11[:, 0]
vib_y_Kr_c3_v11 = dados_Kr_vib_c3_v11[:, 1]
# ****** Fim Caso 3


#sns.set_style()
fig, ax = plt.subplots(4, 3, figsize=(9.8, 8.5), constrained_layout=True)

'''
fig.suptitle(
               'Demonstração das CEPs para os diferentes casos',
               fontweight = 'bold',
               fontsize = 14
            )
'''
# ------------------- Subplot He Caso 1 --------------------------------------
ax[0, 0].plot(vib1_x, vib1_y, 'k')
ax[0, 0].plot(x, y*(43.3641), '#0000ff', lw=2)
ax[0, 0].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[0, 0].set_xlim(3, 7.)
ax[0, 0].set_ylim(-3.5, 2.5)
ax[0, 0].set_title('Metanol + He (Caso 1)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot He Caso 2 --------------------------------------
ax[0, 1].plot(vib_x_c2_v0, vib_y_c2_v0, 'k')
ax[0, 1].plot(vib_x_c2_v1, vib_y_c2_v1, 'k')
ax[0, 1].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[0, 1].set_xlim(1.9, 7.)
ax[0, 1].set_ylim(-5, 3)
ax[0, 1].plot(x1, y1*(43.3641), '#ff0000', lw=2)
ax[0, 1].set_title('Metanol + He (Caso 2)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot He Caso 3 -------------------------------------
ax[0, 2].plot(vib_x_c3_v0, vib_y_c3_v0, 'k')
ax[0, 2].plot(vib_x_c3_v1, vib_y_c3_v1, 'k')
ax[0, 2].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[0, 2].set_xlim(2.4, 7.)
ax[0, 2].set_ylim(-5.9, 5.)
ax[0, 2].plot(x11, y11*(43.3641), '#00cc00', lw=2)
ax[0, 2].set_title('Metanol + He (Caso 3)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ne Caso 1 --------------------------------------
ax[1, 0].plot(vib_x_Ne_c1_v0, vib_y_Ne_c1_v0, 'k')
ax[1, 0].plot(vib_x_Ne_c1_v1, vib_y_Ne_c1_v1, 'k')
ax[1, 0].plot(vib_x_Ne_c1_v2, vib_y_Ne_c1_v2, 'k')
ax[1, 0].plot(vib_x_Ne_c1_v3, vib_y_Ne_c1_v3, 'k')
ax[1, 0].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[1, 0].set_xlim(3., 7.)
ax[1, 0].set_ylim(-6.9, 4.7)
ax[1, 0].plot(x2, y2*(43.3641), '#0000ff', lw=2)
ax[1, 0].set_title('Metanol + Ne (Caso 1)', fontsize = 'small', fontweight = 'bold')
#ax[1, 0].set_ylabel('Energia (meV)', labelpad = 0, fontsize = 'xx-large')
plt.ylabel('Energia (meV)', labelpad = 0, fontsize = 'xx-large')

# ------------------- Subplot Ne Caso 2 --------------------------------------
ax[1, 1].plot(vib_x_Ne_c2_v0, vib_y_Ne_c2_v0, 'k')
ax[1, 1].plot(vib_x_Ne_c2_v1, vib_y_Ne_c2_v1, 'k')
ax[1, 1].plot(vib_x_Ne_c2_v2, vib_y_Ne_c2_v2, 'k')
ax[1, 1].plot(vib_x_Ne_c2_v3, vib_y_Ne_c2_v3, 'k')
ax[1, 1].plot(vib_x_Ne_c2_v4, vib_y_Ne_c2_v4, 'k')
ax[1, 1].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[1, 1].set_xlim(1.9, 7.)
ax[1, 1].set_ylim(-10.2, 5.3)
ax[1, 1].plot(x3, y3*(43.3641), '#ff0000', lw=2)
ax[1, 1].set_title('Metanol + Ne (Caso 2)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ne Caso 3 --------------------------------------
ax[1, 2].plot(vib_x_Ne_c3_v0, vib_y_Ne_c3_v0, 'k')
ax[1, 2].plot(vib_x_Ne_c3_v1, vib_y_Ne_c3_v1, 'k')
ax[1, 2].plot(vib_x_Ne_c3_v2, vib_y_Ne_c3_v2, 'k')
ax[1, 2].plot(vib_x_Ne_c3_v3, vib_y_Ne_c3_v3, 'k')
ax[1, 2].plot(vib_x_Ne_c3_v4, vib_y_Ne_c3_v4, 'k')
ax[1, 2].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[1, 2].set_xlim(2.5, 7.)
ax[1, 2].set_ylim(-10.2, 5.3)
ax[1, 2].plot(x33, y33*(43.3641), '#00cc00', lw=2)
ax[1, 2].set_title('Metanol + Ne (Caso 3)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ar Caso 1 --------------------------------------
ax[2, 0].plot(vib_x_Ar_c1_v0, vib_y_Ar_c1_v0, 'k')
ax[2, 0].plot(vib_x_Ar_c1_v1, vib_y_Ar_c1_v1, 'k')
ax[2, 0].plot(vib_x_Ar_c1_v2, vib_y_Ar_c1_v2, 'k')
ax[2, 0].plot(vib_x_Ar_c1_v3, vib_y_Ar_c1_v3, 'k')
ax[2, 0].plot(vib_x_Ar_c1_v4, vib_y_Ar_c1_v4, 'k')
ax[2, 0].plot(vib_x_Ar_c1_v5, vib_y_Ar_c1_v5, 'k')
ax[2, 0].plot(vib_x_Ar_c1_v6, vib_y_Ar_c1_v6, 'k')
ax[2, 0].plot(vib_x_Ar_c1_v7, vib_y_Ar_c1_v7, 'k')
ax[2, 0].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[2, 0].set_xlim(3.2, 7.)
ax[2, 0].set_ylim(-17.3, 10.3)
ax[2, 0].plot(x4, y4*(43.3641), '#0000ff', lw=2)
ax[2, 0].set_title('Metanol + Ar (Caso 1)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ar Caso 2 --------------------------------------
ax[2, 1].plot(vib_x_Ar_c2_v0, vib_y_Ar_c2_v0, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v1, vib_y_Ar_c2_v1, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v2, vib_y_Ar_c2_v2, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v3, vib_y_Ar_c2_v3, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v4, vib_y_Ar_c2_v4, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v5, vib_y_Ar_c2_v5, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v6, vib_y_Ar_c2_v6, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v7, vib_y_Ar_c2_v7, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v8, vib_y_Ar_c2_v8, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v9, vib_y_Ar_c2_v9, 'k')
ax[2, 1].plot(vib_x_Ar_c2_v10, vib_y_Ar_c2_v10, 'k')
ax[2, 1].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[2, 1].set_xlim(2., 7.)
ax[2, 1].set_ylim(-25.1, 10.3)
ax[2, 1].plot(x5, y5*(43.3641), '#ff0000', lw=2)
ax[2, 1].set_title('Metanol + Ar (Caso 2)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ar Caso 3 --------------------------------------
ax[2, 2].plot(vib_x_Ar_c3_v0, vib_y_Ar_c3_v0, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v1, vib_y_Ar_c3_v1, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v2, vib_y_Ar_c3_v2, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v3, vib_y_Ar_c3_v3, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v4, vib_y_Ar_c3_v4, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v5, vib_y_Ar_c3_v5, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v6, vib_y_Ar_c3_v6, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v7, vib_y_Ar_c3_v7, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v8, vib_y_Ar_c3_v8, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v9, vib_y_Ar_c3_v9, 'k')
ax[2, 2].plot(vib_x_Ar_c3_v10, vib_y_Ar_c3_v10, 'k')
ax[2, 2].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[2, 2].set_xlim(2.7, 7.)
ax[2, 2].set_ylim(-30., 10.3)
ax[2, 2].plot(x55, y55*(43.3641), '#00cc00', lw=2)
ax[2, 2].set_title('Metanol + Ar (Caso 3)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Kr Caso 1 --------------------------------------
ax[3, 0].plot(vib_x_Kr_c1_v0, vib_y_Kr_c1_v0, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v1, vib_y_Kr_c1_v1, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v2, vib_y_Kr_c1_v2, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v3, vib_y_Kr_c1_v3, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v4, vib_y_Kr_c1_v4, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v5, vib_y_Kr_c1_v5, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v6, vib_y_Kr_c1_v6, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v7, vib_y_Kr_c1_v7, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v8, vib_y_Kr_c1_v8, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v9, vib_y_Kr_c1_v9, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v10, vib_y_Kr_c1_v10, 'k')
ax[3, 0].plot(vib_x_Kr_c1_v11, vib_y_Kr_c1_v11, 'k')
ax[3, 0].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[3, 0].set_xlim(3.3, 7.)
ax[3, 0].set_ylim(-22., 12.9)
ax[3, 0].plot(x6, y6*(43.3641), '#0000ff', lw=2)
ax[3, 0].set_title('Metanol + Kr (Caso 1)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Kr Caso 2 --------------------------------------
ax[3, 1].plot(vib_x_Kr_c2_v0, vib_y_Kr_c2_v0, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v1, vib_y_Kr_c2_v1, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v2, vib_y_Kr_c2_v2, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v3, vib_y_Kr_c2_v3, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v4, vib_y_Kr_c2_v4, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v5, vib_y_Kr_c2_v5, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v6, vib_y_Kr_c2_v6, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v7, vib_y_Kr_c2_v7, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v8, vib_y_Kr_c2_v8, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v9, vib_y_Kr_c2_v9, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v10, vib_y_Kr_c2_v10, 'k')
ax[3, 1].plot(vib_x_Kr_c2_v11, vib_y_Kr_c2_v11, 'k')
ax[3, 1].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[3, 1].set_xlim(2.2, 7.)
ax[3, 1].set_ylim(-31.3, 10.5)
ax[3, 1].plot(x7, y7*(43.3641), '#ff0000', lw=2)
ax[3, 1].set_title('Metanol + Kr (Caso 2)', fontsize = 'small', fontweight = 'bold')

ax[3, 1].set_xlabel('Distância ($\AA$)', labelpad = 2, fontsize = 'xx-large')
# ------------------- Subplot Kr Caso 3 --------------------------------------
ax[3, 2].plot(vib_x_Kr_c3_v0, vib_y_Kr_c3_v0, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v1, vib_y_Kr_c3_v1, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v2, vib_y_Kr_c3_v2, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v3, vib_y_Kr_c3_v3, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v4, vib_y_Kr_c3_v4, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v5, vib_y_Kr_c3_v5, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v6, vib_y_Kr_c3_v6, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v7, vib_y_Kr_c3_v7, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v8, vib_y_Kr_c3_v8, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v9, vib_y_Kr_c3_v9, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v10, vib_y_Kr_c3_v10, 'k')
ax[3, 2].plot(vib_x_Kr_c3_v11, vib_y_Kr_c3_v11, 'k')
ax[3, 2].tick_params(labelcolor='k', labelsize='small', width=1.1)
ax[3, 2].set_xlim(2.9, 7.)
ax[3, 2].set_ylim(-32.1, 10.5)
ax[3, 2].plot(x77, y77*(43.3641), '#00cc00', lw=2)
ax[3, 2].set_title('Metanol + Kr (Caso 3)', fontsize = 'small', fontweight = 'bold')



plt.savefig('ryd6_teste5.png', dpi=300, orientation='portrait', transparent=True, format='png')
plt.show()
