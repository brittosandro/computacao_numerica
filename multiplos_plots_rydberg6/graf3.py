import matplotlib.pyplot as plt
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

dados2 = np.loadtxt('Ne_c1.dat')
x2 = dados2[:, 0]
y2 = dados2[:, 9]

dados3 = np.loadtxt('Ne_c2.dat')
x3 = dados3[:, 0]
y3 = dados3[:, 9]

dados33 = np.loadtxt('Ne_c3.dat')
x33 = dados33[:, 0]
y33 = dados33[:, 9]

dados4 = np.loadtxt('Ar_c1.dat')
x4 = dados4[:, 0]
y4 = dados4[:, 9]

dados5 = np.loadtxt('Ar_c2.dat')
x5 = dados5[:, 0]
y5 = dados5[:, 9]

dados55 = np.loadtxt('Ar_c3.dat')
x55 = dados55[:, 0]
y55 = dados55[:, 9]

dados6 = np.loadtxt('Kr_c1.dat')
x6 = dados6[:, 0]
y6 = dados6[:, 9]

dados7 = np.loadtxt('Kr_c2.dat')
x7 = dados7[:, 0]
y7 = dados7[:, 9]

dados77 = np.loadtxt('Kr_c3.dat')
x77 = dados77[:, 0]
y77 = dados77[:, 9]

#sns.set_style()

fig, ax = plt.subplots(4, 3, constrained_layout=True)

fig.suptitle(
               'Demonstração das CEPs para os diferentes casos',
               fontweight = 'bold',
               fontsize = 14
            )

# ------------------- Subplot He Caso 1 --------------------------------------
ax[0, 0].plot(x, y*(43.3641), '#8A2BE2', lw=2)
ax[0, 0].set_title('Metanol + He (Caso 1)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot He Caso 2 --------------------------------------
ax[0, 1].plot(x1, y1*(43.3641), lw=2)
ax[0, 1].set_title('Metanol + He (Caso 2)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot He Caso 3 --------------------------------------
ax[0, 2].plot(x11, y11*(43.3641), '#006400', lw=2)
ax[0, 2].set_title('Metanol + He (Caso 3)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ne Caso 1 --------------------------------------
ax[1, 0].plot(x2, y2*(43.3641), '#8A2BE2', lw=2)
ax[1, 0].set_title('Metanol + Ne (Caso 1)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ne Caso 2 --------------------------------------
ax[1, 1].plot(x3, y3*(43.3641), lw=2)
ax[1, 1].set_title('Metanol + Ne (Caso 2)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ne Caso 3 --------------------------------------
ax[1, 2].plot(x33, y33*(43.3641), '#006400', lw=2)
ax[1, 2].set_title('Metanol + Ne (Caso 3)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ar Caso 1 --------------------------------------
ax[2, 0].plot(x4, y4*(43.3641), '#8A2BE2', lw=2)
ax[2, 0].set_title('Metanol + Ar (Caso 1)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ar Caso 2 --------------------------------------
ax[2, 1].plot(x5, y5*(43.3641), lw=2)
ax[2, 1].set_title('Metanol + Ar (Caso 2)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Ar Caso 3 --------------------------------------
ax[2, 2].plot(x55, y55*(43.3641), '#006400', lw=2)
ax[2, 2].set_title('Metanol + Ar (Caso 3)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Kr Caso 1 --------------------------------------
ax[3, 0].plot(x6, y6*(43.3641), '#8A2BE2')
ax[3, 0].set_title('Metanol + Kr (Caso 1)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Kr Caso 2 --------------------------------------
ax[3, 1].plot(x7, y7*(43.3641), lw=2)
ax[3, 1].set_title('Metanol + Kr (Caso 2)', fontsize = 'small', fontweight = 'bold')

# ------------------- Subplot Kr Caso 3 --------------------------------------
ax[3, 2].plot(x77, y77*(43.3641), '#006400', lw=2)
ax[3, 2].set_title('Metanol + Kr (Caso 3)', fontsize = 'small', fontweight = 'bold')

plt.savefig('ryd6_teste1.png', dpi=300, orientation='portrait', transparent=True, format='png')
plt.show()
