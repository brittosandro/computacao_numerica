import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

dados = np.loadtxt('He_c1.dat')
x = dados[:, 0]
y = dados[:, 9]

dados1 = np.loadtxt('He_c2.dat')
x1 = dados1[:, 0]
y1 = dados1[:, 9]

dados2 = np.loadtxt('Ne_c1.dat')
x2 = dados2[:, 0]
y2 = dados2[:, 9]

dados3 = np.loadtxt('Ne_c2.dat')
x3 = dados3[:, 0]
y3 = dados3[:, 9]

dados4 = np.loadtxt('Ar_c1.dat')
x4 = dados4[:, 0]
y4 = dados4[:, 9]

dados5 = np.loadtxt('Ar_c2.dat')
x5 = dados5[:, 0]
y5 = dados5[:, 9]

dados6 = np.loadtxt('Kr_c1.dat')
x6 = dados6[:, 0]
y6 = dados6[:, 9]

dados7 = np.loadtxt('Kr_c2.dat')
x7 = dados7[:, 0]
y7 = dados7[:, 9]


sns.set_style()
fig = plt.figure(figsize=(8, 18))


sub1 = fig.add_subplot(4, 2, 1)
sub1.plot(x, y*(43.3641))
sub1.set_title('Metanol + He (Caso 1)', fontsize = 'large', fontweight = 'bold')

sub2 = fig.add_subplot(4, 2, 2)
sub2.plot(x1, y1*(43.3641))
sub2.set_title('Metanol + He (Caso 2)', fontsize = 'large', fontweight = 'bold')


sub3 = fig.add_subplot(4, 2, 3)
sub3.plot(x2, y2*(43.3641))
sub3.set_title('Metanol + Ne (Caso 1)', fontsize = 'large', fontweight = 'bold')


sub4 = fig.add_subplot(4, 2, 4)
sub4.plot(x3, y3*(43.3641))
sub4.set_title('Metanol + Ne (Caso 2)', fontsize = 'large', fontweight = 'bold')


sub5 = fig.add_subplot(4, 2, 5)
sub5.set_title('Metanol + Ar (Caso 1)', fontsize = 'large', fontweight = 'bold')
sub5.plot(x4, y4*(43.3641))

sub6 = fig.add_subplot(4, 2, 6)
sub6.set_title('Metanol + Ar (Caso 2)', fontsize = 'large', fontweight = 'bold')
sub6.plot(x5, y5*(43.3641))

sub7 = fig.add_subplot(4, 2, 7)
sub7.set_title('Metanol + Kr (Caso 1)', fontsize = 'large', fontweight = 'bold')
sub7.plot(x6, y6*(43.3641))

sub8 = fig.add_subplot(4, 2, 8)
sub8.set_title('Metanol + Kr (Caso 2)', fontsize = 'large', fontweight = 'bold')
sub8.plot(x7, y7*(43.3641))

#fig.tight_layout()

plt.show()
