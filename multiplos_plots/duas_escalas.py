import matplotlib.pyplot as plt
import numpy as np

fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(8, 6), sharex=True)

fig.subplots_adjust(
                     left=0.08,
                     right=0.86,         # {Define as distâncias entre os extremos}
                     bottom=0.08,
                     top=.94,
                     hspace=0.33,   # Organiza espaçoes entre os subplots
                     wspace=0.27    # Organiza espaçoes entre os subplots
                   )

titulo1 = ['Linear', 'Quadrado', 'Cubo', 'Quarta']
titulo2 = ['Seno', 'Cosseno', 'Exponencial', 'Logaritmo']

x = np.linspace(2, 5, 20)
y_vals = [x, x*x, x**3, x**4]
y_novo = [np.sin(x), np.cos(x), np.exp(x), np.log(x)]

for ax, titulo, y, yy in zip(axes.flatten(), titulo1, y_vals, y_novo):
    axx = ax.twinx()
    ax.plot(x, y)
    axx.plot(x, yy, color='r')
    ax.set_title(titulo)

plt.show()
