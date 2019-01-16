import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

x = np.linspace(-6.5, 6.5, 100)
y = np.exp(-0.2 * (x-0.5)**2)

fonte = {
         'family': 'serif',
         'color':  'k',
         'weight': 'normal',
         'size': 16,
        }

sns.set(style="ticks")
plt.plot(x, y, linewidth=2, label='Gaussiana')
plt.legend(loc='upper left', fontsize='large', shadow=False, frameon=False)
plt.xlabel('x', fontdict=fonte)
plt.ylabel('g(x)', fontdict=fonte)

plt.savefig(
             'gaussiana.png',
             dpi=300,
             orientation = 'portrait',
             transparent = True,
             format='png'
            )

plt.show()
