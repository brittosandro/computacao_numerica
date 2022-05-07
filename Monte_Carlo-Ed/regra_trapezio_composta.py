import numpy as np


def func1(x):
    return (x**3) * (np.log(x))

def regra_trapezio_comp(x, y):
    rtrap = []
    n = len(x)
    for i in range(n):
        if i == 0 or i == n-1:
            rtrap.append(y[i])
        else:
            rtrap.append(2*y[i])
    I_1 = (h/2)*sum(rtrap)

    return I_1


intervalo_ini = 1
intervalo_final = 3
sub_divisoes = 4
h = (intervalo_final - intervalo_ini) / sub_divisoes

x = [i for i in np.arange(intervalo_ini, intervalo_final+0.1, h)]
y = [round(func1(i), 4) for i in x]

print(x)
print(y)
I_1 = regra_trapezio_comp(x, y)


print(I_1)
