import time

def fatorial(n):
    return 1 if n < 2 else n*fatorial(n-1)

if __name__ == "__main__":
    n = 300
    t1 = time.time()
    for i in range(0, n+1):
        print(fatorial(i))
    t2 = time.time() - t1

    print('\n##########################################################')
    print('#O tempo que o Python levou para calcular foi de:        #')
    print('#{:3f}                                                #'.format(t2))
    print('##########################################################')
