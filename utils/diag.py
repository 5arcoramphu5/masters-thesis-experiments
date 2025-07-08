from sympy import *

def print_complex(x):
    print('Complex(', re(x), ',', im(x), ')', end=' ')

def print_matrix(label, M):
    print(label)
    for i in range(4):
        print('{', end='')
        for j in range(4):
            print_complex(M[i, j])
            print(', ', end=' ')
        print('},')


M = Matrix([[0, 0,  1, 0],
            [0, 0, 0, 1],
            [0.75, 0, 0, 2],
            [0, 2.25, -2, 0]])

orderChange = Matrix([
    [0, 1, 0, 0],
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [1, 0, 0, 0]
])

P, D = M.diagonalize()  
P = P * orderChange

print_matrix("D:", D)
print_matrix("", N(P.inv() * M * P))
print_matrix("P:", P)
print_matrix("P inv:", P.inv())
print_matrix("P inverse * P:", N(P.inv() * P))