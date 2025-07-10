from math import sqrt, prod, log2, log10, log, floor, ceil, nan

########
# Elliptic Functions
########

# Frenkel: A miracle describing an infinite sequence of numbers a(p) using Harmonic Analysis
#
# validate by computation: The counting problem solution for a particular elliptic
# curve modulo p equals the coefficients of a generating function: For each q raised
# to the p. The elliptic curve in consideration is y**2 + y = x**3 - x**2.
# 
# From numberphile: Frenkel on the Langlands program, https://www.youtube.com/watch?v=4dyytPboqvE&t=910s
#
# polynomial = q (1-q)**2 (1 - q**11)**2 (1-q**2)**2 (1-q**22)**2 (1-q**3)**2 (1-q**33)**2 (1-q**4)**2
#   giving coefficients of q, q**2, q**3, ... as 1, -2, -1, 2, 1, 2, -2, -2, -2, 1, -2, 4, ...
# This can be checked for primes 2, 3, 5, 7, 11, 13 (or more were we inclined to do some algebra)


def EvaluateEF(a0, a1, a2, b0, b1, b2, b3, x, y, p):
    '''Bool evaluate elliptic function mod p equality: quad in y (a) vs cubic in x (b)'''
    y_sum = a0 + a1*y + a2*(y**2) 
    x_sum = b0 + b1*x + b2*(x**2) + b3*(x**3)          # odd error: make final exponent 2 to recover count = -p
    return (y_sum % p) == (x_sum % p)

    
def EllipticTest(x, y, p):
    '''Bool evaluate elliptic function mod p equality: quad in y (a) vs cubic in x (b)'''
    y_sum = y**2 + y 
    x_sum = x**3 - x**2
    return (y_sum % p) == (x_sum % p)


########
# main is used to test selected functions
########

if __name__ == '__main__':
    '''elliptic module main: test functions. Run from the command line "python elliptic.py"'''
    print("Diagnostic printouts for elliptic module functions")
    print()
    print('EvaluateEF(1, 1, 1, 1, 1, 1, 1, 2, 3, 11):', EvaluateEF(1, 1, 1, 1, 1, 1, 1, 2, 3, 11))
    print()
