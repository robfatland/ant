# analytic number theory utilities
#
# Conventions:
#   - Function names use CamelCase (no underscores)
#   - All functions raise ValueError for invalid arguments
#   - Some functions return lists indexed from 0; so result[0] = f(1)
#   - This module requires mpmath for zeta function evaluation


from math import sqrt, prod, log2, log, floor, ceil
from math import gcd as _gcd
from functools import reduce
import mpmath


######
# Constants
######


gamma = 0.57721566490153286060651209008240243104215933593992


def EulersConstant(n):
    '''Calculate an approximation of Euler's constant from n terms: Partial harmonic sum - log n'''
    if n < 1:
        raise ValueError(f"EulersConstant(n): n must be a positive integer, got {n}")
    ph_sum = 0
    for i in range(1, n + 1):
        ph_sum += 1 / i
    return ph_sum - log(n)


def Gamma(n):
    '''Alias for EulersConstant(n)'''
    return EulersConstant(n)


######
# Basic arithmetic
######


def Even(a):
    '''Even(n) returns True or False'''
    return a % 2 == 0


def Odd(a):
    '''Odd(n) returns True or False'''
    return a % 2 != 0


def Divides(a, b):
    '''Divides(a, b): does a divide b? Returns True or False'''
    if a == 0:
        raise ValueError("Divides(a, b): a must be nonzero")
    return b % a == 0


def Prime(n):
    '''Prime(n) returns True if n is prime, False otherwise'''
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False
    for i in range(3, int(sqrt(n)) + 1, 2):
        if n % i == 0:
            return False
    return True


def Factor(n):
    '''Factor(n) returns a sorted list: prime factorization of n (with repetition)'''
    if n < 1:
        raise ValueError(f"Factor(n): argument must be a positive integer, got {n}")
    if n == 1:
        return [1]
    if n < 4:
        return [n]
    n_redux = n
    f, d, u = [], 2, int(sqrt(n_redux))
    while True:
        while n_redux % d == 0:
            f.append(d)
            n_redux = n_redux // d
        d = 3 if d == 2 else d + 2
        if d > u:
            break
    if n_redux > 1:
        f.append(n_redux)
    return f


def ExponentFactor(n):
    '''ExponentFactor(n) returns a list of (prime, exponent) tuples'''
    if n < 1:
        raise ValueError(f"ExponentFactor(n): argument must be a positive integer, got {n}")
    if n == 1:
        return [(1, 1)]
    f = sorted(Factor(n))
    p, ptr = [(f[0], 1)], 0
    for i in f[1:]:
        if i == p[ptr][0]:
            p[ptr] = (i, p[ptr][1] + 1)
        else:
            p.append((i, 1))
            ptr += 1
    return p


def UniquePrimeFactors(n):
    '''UniquePrimeFactors(n) returns a list of unique prime factors'''
    if n < 1:
        raise ValueError(f"UniquePrimeFactors(n): argument must be a positive integer, got {n}")
    if n == 1:
        return [1]
    f, d, u = [], 2, int(sqrt(n))
    while True:
        while n % d == 0:
            if d not in f:
                f.append(d)
            n = n // d
        d = 3 if d == 2 else d + 2
        if d > u:
            break
    if n > 1 and n not in f:
        f.append(n)
    return f


def Gcd(a, b):
    '''Gcd(a, b) returns the greatest common divisor using the Euclidean algorithm'''
    if a < 1 or b < 1:
        raise ValueError(f"Gcd(a, b): both arguments must be positive integers, got ({a}, {b})")
    return _gcd(a, b)


def RelativelyPrime(a, b):
    '''RelativelyPrime(a, b) returns True if gcd(a, b) == 1'''
    return Gcd(a, b) == 1


def ListProduct(l):
    '''ListProduct(l) returns the product of all elements in l'''
    return prod(l)


def Divisors(n):
    '''
    Divisors(n) returns a sorted list of all divisors of n.
    Built from the exponent factorization: generates all combinations of prime powers.
    '''
    if n < 1:
        raise ValueError(f"Divisors(n): argument must be a positive integer, got {n}")
    if n == 1:
        return [1]
    ef = ExponentFactor(n)
    # Start with [1] and iteratively multiply in each prime's powers
    divs = [1]
    for (p, a) in ef:
        if p == 1:
            continue
        new_divs = []
        for d in divs:
            pk = 1
            for _ in range(a + 1):
                new_divs.append(d * pk)
                pk *= p
        divs = new_divs
    return sorted(divs)


DivisorsLessN = lambda n: Divisors(n)[:-1]
DivisorsLessN.__doc__ = '''DivisorsLessN(n) returns all divisors of n except n itself'''


#######
# Arithmetic functions
#######


def Totient(n):
    '''Totient(n) returns integer count of numbers <= n that are relatively prime to n'''
    if n < 1:
        raise ValueError(f"Totient(n): argument must be a positive integer, got {n}")
    totient = 1
    for i in range(2, n):
        if RelativelyPrime(n, i):
            totient += 1
    return totient


def GeneralizedTotient(k, n):
    '''GeneralizedTotient(k, n) returns sum of k-th powers of { m : 1 <= m < n, (m, n) = 1 }'''
    if k < 0:
        raise ValueError(f"GeneralizedTotient(k, n): k must be non-negative, got {k}")
    if n < 1:
        raise ValueError(f"GeneralizedTotient(k, n): n must be a positive integer, got {n}")
    gtotient = 1
    for i in range(2, n):
        if RelativelyPrime(n, i):
            gtotient += i**k
    return gtotient


def ExtendedTotient(x, n):
    '''ExtendedTotient(x, n) returns count of integers m <= floor(x) with (m, n) = 1'''
    if x < 1:
        return 0
    if n < 1:
        raise ValueError(f"ExtendedTotient(x, n): n must be a positive integer, got {n}")
    extended_totient = 0
    for i in range(1, floor(x) + 1):
        if RelativelyPrime(i, n):
            extended_totient += 1
    return extended_totient


def Mobius(n):
    '''Mobius(n) returns mu(n)'''
    if n < 1:
        raise ValueError(f"Mobius(n): argument must be a positive integer, got {n}")
    if n == 1:
        return 1
    p = ExponentFactor(n)
    for i in p:
        if i[1] > 1:
            return 0
    return int((-1)**len(p))


def TotientMobiusProduct(n):
    '''TotientMobiusProduct(n) returns Totient(n) * Mobius(n)'''
    return Totient(n) * Mobius(n)


def Nu(n):
    '''Nu(n) returns integer count of unique prime factors of n'''
    if n < 1:
        raise ValueError(f"Nu(n): argument must be a positive integer, got {n}")
    if n == 1:
        return 0
    return len(UniquePrimeFactors(n))


########
# Zeta function (via mpmath for precision and analytic continuation)
########


def Zeta(s, dps=15):
    '''
    Zeta(s) returns the Riemann zeta function evaluated at s.
    Uses mpmath for arbitrary precision and full analytic continuation.
    Optional dps parameter sets decimal places of precision (default 15).
    '''
    old_dps = mpmath.mp.dps
    mpmath.mp.dps = dps
    try:
        result = complex(mpmath.zeta(s))
        # Return float if the result is real
        if result.imag == 0:
            return result.real
        return result
    finally:
        mpmath.mp.dps = old_dps


########
# Dirichlet convolution
########


def Dirichlet(fn1, fn2, n):
    '''Dirichlet(fn1, fn2, n) returns the Dirichlet product of two functions evaluated at n'''
    if n < 1:
        raise ValueError(f"Dirichlet(fn1, fn2, n): n must be a positive integer, got {n}")
    if n == 1:
        return fn1(1) * fn2(1)
    total = 0
    for d in Divisors(n):
        total += fn1(d) * fn2(n // d)
    return total


def ListDirichlet(l1, l2, n):
    '''
    ListDirichlet(l1, l2, n) returns the Dirichlet product of two list-functions
    evaluated at n. Lists are indexed [0, ..., n-1] so index n-1 corresponds to argument n.
    '''
    if n < 1:
        raise ValueError(f"ListDirichlet(l1, l2, n): n must be a positive integer, got {n}")
    if len(l1) < n or len(l2) < n:
        raise ValueError(f"ListDirichlet: lists must have length >= {n}")
    dirichlet_sum = 0
    for d in Divisors(n):
        dirichlet_sum += l1[d - 1] * l2[n // d - 1]
    return dirichlet_sum


def ScanningListDirichlet(l1, l2):
    '''
    ScanningListDirichlet(l1, l2) returns a list of n Dirichlet products
    indexed 0, 1, ..., n-1 where n = len(l1) = len(l2).
    '''
    if len(l1) != len(l2):
        raise ValueError("ScanningListDirichlet: lists must have equal length")
    n = len(l1)
    if n < 1:
        raise ValueError("ScanningListDirichlet: lists must be non-empty")
    return [ListDirichlet(l1, l2, i + 1) for i in range(n)]


######
# Standard arithmetic functions (Apostol notation)
######


def I(n):
    '''I(n): the identity for Dirichlet multiplication. Returns 1 if n==1, else 0.'''
    if n < 1:
        raise ValueError(f"I(n): argument must be a positive integer, got {n}")
    return 1 if n == 1 else 0


def U(n):
    '''U(n): the unit function, returns 1 for all n in Z+.'''
    if n < 1:
        raise ValueError(f"U(n): argument must be a positive integer, got {n}")
    return 1


def N(n):
    '''N(n): the identity function, returns n.'''
    if n < 1:
        raise ValueError(f"N(n): argument must be a positive integer, got {n}")
    return n


def Nk(n, k):
    '''Nk(n, k) returns n raised to the k power'''
    if n < 1:
        raise ValueError(f"Nk(n, k): n must be a positive integer, got {n}")
    return n**k


# For convenience in Dirichlet() calls which expect single-argument functions,
# use closures rather than a family of named functions:
#   NkFunc(k) returns a function n -> n**k
#   GeneralizedTotientFunc(k) returns a function n -> GeneralizedTotient(k, n)

NkFunc = lambda k: (lambda n: Nk(n, k))
NkFunc.__doc__ = '''NkFunc(k) returns a closure: a function of n that computes n**k'''

GeneralizedTotientFunc = lambda k: (lambda n: GeneralizedTotient(k, n))
GeneralizedTotientFunc.__doc__ = '''GeneralizedTotientFunc(k) returns a closure: a function of n that computes GeneralizedTotient(k, n)'''


def D(n):
    '''D(n) returns the number of divisors of n (i.e. sigma_0)'''
    return len(Divisors(n))


#######
# Arithmetic Mean
#######


def ArithmeticMean(fn, n):
    '''For upper limit n calculate the arithmetic mean of fn(i) for i = 1, ..., n'''
    if n < 1:
        raise ValueError(f"ArithmeticMean(fn, n): n must be a positive integer, got {n}")
    return sum([fn(i) for i in range(1, n + 1)]) / n


def AMList(fn, n):
    '''Return a list of arithmetic means for fn over 1, 2, 3, ..., n'''
    return [ArithmeticMean(fn, i) for i in range(1, n + 1)]


########
# main: test selected functions
########

if __name__ == '__main__':
    from elliptic import EllipticTest

    print("Diagnostic prints for ant module functions")
    print("=" * 60)
    print()
    print('Even(17):', Even(17))
    print('Odd(17):', Odd(17))
    print('Divides(7, 63):', Divides(7, 63))
    print('Prime(2):', Prime(2), '(should be True)')
    print('Prime(100), Prime(101):', Prime(100), Prime(101))
    print('Factor(180):', Factor(180))
    print('ExponentFactor(300):', ExponentFactor(300))
    print('UniquePrimeFactors(2*3*5*7*9*11):', UniquePrimeFactors(2 * 3 * 5 * 7 * 9 * 11))
    print('UniquePrimeFactors(600000):', UniquePrimeFactors(600000))
    print('Gcd(70, 1000):', Gcd(70, 1000))
    print('RelativelyPrime(20, 65):', RelativelyPrime(20, 65))
    print('RelativelyPrime(17, 35):', RelativelyPrime(17, 35))
    print('ListProduct([1, 2, 3, 4]):', ListProduct([1, 2, 3, 4]))
    print('Divisors(210):', Divisors(210))
    print('DivisorsLessN(210):', DivisorsLessN(210))
    print('Totient(31):', Totient(31))
    print('GeneralizedTotient(2, 5):', GeneralizedTotient(2, 5))
    print('ExtendedTotient(17, 12):', ExtendedTotient(17, 12))
    print('Mobius(1..8):')
    for i in range(1, 9):
        print(f'       {i}: {Mobius(i)}')
    print('TotientMobiusProduct(15..22):')
    for i in range(15, 23):
        print(f'       {i}: {TotientMobiusProduct(i)}')
    print('Nu(210):', Nu(210))
    print('Dirichlet(Totient, Mobius, 18):', Dirichlet(Totient, Mobius, 18))
    print('ListDirichlet([1,2,3,4,6,6,7], [7,6,5,4,3,2,1], 7):',
          ListDirichlet([1, 2, 3, 4, 6, 6, 7], [7, 6, 5, 4, 3, 2, 1], 7))
    print('ScanningListDirichlet([1,2,3,4,6,6,7], [7,6,5,4,3,2,1]):',
          ScanningListDirichlet([1, 2, 3, 4, 6, 6, 7], [7, 6, 5, 4, 3, 2, 1]))
    print('I(10):', I(10))
    print('U(20):', U(20))
    print('N(73):', N(73), '(should be 73)')
    print('Nk(4, 3):', Nk(4, 3), '(should be 64)')
    print('NkFunc(3)(4):', NkFunc(3)(4), '(should be 64)')
    print()
    print('Zeta(2):', Zeta(2), f'(should be ~{3.14159**2/6:.10f})')
    print('Zeta(0):', Zeta(0), '(should be -0.5)')
    print('Zeta(-1):', Zeta(-1), '(should be -1/12)')
    print()

    print('Elliptic curve test: comparing generating function coefficients')
    print('  with counting problem solutions for y^2 + y = x^3 - x^2 mod p')
    p_less_solns = []
    for p in [2, 3, 5, 7, 11, 13]:
        count_solutions = 0
        for x in range(p):
            for y in range(p):
                if EllipticTest(x, y, p):
                    count_solutions += 1
        p_less_solns.append(p - count_solutions)
    print('  Computed: ', p_less_solns)
    print('  Expected: ', [-2, -1, 1, -2, 1, 4])
