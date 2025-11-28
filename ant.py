# analytic number theory utilities
#   notes
#     - Some arithmetic function functions() return a list. Such a list does not index by integer argument
#       as this is too much hassle. So (mixing notation) mu[0] = mu(1) = 1 and mu[1] = mu(2) = -1


from math import sqrt, prod, log2, log10, log, floor, ceil, nan


######
# Test Python formalism: On import gives 'ant'; on direct run gives '__main__'.
######
# print("Within ant.py: __name__ is", __name__)


######
# Constants
######

# Euler's constant would be written ant.gamma in Notebook code


gamma = 0.57721566490153286060651209008240243104215933593992


def EulersConstant(n):
    '''Calculate an approximation of Euler's constant from n terms: Partial harmonic sum - log n'''
    ph_sum = 0
    for i in range(1, n+1): ph_sum += 1/i
    return ph_sum - log(n)


def Gamma(n): return EulersConstant(n)


######
# Basic arithmetic 
######


def even(a):
    '''even(n) returns True or False'''
    return False if a%2 else True


def odd(a): 
    '''odd(n) returns True or False'''
    return True if a%2 else False


def divides(a, b): 
    '''divides(a, b) "does a divide b?" returns True or False'''
    return not bool(b % a)


def prime(n):
    '''prime(n) returns True or False'''
    if not n%2: return False
    for i in range(3, int(sqrt(n)) + 1, 2): 
        if not n%i: return False
    return True


def is_prime(n): 
    '''is_prime(n) returns True or False'''
    return(prime(n))


def IsPrime(n):
    '''IsPrime(n) returns True or False'''
    return(prime(n))


def factor(n):
    '''factor(n) returns a sorted list: prime factorization (or [1])'''
    if n < 1: return('factor(n) error: argument not a positive integer')
    if n < 4: return [n]
    n_redux = n
    f, d, u = [], 2, int(sqrt(n_redux))
    while True:                                # ...continue til possible divisors exhausted
        while not n_redux % d:                 # ...all factors of value d
            f.append(d)
            n_redux = n_redux // d             # balance of original n remaining
        d = 3 if d == 2 else d + 2             # check 2 and odd numbers up to floor sqrt (n)
        if d > u: break                        # end of outer while
    if n_redux > 1: f.append(n_redux)          # include any residual > 1
    return f


def exponentfactor(n):
    '''exponentfactor(n > 1) returns a list of (prime, exponent) tuples'''
    if n < 1: print('exponentfactor(n) error: argument not a positive integer')
    if n == 1: return [(1, 1)]
    f = sorted(factor(n))           # f is a list
    p, ptr = [(f[0], 1)], 0         # p begins as first factor raised to the 1
    for i in f[1:]:                 # this loop will run if n is composite; 2nd factor onward
        if i == p[ptr][0]: p[ptr] = (i, p[ptr][1] + 1)       # immutable tuple construct
        else:              p.append((i, 1)); ptr += 1        # moving on to a new prime
    return(p)


def uniqueprimefactors(n):
    '''uniqueprimefactors(n > 1) returns a list of unique prime factors'''
    if n < 1: return('uniqueprimefactors(n) error:  argument not a positive integer')
    if n == 1: return [1]
    f, d, u = [], 2, int(sqrt(n))
    while True:
        while not n % d:
            if not d in f: f.append(d)
            n = n // d
        d = 3 if d == 2 else d + 2             # lazy: just check 2 and odds
        if d > u: break
    if n > 1 and not n in f: f.append(n)
    return f


def uniquefactors(n): 
    '''uniquefactors(n > 1) returns a list of unique prime factors'''
    return uniqueprimefactors(n)

        
def gcd(a, b):
    '''gcd(a, b) arguments in Z+ returns (a, b)'''
    if a < 1 or b < 1: return('gcd(a, b) error: at least one argument is less than 1')
    af, bf, g = factor(a), factor(b), 1
    for f in af:
        if f in bf: bf.remove(f); g *= f
    return g


def relativelyprime(a, b):
    '''relativelyprime(a, b) returns True or False'''
    return True if gcd(a, b) == 1 else False


def relprime(a, b): 
    '''relprime(a, b) returns True or False'''
    return relativelyprime(a, b)


def is_relativelyprime(a, b): 
    '''is_relativelyprime(a, b) returns True or False'''
    return relativelyprime(a, b)


def listproduct(l):
    '''listproduct(l) returns the product of all elements in l'''
    return prod(l)


def divisors(n):
    '''
    divisors(n) returns a sorted list of all divisors of n starting
    at element 0 having value 1.
    '''
    if n < 1: return
    d, f = [1], factor(n)                      # factor() is the list of primes with repetitions
    nf = len(f)                                # 60 will have nf = 4: [2, 2, 3, 5]
    for i in range(1, 2**nf):                  #   binomial theorem
        k = boolkey(i)                         # convert i to a binary list of booleans
        this_list = [1]                        # this_list will be a compiled list of factors from f[]
        for j in range(len(k)):                # scan through the boolean list k
            if k[j]: this_list.append(f[j])    #   ...appending factors from True values
        d.append(prod(this_list))              # append to the result this list-product
    return sorted(list(set(d)))


def divisors_less_n(n):
    '''
    divisors_less_n(n) returns a sorted list of all divisors of n except n itself
    starting at element 0 equals 1.
    '''
    if n < 1: return
    d, f = [1], factor(n)                      # factor() is the list of primes with repetitions
    nf = len(f)                                # 60 will have nf = 4: [2, 2, 3, 5]
    for i in range(1, 2**nf):                  #   binomial theorem
        k = boolkey(i)                         # convert i to a binary list of booleans
        this_list = [1]                        # this_list will be a compiled list of factors from f[]
        for j in range(len(k)):                # scan through the boolean list k
            if k[j]: this_list.append(f[j])    #   ...appending factors from True values
        d.append(prod(this_list))              # append to the result this list-product
    almost_done = sorted(list(set(d)))
    del(almost_done[-1])
    return(almost_done)


def boolkey(n):
    '''
    boolkey(n) returns list(bool(binary(n))) in ascending-power order
    '''
    key = [False]*(int(log2(n))+1)
    while n > 0:
        p = int(log2(n))
        key[p] = True
        n -= 2**p
    return key


#######
# Functions part 1
#######


def Totient(n):
    '''Totient(n) returns integer count of numbers < n that are relatively prime to n'''
    if n < 1: return 'Totient(n) error: argument not in Z+'
    totient = 1
    for i in range(2, n):
        if relativelyprime(n, i): totient += 1
    return totient


def totient(n): return Totient(n)


def GeneralizedTotient(k, n):
    '''GeneralizedTotient(k, n) returns integer sum of k'th powers of { m | (m < n, n) = 1 }'''
    if k < 0: return 'GeneralizedTotient(k, n) error: k not in 0, 1, ... '
    if n < 1: return 'GeneralizedTotient(k, n) error: n not in Z+'
    gtotient = 1
    for i in range(2, n):    # Only active for n > 2 notice
        if relativelyprime(n, i): gtotient += i**k
    return gtotient


def GeneralizedTotient0(n): return GeneralizedTotient(0, n)
def GeneralizedTotient1(n): return GeneralizedTotient(1, n)
def GeneralizedTotient2(n): return GeneralizedTotient(2, n)
def GeneralizedTotient3(n): return GeneralizedTotient(3, n)
def GeneralizedTotient4(n): return GeneralizedTotient(4, n)
def GeneralizedTotient5(n): return GeneralizedTotient(5, n)
def GeneralizedTotient6(n): return GeneralizedTotient(6, n)


def ExtendedTotient(x, n): 
    '''ExtendedTotient(float x, int n) returns integer count of rp-to-n for m <= floor(x)'''
    if x < 1: return 0
    if n < 1: return 'ExtendedTotient(float x, int n) error: n not in Z+'
    extended_totient = 0
    for i in range(1, floor(x) + 1):
        if relativelyprime(i, n): extended_totient += 1
    return extended_totient       


def Mobius(n):
    '''Mobius(n) returns integer mu(n)'''
    if n < 1: return 'Mobius(n) error: n not in Z+'
    if n == 1: return 1
    p = exponentfactor(n)
    for i in p: 
        if i[1] > 1: return 0
    return int((-1)**len(p))


def mobius(n): return Mobius(n)


def TotientMobiusProduct(n):
    '''TotientMobiusProduct(n) returns Totient(n) * Mobius(n)'''
    return Totient(n) * Mobius(n)


def Nu(n):
    '''Nu(n) returns integer count of unique prime factors of n'''
    if n < 1:  print('Nu(n) error: n must be in Z+')
    if n == 1: return 0
    return len(uniquefactors(n))
    
    
########
# Functions part 2: Zeta on the complex plain
########

def zeta(s):
    # Handles s types int, float or complex 
    if type(s) == int: s = float(s) 
    if type(s) == float:
        if s <= 0.: return 0.                    # stopgap
        if s < 1.: return zeta_s_on_0_1(s, 10000)
        if s == 1.: return float('inf')          # zeta(1) diverges to +-infinity
        if s > 1.: return zeta_s_gt_1(s, 10000)
    elif type(s) == complex:
        return 0.    # stopgap
    else:
        return 'type not recognized in zeta(s)'

def zeta_s_gt_1(s, n):
    '''
    Evaluate zeta > 1, where n is a precision parameter: Number of terms to sum,
    accurate to about one part in a million for n = 1000000. For further options: 
    See Wikipedia on Rational zeta series.
    n.b.: Basel problem s = 2 gives pi squared over six.
    '''
    zeta_sum = 0
    for i in range(1, n+1): zeta_sum += pow(1/i, s)
    return zeta_sum


# Incomplete, hardcoded
def zeta_s_on_0_1(s, n):
    '''
    Modify to return the limit of expr(x) as x goes to infinity: See Tommy p.55 for expr.
    This returns the value for s = 1/2: -1.46.
    '''
    return -1.46035450880958681288       # stopgap


########
# Dirichlet convolution
########
# 
# This code does Dirichlet convolutions
#   Dirichlet(fn1, fn2, n) uses two passed function names to calculate fn1 * fn2 at n
#   ListDirichlet(l1, l2, n) returns a Dirichlet product at n:
#     Assumption: Functions are passed as two lists (contrast with above)
#     Assumption: len(L1) = len(L2) = n
#     Convention: l1[0] is the function L1 evaluated at n = 1
#   ScanningListDirichlet(L1, L2, n) 
#     ...uses the same function-as-list approach for L1 and L2
#     ...returns a list of Dirichlet products for 1, 2, ..., n

 
def Dirichlet(fn1, fn2, n):
    '''Dirichlet(fn1, fn2, n) returns D'product of two functions evaluated at n'''
    if n < 1: print("Dirichlet(fn1, fn2, n) error: n not in Z+")
    if n == 1: return fn1(1)*fn2(1)
    total = 0
    for d in divisors(n): total += fn1(d) * fn2(n//d)
    return total


def dirichlet(fn1, fn2, n): 
    '''
    Return the Dirichlet product of two functions evaluated at n
    '''
    return Dirichlet(fn1, fn2, n)



def ListDirichlet(l1, l2, n):
    '''
    ListDirichlet(l1, l2, n) returns the Dirichlet product of two (list) functions 
    evaluated at n. The functions l1 and l2 are indexed [0, ..., n-1] so n corresponds
    to list index n-1.
    '''
    if n < 1:                      return False
    if len(l1) < n or len(l2) < n: return False
    dirichlet_sum = 0
    for d in divisors(n): dirichlet_sum += l1[d-1] * l2[n//d - 1]
    return dirichlet_sum


def ScanningListDirichlet(l1, l2):
    '''
    ScanningListDirichlet(l1, l2) returns a list of n Dirichlet products indexed
    0, 1, ... n-1. Here n = len(l1) presumed to equal len(l2). 
    '''
    Products = []
    if not len(l1) == len(l2): return False
    n = len(l1)
    if n < 1: return False
    for i in range(n): Products.append(ListDirichlet(l1, l2, i + 1))
    return(Products)


    
######
# Functions from number theory
######


def I(n):
    '''I(n) returns 1 if n is 1, 0 otherwise (the identity function for D'multiplication)'''
    if n < 1: return 'I(n) error: n not in Z+'
    if n == 1: return 1
    return 0


def U(n):
    '''U(n) returns 1 for all n (the unit function, Mobius inverse)'''
    if n < 1: return 'U(n) error: n not in Z+'
    return 1


def N(n):
    '''N(n) returns n'''
    if n < 1: return 'N(n) error: n not in Z+'
    return n


def Na(n, a):
    '''Na(n, a) returns n raised to the a power'''
    if n < 1: return 'Na(n) error: n not in Z+'
    return n**a

def Nk(n, k):
    '''Nk(n, k) returns n raised to the k power'''
    if n < 1: return 'Nk(n, k) error: n not in Z+'
    return n**k

# Convenience example power functions using exponent k
def Nk0(n): return Nk(n, 0)
def Nk1(n): return Nk(n, 1)
def Nk2(n): return Nk(n, 2)
def Nk3(n): return Nk(n, 3)
def Nk4(n): return Nk(n, 4)
def Nk5(n): return Nk(n, 5)
def Nk6(n): return Nk(n, 6)


def d(n):
    '''d(n) returns the number of divisors of n'''
    return len(divisors(n))


#######
# Arithmetic Mean (AM)
#######


def ArithmeticMean(fn, n):
    '''For upper limit n calculate the arithmetic mean of fn(i)'''
    return sum([fn(i) for i in range(1, n+1)])/n


def AMList(fn, n):
    '''Return a list of arithmetic means for fn over 1, 2, 3, ..., n'''
    return [ArithmeticMean(fn, i) for i in range(1, n+1)]


########
# main is used to test selected functions
########

if __name__ == '__main__':
    '''ant module main: test the functions. Run this from the command line "python ant.py"'''
    print("Diagnostic prints for all ant module functions")
    print("...not comprehensive but these act as a first cut at validation...")
    print()
    print('even(17):', even(17))
    print('odd(17):', odd(17))
    print('divides(7, 63):', divides(7, 63))
    print('prime(100), prime(101):', prime(100), prime(101))
    print('is_prime(100), is_prime(101):', is_prime(100), is_prime(101))
    print('IsPrime(100), IsPrime(101):', IsPrime(100), IsPrime(101))
    print('factor(180):', factor(180))
    print('exponentfactor(300):', exponentfactor(300))
    print('uniqueprimefactors(2*3*5*7*9*11):', uniqueprimefactors(2*3*5*7*9*11))
    print('uniqueprimefactors(600000):', uniqueprimefactors(600000))
    print('gcd(70, 1000):', gcd(70, 1000))
    print('relativelyprime(20, 65):', relativelyprime(20, 65))
    print('relprime(17, 35):', relprime(17, 35)) 
    print('is_relativelyprime(17, 35):', is_relativelyprime(17, 35))
    print('listproduct([1, 2, 3, 4]):', listproduct([1, 2, 3, 4]))
    print('divisors(210):', divisors(210)) 
    print('divisors_less_n(210):', divisors_less_n(210)) 
    print('boolkey(26): ', str(boolkey(26)), 'where correct is [F, T, F, T, T]')
    print('Totient(31):', Totient(31)) 
    print('GeneralizedTotient(2, 5):', GeneralizedTotient(2, 5))
    print('ExtendedTotient(17, 12):', ExtendedTotient(17, 12))  
    print('Mobius(1, 2, 3, 4, 5, 6, 7, 8):')
    for i in range(1, 9): print('       ', i, Mobius(i))
    print('TotientMobiusProduct(15...22):')
    for i in range(15, 23): print('      ', i, TotientMobiusProduct(i))
    print('Nu(210):', Nu(210)) 
    print('Dirichlet(Totient, Mobius, 18):', Dirichlet(Totient, Mobius, 18))
    print('ListDirichlet([1, 2, 3, 4, 6, 6, 7], [7, 6, 5, 4, 3, 2, 1], 7):',
          ListDirichlet([1, 2, 3, 4, 6, 6, 7], [7, 6, 5, 4, 3, 2, 1], 7))
    print('ScanningListDirichlet([1, 2, 3, 4, 6, 6, 7], [7, 6, 5, 4, 3, 2, 1]):',
          ScanningListDirichlet([1, 2, 3, 4, 6, 6, 7], [7, 6, 5, 4, 3, 2, 1]))
    print('I(10):', I(10))
    print('U(20):', U(20))
    print('N(73):', N(73), '(should be identity)')
    print('Na(4, 3):', Na(4, 3), '(should be a**b)')
    
    print('\nCompare generating function coefficient values with counting problem solutions')
    p_less_solns = []
    for p in [2, 3, 5, 7, 11, 13]:
        count_solutions = 0
        for x in range(p):
            for y in range(p):
                if EllipticTest(x, y, p): count_solutions += 1
        p_less_solns.append(p - count_solutions)
    print(p_less_solns)
    print([-2, -1, 1, -2, 1, 4])
    
                
                






