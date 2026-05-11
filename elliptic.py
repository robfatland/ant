# Elliptic curve utilities
#
# Motivation (Frenkel / Langlands program):
#   The number of solutions to an elliptic curve mod p can be predicted by
#   the coefficients of a modular form. Specifically for y^2 + y = x^3 - x^2,
#   the generating function is:
#     q * prod_{n=1}^{inf} (1-q^n)^2 * (1-q^{11n})^2
#   whose coefficients of q^p (for prime p) give a(p) = p - (number of solutions mod p).
#
#   Reference: Frenkel on the Langlands program
#   https://www.youtube.com/watch?v=4dyytPboqvE&t=910s
#
# Note on EvaluateEF:
#   This is a general-purpose evaluator for curves of the form
#     a0 + a1*y + a2*y^2 = b0 + b1*x + b2*x^2 + b3*x^3  (mod p)
#   It exists for potential exploration of other curves beyond the specific
#   y^2 + y = x^3 - x^2 case handled by EllipticTest. The original code had
#   a comment about an "odd error" when changing the cubic exponent to 2 —
#   this is expected: reducing the cubic to a quadratic changes the curve
#   entirely, collapsing it to a conic, which has exactly p solutions mod p
#   (hence count - p = 0, not the modular form coefficients).


def EvaluateEF(a0, a1, a2, b0, b1, b2, b3, x, y, p):
    '''
    EvaluateEF(a0, a1, a2, b0, b1, b2, b3, x, y, p)
    Returns True if a0 + a1*y + a2*y^2 == b0 + b1*x + b2*x^2 + b3*x^3 (mod p).
    General evaluator for Weierstrass-form elliptic curves.
    '''
    if p < 2:
        raise ValueError(f"EvaluateEF: p must be >= 2, got {p}")
    y_sum = a0 + a1 * y + a2 * (y**2)
    x_sum = b0 + b1 * x + b2 * (x**2) + b3 * (x**3)
    return (y_sum % p) == (x_sum % p)


def EllipticTest(x, y, p):
    '''
    EllipticTest(x, y, p)
    Returns True if y^2 + y == x^3 - x^2 (mod p).
    This is the specific curve tied to the Ramanujan tau / weight-2 level-11 modular form.
    '''
    if p < 2:
        raise ValueError(f"EllipticTest: p must be >= 2, got {p}")
    y_sum = y**2 + y
    x_sum = x**3 - x**2
    return (y_sum % p) == (x_sum % p)


def CountSolutions(p):
    '''
    CountSolutions(p) returns the number of (x, y) pairs in [0, p) x [0, p)
    satisfying y^2 + y == x^3 - x^2 (mod p).
    '''
    if p < 2:
        raise ValueError(f"CountSolutions: p must be >= 2, got {p}")
    count = 0
    for x in range(p):
        for y in range(p):
            if EllipticTest(x, y, p):
                count += 1
    return count


def CoefficientA(p):
    '''
    CoefficientA(p) returns a(p) = p - CountSolutions(p).
    This should match the p-th coefficient of the associated modular form.
    '''
    return p - CountSolutions(p)


########
# main: test elliptic functions
########

if __name__ == '__main__':
    print("Diagnostic prints for elliptic module functions")
    print("=" * 60)
    print()
    print("EvaluateEF(0, 1, 1, 0, 0, -1, 1, 2, 3, 11):",
          EvaluateEF(0, 1, 1, 0, 0, -1, 1, 2, 3, 11))
    print()
    print("Counting solutions and computing a(p) = p - #solutions")
    print("for y^2 + y = x^3 - x^2 mod p:")
    print()
    print(f"  {'p':>4}  {'#solutions':>10}  {'a(p)':>6}")
    print(f"  {'---':>4}  {'----------':>10}  {'----':>6}")
    expected = {2: -2, 3: -1, 5: 1, 7: -2, 11: 1, 13: 4}
    for p in [2, 3, 5, 7, 11, 13]:
        count = CountSolutions(p)
        a_p = p - count
        status = "OK" if a_p == expected[p] else "MISMATCH"
        print(f"  {p:>4}  {count:>10}  {a_p:>6}  {status}")
    print()
    print("Expected a(p) values: [-2, -1, 1, -2, 1, 4]")
