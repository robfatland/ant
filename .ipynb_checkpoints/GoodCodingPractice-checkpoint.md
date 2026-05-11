# Good Coding Practice

A reference for the `ant` project. Written for a programmer who knows what they're doing
but never formalized certain habits.

---

## Error Handling: Exceptions over Print Statements

### The Problem with Ad Hoc Error Handling

You've probably written code like this a hundred times:

```python
def factor(n):
    if n < 1:
        print("error: n must be positive")
        return None
```

This fails silently. The caller gets `None` back, tries to iterate over it, and you get
a confusing `TypeError` three function calls later with no indication of where the real
problem was. Or worse: you return an error *string* and the caller treats it as valid data.

### The Fix: Raise Exceptions

```python
def Factor(n):
    if n < 1:
        raise ValueError(f"Factor(n): argument must be a positive integer, got {n}")
    ...
```

What this buys you:

1. **Immediate failure at the source.** The traceback points to the exact line where the
   contract was violated, not some downstream consequence.

2. **The caller can choose.** They can catch it with `try/except` if they want to handle
   it gracefully, or let it propagate if it represents a genuine bug.

3. **No ambiguous return types.** Your function always returns a list (or int, or whatever
   it promises). Never a string masquerading as an error message.

### Which Exception to Use

- `ValueError` — the argument is the wrong value (negative n, zero divisor, etc.)
- `TypeError` — the argument is the wrong type (string where int expected)
- `ZeroDivisionError` — self-explanatory
- `NotImplementedError` — stub for future work

For a math library like `ant.py`, `ValueError` covers 95% of cases. You're validating
that arguments are in the domain of the function.

### Pattern

```python
def MyFunction(n):
    '''Docstring explaining what it does'''
    if n < 1:
        raise ValueError(f"MyFunction(n): n must be a positive integer, got {n}")
    # ... actual logic ...
    return result
```

Include the function name and the offending value in the message. When you see the
traceback at 2am, you'll thank yourself.

### Catching Exceptions (Caller Side)

```python
try:
    result = Factor(-5)
except ValueError as e:
    print(f"Handled gracefully: {e}")
    result = None  # or some fallback
```

In notebook code you'll rarely write try/except — you *want* the exception to surface
so you can see what went wrong. The try/except pattern is for production code or loops
where one bad input shouldn't kill the whole run.

---

## API Design: Naming Conventions

### The Rule for This Project

All public function names use **CamelCase** (also called PascalCase):

```python
Factor(n)           # not factor(n)
ExponentFactor(n)   # not exponent_factor(n)
RelativelyPrime(a, b)  # not relativelyprime(a, b) or is_relatively_prime(a, b)
```

No underscores in function names. No lowercase-only names. One convention, everywhere.

### Why This Matters

An API is a contract. When you `import ant` in a notebook, you need to be able to *guess*
the function name without looking it up. If some functions are `totient`, some are `Totient`,
and some are `is_prime` — you can't guess. You have to check every time.

With a single convention, the rule is simple: "It's always CamelCase." You know it's
`RelativelyPrime`, not `relprime` or `is_relativelyprime`.

### Aliases Are Debt

The old code had:

```python
def prime(n): ...
def is_prime(n): return prime(n)
def IsPrime(n): return prime(n)
```

Three names for one function. This feels helpful ("the user can call whichever they
remember") but it triples the surface area of the API and makes autocomplete useless.
Pick one name. Use it. Delete the rest.

### When You Have a Family of Related Functions

Don't create `Nk0`, `Nk1`, `Nk2`, ..., `Nk6`. Use a factory (see Lambda section below):

```python
NkFunc = lambda k: (lambda n: Nk(n, k))

# Usage in notebook:
cube = NkFunc(3)
cube(4)  # returns 64
Dirichlet(Totient, NkFunc(2), 18)  # pass n^2 as a single-argument function
```

---

## Lambda Functions and Closures

### What a Lambda Is

A lambda is an anonymous function — a function without a name, defined inline:

```python
square = lambda n: n**2
```

This is exactly equivalent to:

```python
def square(n):
    return n**2
```

Use lambdas when the function is short (one expression) and you don't need a docstring
at the point of definition.

### Closures: Functions That Remember

A closure is a function that captures variables from its enclosing scope:

```python
NkFunc = lambda k: (lambda n: n**k)
```

`NkFunc(3)` returns a *new function* that remembers k=3. That returned function takes
one argument `n` and computes `n**3`. This eliminates the need for:

```python
def Nk0(n): return Nk(n, 0)
def Nk1(n): return Nk(n, 1)
def Nk2(n): return Nk(n, 2)
# ... ad nauseam
```

### The DivisorsLessN Pattern

```python
DivisorsLessN = lambda n: Divisors(n)[:-1]
```

This is a one-liner that delegates to `Divisors` and slices off the last element.
It avoids duplicating the entire divisor-generation logic. If you fix a bug in `Divisors`,
`DivisorsLessN` gets the fix for free.

You can attach a docstring after the fact:

```python
DivisorsLessN.__doc__ = '''DivisorsLessN(n) returns all divisors of n except n itself'''
```

### When NOT to Use Lambda

- When the function is more than one expression
- When you need proper error handling inside it
- When the logic is complex enough to deserve a name and docstring at definition time

Rule of thumb: if you can't read the lambda in one glance, make it a `def`.

---

## mpmath: When the Standard Library Isn't Enough

### What mpmath Is

`mpmath` is a pure-Python library for arbitrary-precision floating-point arithmetic.
It provides:

- Arbitrary decimal precision (set `mpmath.mp.dps = 50` for 50 digits)
- Full implementations of special functions (zeta, gamma, hypergeometric, elliptic, ...)
- Analytic continuation (zeta at negative integers, complex arguments, etc.)
- Numerical integration, differentiation, root-finding

### When It Supersedes `math`

The standard `math` library gives you hardware-precision (64-bit float, ~15 digits)
implementations of elementary functions: `sqrt`, `log`, `sin`, `exp`, etc.

Use `mpmath` when you need:

1. **More precision.** `math.pi` gives 15 digits. `mpmath.pi` gives as many as you ask for.

2. **Special functions.** `math` has no zeta function, no Bernoulli numbers, no
   polylogarithm, no hypergeometric functions. `mpmath` has all of these.

3. **Complex arguments.** `math.sqrt(-1)` raises an error. `mpmath.sqrt(-1)` returns `1j`.

4. **Analytic continuation.** The naive zeta series only converges for Re(s) > 1.
   `mpmath.zeta(s)` works for all complex s ≠ 1, using the functional equation and
   other techniques internally.

### Usage in This Project

```python
import mpmath

# Set precision
mpmath.mp.dps = 25  # 25 decimal places

# Evaluate zeta
mpmath.zeta(2)       # pi^2/6 to 25 digits
mpmath.zeta(0)       # -1/2 (analytic continuation)
mpmath.zeta(-1)      # -1/12 (Ramanujan's famous result)
mpmath.zeta(0.5 + 14.134725j)  # near a zero of zeta
```

### Installation

```
pip install mpmath
```

It's a dependency of `sympy`, so if you have sympy installed you already have mpmath.

### The Old Approach vs. The New

Old (`ant.py` before cleanup):
```python
def zeta_s_gt_1(s, n):
    zeta_sum = 0
    for i in range(1, n+1): zeta_sum += pow(1/i, s)
    return zeta_sum
```

This sums 10000 terms of 1/n^s. For s=2 the error is roughly 1/10000 ≈ 0.0001.
For s close to 1 it's much worse. And it can't handle s < 1 at all.

New:
```python
def Zeta(s, dps=15):
    mpmath.mp.dps = dps
    return float(mpmath.zeta(s))
```

Full analytic continuation, arbitrary precision, one line.

---

## Overall Architecture: Is a Two-File Library the Right Track?

Yes. The split makes sense:

- `ant.py` — general number-theoretic utilities tied to Apostol's text
- `elliptic.py` — specialized elliptic curve / modular form exploration

As the project grows (Chapter 3 onward: averages, summation formulas, Chebyshev),
you might eventually want:

```
ant.py          — core arithmetic functions, Dirichlet machinery
analytic.py     — asymptotic estimates, summation formulas, O-notation tools
elliptic.py     — elliptic curves and modular forms
```

But don't split prematurely. The current two-file structure is fine until `ant.py`
grows past ~500 lines or you find yourself scrolling past unrelated code to find
what you need. At that point, factor out a new module.

The key discipline: every function lives in exactly one place, has exactly one name,
raises exceptions on bad input, and has a docstring. Everything else is style.
