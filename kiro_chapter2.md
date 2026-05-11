# Kiro's Solutions: Problems 2.24 and 2.34

## Problem 2.24: Formal power series have no zero divisors

**Statement:** If the product of two formal power series $A(x)$ and $B(x)$ is the zero series,
prove that at least one factor is the zero series.

**Proof:**

Let

$$A(x) = \sum_{n=0}^{\infty} a_n x^n, \quad B(x) = \sum_{n=0}^{\infty} b_n x^n$$

and suppose $A(x) \cdot B(x) = 0$, meaning every coefficient of the product is zero:

$$\sum_{k=0}^{n} a_k \, b_{n-k} = 0 \quad \text{for all } n \geq 0$$

Assume for contradiction that neither $A$ nor $B$ is the zero series. Then there exist
smallest indices $r$ and $s$ such that $a_r \neq 0$ and $b_s \neq 0$ (with all earlier
coefficients being zero in their respective series).

Now examine the coefficient of $x^{r+s}$ in the product:

$$\sum_{k=0}^{r+s} a_k \, b_{r+s-k} = 0$$

- For $k < r$: $a_k = 0$ by minimality of $r$.
- For $k > r$: we need $b_{r+s-k}$ where $r+s-k < s$, so $b_{r+s-k} = 0$ by minimality of $s$.

The only surviving term is $k = r$:

$$a_r \cdot b_s = 0$$

But $a_r \neq 0$ and $b_s \neq 0$, and we are working over the integers (or more generally
any integral domain), which have no zero divisors. Contradiction.

Therefore at least one of $A$ or $B$ must be the zero series. $\square$

**Remarks:**

The technique here is sometimes called the "lowest-order term" or "leading term" argument.
You find the first nonzero coefficient in each series and show their product cannot vanish.
This is the same reason polynomial rings over integral domains are themselves integral domains —
the argument generalizes from finite to infinite degree without difficulty.

Note that the parenthetical in the problem text — "the ring of formal power series has no
zero divisors" — is exactly what we proved. A *zero divisor* in a ring is a nonzero element
whose product with some other nonzero element is zero. We showed no such pair exists among
formal power series with integer coefficients.

---

## Problem 2.34: $\lambda(n)|\mu(n)| = \mu(n)$

**Statement:** Prove that $\lambda(n) \cdot |\mu(n)| = \mu(n)$ for all $n \geq 1$.

(Here $\lambda$ is Liouville's function and $\mu$ is the Möbius function.)

**Proof:**

We proceed by cases on the prime factorization of $n$.

**Case 1: $n = 1$.**

$\lambda(1) \cdot |\mu(1)| = 1 \cdot 1 = 1 = \mu(1)$. $\checkmark$

**Case 2: $n$ is squarefree, $n > 1$.**

Write $n = p_1 p_2 \cdots p_r$ (distinct primes). Then:

- $\Omega(n) = r$, so $\lambda(n) = (-1)^r$
- $n$ is squarefree, so $|\mu(n)| = 1$
- $\mu(n) = (-1)^r$

Therefore $\lambda(n) \cdot |\mu(n)| = (-1)^r \cdot 1 = (-1)^r = \mu(n)$. $\checkmark$

**Case 3: $n$ is not squarefree.**

Some $p^2 \mid n$, so $\mu(n) = 0$ and $|\mu(n)| = 0$.

The left side is $\lambda(n) \cdot 0 = 0 = \mu(n)$. $\checkmark$

These three cases are exhaustive, so the identity holds for all $n \geq 1$. $\square$

**Remarks:**

The identity says: $\mu(n)$ is just Liouville's function restricted to squarefree integers
(and zero elsewhere). This makes sense intuitively — for squarefree $n$, the number of
prime factors counted with multiplicity ($\Omega$) equals the number of distinct prime
factors ($\omega$), so $\lambda$ and $\mu$ agree. The $|\mu(n)|$ factor acts as the
"squarefree indicator," killing everything else.

An equivalent way to state this: $\mu = \lambda \cdot |\mu|$, which is a pointwise product
identity (not a Dirichlet convolution).

---

## Kun?

1. **On 2.24:** I assumed the coefficient ring is the integers (or at least an integral domain).
   Apostol's formal power series in Chapter 2 — is he working over $\mathbb{Z}$, or does he
   leave the coefficient ring general? If it's over an arbitrary commutative ring, the result
   fails (e.g., formal power series over $\mathbb{Z}/6\mathbb{Z}$ *do* have zero divisors).
   Does your edition specify?

2. **On 2.34:** I'm not 100% certain this is the exact statement in your copy of Apostol.
   The problem slot in your notebook is just "$\dots$" with no statement. Can you confirm
   the problem text? If it's something else (perhaps involving the generalized Möbius
   function $\mu_k$ since that's the surrounding context), I'll redo it.

3. **On Chip's solutions:** You mentioned struggling with his approach to these two. Was it
   the *technique* that was unclear, or the *notation/presentation*? For 2.24 in particular,
   some solution sets use a more algebraic/ring-theoretic framing that can obscure what is
   really a simple "first nonzero coefficient" argument.
