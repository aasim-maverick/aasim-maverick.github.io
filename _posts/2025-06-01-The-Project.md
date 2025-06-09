---
title: FLINT-Backed Multivariate Polynomials
date: 2025-01-19 00:00:00 +0000
categories: [blog]
tags: [sympy, gsoc25]
math: true
mermaid: true
comments: true
---


## **The Project: FLINT-Backed Multivariate Polynomials in SymPy**

Let’s kick off at the true “core” of `SymPy`. An intuitive guess might point one to the `sympy/core` directory, and they wouldn’t be wrong, but the idea of “core” really spans three subsystems, as I learned from my mentor Oscar Benjamin:

- **Numerical evaluation subsystem** (`sympy.core.evalf`)
- **Symbolic manipulation system** (`sympy.core.basic`)
- **Computational algebra subsystem** (`sympy.polys`)

If you’re new to `SymPy`, Oscar’s [blog](https://oscarbenjamin.github.io/blog/czi/post1.html) is a fantastic primer on how those three pieces fit together. I strongly recommend all contributors, new or old, to take some time out to read that. 

---

This project zeroes in on that third pillar, the **`polys`** module, and aims to supercharge its performance by tapping into `FLINT`’s C-level speed. Since `polys` underpins so much of SymPy’s heavy lifting (think factorization, GCDs, polynomial comparison, and more), any slowdown here ripples across the entire library.

In `SymPy` today:

- **Sparse multivariate polynomials** live in `sympy/polys/rings.py` via the `PolyRing` and `PolyElement` classes.
- **Dense polynomials** are handled by `DMP` in `sympy/polys/polyclasses.py`.

Almost[^1] all of that logic is pure Python, which is great for flexibility but can become a bottleneck when you’re dealing with high-dimensional polynomial arithmetic. By layering in FLINT’s blazing-fast C routines, wrapped through `python-flint`, we’ll keep the same clean API but give `SymPy` the speed boost it needs.

In short, this isn’t just about one module: it’s about empowering every feature that relies on fast, robust polynomial math across `SymPy`.

---

## **What Is the Project?**

This project aims to integrate `FLINT`-backed multivariate polynomials into `SymPy`, using the `python-flint` library. While `SymPy` already uses `FLINT` for univariate dense polynomials through DUP_Flint, multivariate support has been missing.

This project brings `FLINT`'s blazing-fast C routines into `SymPy` for both sparse and dense multivariate polynomial computations over integer (ℤ) and rational (ℚ) domains.

---

## **Why it matters**

Multivariate polynomials are used everywhere in `SymPy`:

- **Equation solving** (`solve`, `linsolve`)

- **Expression simplification** (`simplify`, `cancel`)

- **Symbolic integration** (e.g. `heurisch`)

- **Matrix computations**

- **Groebner basis calculations**

- Lots and lots of other important things

But `SymPy` currently relies on pure Python implementations for multivariate cases, leading to performance bottlenecks.

---

## **Dense vs Sparse: A Quick Primer**

A dense representation of a polynomial sotres every coefficent in a list, even zeroes, in an ascending or descending order of degree of the monomials. 

For example, let's say we have the following univariate polynomial:

$5x^5 +  2x^2 + 3^x + 1$

In `SymPy`, this can be represented internally in the dense representation:

```python
from sympy import Poly
from sympy.abc import x
p = Poly(5*x**5 + 2*x**2 + 3*x + 1)
print(p.rep.to_list())  # Outputs [5, 0, 0, 2, 3, 1]
```

> `p.rep` gives you access to the internal `DMP` representation of the `Poly` object `p`. In this univariate case, it specifically refers to a `DUP_Flint` instance, which is a **subclass** of the **base** `DMP`. If you want to convert it to a regular `Python` **list**, it's as simple as calling the `to_list()` method.
{: .prompt-info }

For the multivariate case's dense representation, we just apply the univariate dense representation recursively giving a nested list structure representing the polynomial.

For example, consider the multivariate polynomial:

$x^2 + x^2y^3 + xy^2 + 1$

Rewriting it as a polynomial in `x` with coefficients in `y`: 

$p(x, y) = (y^3 + 1)x^2 + (y^2)x + 1$

```python
from sympy import *
x, y = symbols('x, y')
p = x**2 + x**2*y**3 + x*y**2 + 1
print(p.collect(x)) # Outputs p as function in x with coefficients in y
print(p.as_poly(x, y).rep.to_list() ) # Outputs [[1, 0, 0, 1], [1, 0, 0], [1]]
```

Here,
- [[1, 0, 0, 1], [1, 0, 0], [1]] represents the structure where
the outermost list corresponds to powers of x, and

- Each nested list represents the coefficients, which are themselves polynomials
in y.

One can eventually see just how wasteful can this representation get when there are a lot of 0 terms and only a few non-zero ones in the polynomial. Eg. $x^{1000} + 1$ will contain a list of length 1001 and only 2 of those entries will be non zero(i.e. the useful qunantities).

```python
from sympy import *
x = symbols('x')
p = Poly(x**1000 + 1)
print(len(p.rep.to_list())) # 1001

dense_repr = p.rep.to_list()

for i in range(len(dense_repr)):
    if dense_repr[i]:
        print(dense_repr[i]) # Prints the coefficients 1 and 1 for x**1000 and x**0
```

This gets extremely inefficient, especially when the case is multivariate and there is a lot of nesting, the polynomial operations designed/programmed for the dense representation iterate through these nested lists and that is extremely slow. 

To overcome the above mentioned problem, the sparse polynomial representation is used, and that has entirely separate ways of programming the same operations. 

Say we have, $x^100 + x^2 + 1$, as discussed, representing this in the dense format is wasteful since we would need a list with 101 entries with 98 of them just 0s. 
Instead of storing all coefficients (including zeros), we only store nonzero terms
along with their corresponding degrees. This is done using a dictionary(`dict`) structure in
`SymPy`, where:

- Keys are the exponents (degrees of terms).
- Values are the corresponding coefficients.

For the above polynomial, we can efficiently represent it as:

```python
from sympy import *
from sympy.abc import x, y
p = Poly(x**100 + x**2 + 1)
print(p.as_dict()) # Outputs {(0,): 1, (2,): 1, (100,): 1}
```

Here, instead of storing 101 entries, we only store 3, making it significantly more
efficient.

For multivariate polynomials, we extend this idea by using tuples as keys to store the
exponents for each variable.
For example, consider:

$p(x, y)$ = x^{100}y + x^2y^2 + 1$

This has the sparse representation:

```python
p = Poly(x**100*y + x**2*y**2 + 1)
print(p.as_dict()) # Outputs {(0, 0): 1, (2, 2): 1, (100, 1): 1}
```

Here:
- `(0, 0)`: 1 represents the constant term .
- `(2, 2)`: 1 represents $x^2y^2$.
- `(100, 1)`: 1 represents $x^{100}y$.

---

## **Project Vision**

The problem here stated in the most succinct manner is that `SymPy` lacks fast multivariate polynomial support over ℤ/ℚ. Existing dense and sparse polynomial classes are Python-backed.

Hence...

## **The Goal**

Implement and integrate `FLINT` versions of the multivariate distributed polynomial rings (`PolyRing`), their sparse polynomial elements (`PolyElement`) and a new subclass of `DMP` to back multivariate polynomials with `FLINT` just like its univariate counterpart (`DUP_Flint`) already implemented and integrated. All of this over the integer(ℤ) and rational(ℚ) coefficent domain

Hence...

### **The Deliverables**

- **`FlintPolyRing`**: `FLINT`-backed sparse polynomial ring.

- **`FlintPolyElement`**: Sparse polynomial elements.

- **`DMP_Flint`**: `FLINT`-backed dense multivariate polynomials.


These will plug into the existing hierarchy, and activate only when `python-flint` and supported domains are detected on the user's side.

---

## Proof of concept

I have done some early work in the direction of this project and have coded some prototypes for the deliverables I enlisted above. Benchmarking from these prototypes have shown exciting results.

For instance, using `DMP.factor_list()` on a large polynomial:

On `master` (using `DMP_Python`):

```python
In [1]: p = Poly((x + y + z)**100*(x + y)*(y + z)*(z + x)) # A fairly huge polynomial

In [2]: type(p.rep)
Out[2]: sympy.polys.polyclasses.DMP_Python

In [3]: %time p.rep.factor_list()
CPU times: user 4min 44s, sys: 1.3 s, total: 4min 45s
Wall time: 4min 45s
Out[3]: 
(1,
 [(DMP_Python([[[1], [1, 0]]], ZZ), 1),
  (DMP_Python([[[1]], [[1], []]], ZZ), 1),
  (DMP_Python([[[1]], [[1, 0]]], ZZ), 1),
  (DMP_Python([[[1]], [[1], [1, 0]]], ZZ), 100)])
```

On my branch (with `DMP_Flint`):

```python
In [1]: p = Poly((x + y + z)**100*(x + y)*(y + z)*(z + x)) # A fairly huge polynomial

In [2]: type(p.rep)
Out[2]: sympy.polys.polyclasses.DMP_Flint

In [3]: %time p.rep.factor_list()
CPU times: user 12.7 ms, sys: 0 ms, total: 12.7 ms
Wall time: 14.3 ms
Out[3]: 
(1,
 [(DMP_Flint([[[1], [1, 0]]], ZZ), 1),
  (DMP_Flint([[[1]], [[1, 0]]], ZZ), 1),
  (DMP_Flint([[[1]], [[1], []]], ZZ), 1),
  (DMP_Flint([[[1]], [[1], [1, 0]]], ZZ), 100)])
```

That is about **~22000× times faster**! 

Although, my personal favourite speedup is in the `GCD` computation whose slowness is a significant and known issue in itself in `SymPy`.

Again, on `master` (using `DMP_Python`):

```python
In [1]: p = Poly((x + y + z)**100*(x + y)*(y + z)*(z + x)) # A fairly huge polynomial

In [2]: q = Poly((x + y)*(y + z)*(z + x))

In [3]: type(p.rep), type(q.rep) # Both DMP_Python objects
Out[3]: (sympy.polys.polyclasses.DMP_Python, sympy.polys.polyclasses.DMP_Python)

In [4]: %time p.rep.gcd(q.rep)
CPU times: user 800 ms, sys: 946 μs, total: 801 ms
Wall time: 800 ms
Out[4]: DMP_Python([[[1], [1, 0]], [[1], [2, 0], [1, 0, 0]], [[1, 0], [1, 0, 0], []]], ZZ)
```

Now on my branch (with `DMP_Flint`):

```python
In [1]: p = Poly((x + y + z)**100*(x + y)*(y + z)*(z + x)) # A fairly huge polynomial

In [2]: q = Poly((x + y)*(y + z)*(z + x))

In [3]: type(p.rep), type(q.rep) # Both DMP_Flint objects
Out[3]: (sympy.polys.polyclasses.DMP_Flint, sympy.polys.polyclasses.DMP_Flint)

In [4]: %time p.rep.gcd(q.rep)
CPU times: user 1.77 ms, sys: 0 ns, total: 1.77 ms
Wall time: 1.78 ms
Out[4]: DMP_Flint([[[1], [1, 0]], [[1], [2, 0], [1, 0, 0]], [[1, 0], [1, 0, 0], []]], ZZ, 2)
```

That's about **~450× times faster**!


Similarly, impressive speedups have been seen in the sparse case too:

Again, on `master` (using `PolyRing` and `PolyElement`):

```python
In [1]: R = PolyRing("x, y, z", ZZ, lex)

In [2]: x, y, z = R.gens

In [3]: p = (x + y + z)**100*(x + y)*(y + z)*(z + x) # A fairly huge polynomial

In [4]: q = (x + y)*(y + z)*(z + x)

In [5]: type(p), type(q) # Both PolyElement objects
Out[5]: (sympy.polys.rings.PolyElement, sympy.polys.rings.PolyElement)

In [6]: %time p.gcd(q)
CPU times: user 442 ms, sys: 1.28 ms, total: 444 ms
Wall time: 443 ms
Out[6]: x**2*y + x**2*z + x*y**2 + 2*x*y*z + x*z**2 + y**2*z + y*z**2
```

Now on my branch (with `FlintPolyRing` and `FlintPolyElement`):

```python
In [1]: R = FlintPolyRing("x, y, z", ZZ, lex)

In [2]: x, y, z = R.gens

In [3]: p = (x + y + z)**100*(x + y)*(y + z)*(z + x) # A fairly huge polynomial

In [4]: q = (x + y)*(y + z)*(z + x)

In [5]: type(p), type(q) # Both FlintPolyElement objects
Out[5]: (GSoC25.flintrings.FlintPolyElement, GSoC25.flintrings.FlintPolyElement)

In [6]: %time p.gcd(q)
CPU times: user 1.78 ms, sys: 0 ns, total: 1.78 ms
Wall time: 1.79 ms
Out[6]: x^2*y + x^2*z + x*y^2 + 2*x*y*z + x*z^2 + y^2*z + y*z^2
```

That's about **~250× times faster** !

---

## Closing Thoughts

This project is not just a performance upgrade. It’s about making `SymPy` reach **"state of the art"** speeds, backed by the best tools in computational algebra.

`SymPy` users won’t have to change a single line of code. But they’ll feel the difference.

With FLINT in the engine, I aim to give `SymPy` a jetpack it deserves.

---

For anyone who is interested in the code for the prototypes I made for the deliverables, or needs to understand this project in a greater detail than given in this blog, you can check out my project prosal :

📄 [Project Proposal](https://drive.google.com/file/d/1KCJmx-uxM1QvxGpM9gCsw7RvQyW2MPQg/view?usp=sharing)

---

> [^1]: The **`DMP`** (Dense Multivariate Polynomial) class returns a **`DUP_Flint`**  object when the coefficient domain is either integer (ℤ), rational (ℚ) or a finite field and the case is univariate. 

