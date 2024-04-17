# Advanced Algebra: Finite Field Implementation (Course 83-804, 2024)

## Overview
This project implements arithmetic in finite fields of cardinality \(p^n\) for a given prime \(p\) and natural number \(n\). An \(n\)-dimensional extension of the prime field \(k = \mathbb{F}_p\) is constructed using an irreducible polynomial \(f(x)\) in \(k[x]\).

## Features
- Represents finite field elements as vectors, corresponding to polynomials of degree less than \(n\).
- Implements basic arithmetic over the prime field \(k\) and the extended field \(l\).
- Realizes non-zero elements in \(l\) as matrices in \(GL_n(k)\).
- Provides pretty-printing of finite field elements and implements arithmetic operations.



## Implementation
1. **Prime Field Arithmetic**: Define a class `PrimeFieldElement` with methods for addition, subtraction, multiplication, division, and inversion modulo \(p\) using the extended Euclidean algorithm.
2. **Finite Field Construction**: Create the `FiniteField` class to represent \(l\), taking parameters \(p\) and \(f(x)\), and validate their properties.
3. **Finite Field Elements**: Define the `FiniteFieldElement` class for elements in \(l\), linked to a specific instance of `FiniteField`.
4. **Matrix Representation**: Embed elements of \(l\) into \(GL_n(k)\), using the polynomial basis \(\{1, x, ..., x^{n-1}\}\).
5. **Arithmetic Operations**: Overload operators to support addition, subtraction, multiplication, and division in `FiniteFieldElement`.
6. **Exponentiation**: Implements efficient exponentiation by squaring for elements in \(k\) and \(l\).
7. **Multiplicative Order**: Add a method in `FiniteFieldElement` to calculate the multiplicative order of an element.
8. **Generator Finding**: Include functionality in `FiniteField` to find a generator for the multiplicative group \(l^\times\).
9. **BSGS Algorithm**: Implement the Baby-Step Giant-Step (BSGS) algorithm for the discrete logarithm problem in \(l\).
## NOTE:
In the implementation the finitie field element is in the form of [a0, a1, . . . , an−1] 
