#                                                     Advanced Algebra: Finite Field Implementation (Course 83-804, 2024)

 ###                                              By Ayelet Aharon (200641058), Erel Naor (206913600) and Sigal Grabois(319009304)

## Overview
This project implements arithmetic in finite fields of cardinality \(p^n\) for a given prime \(p\) and natural number \(n\). An \(n\)-dimensional extension of the prime field \(k = \mathbb{F}_p\) is constructed using an irreducible polynomial \(f(x)\) in \(k[x]\).

## Features
- **Vector Representation**: Represents elements of finite fields as vectors, corresponding to polynomials with degrees less than \(n\).
- **Comprehensive Arithmetic**: Implements both basic and complex arithmetic operations over the prime field \(k\) and its extension \(l\).
- **Matrix Realization**: Realizes non-zero elements in \(l\) as matrices within \(GL_n(k)\), which is crucial for performing complex algebraic computations.
- **Enhanced Display**: Pretty-prints finite field elements in various formats, including polynomials and matrices, to facilitate easier debugging and understanding.

## Implementation
1. **Prime Field Arithmetic**: The `PrimeFieldElement` class manages operations such as addition, subtraction, multiplication, division, and inversion modulo \(p\), utilizing the extended Euclidean algorithm for modular inversion.
2. **Field Construction**: The `FiniteField` class initializes field \(l\) with parameters \(p\) and \(f(x)\), ensuring the polynomial's irreducibility for valid field construction.
3. **Element Management**: The `FiniteFieldElement` class allows manipulation of elements in \(l\), with each linked to a specific instance of `FiniteField`.
4. **Matrix Representation**: Implements an embedding of \(l\) into \(GL_n(k)\) using a polynomial basis, which is pivotal for complex algebraic operations.
5. **Operator Overloading**: Supports intuitive arithmetic operations within `FiniteFieldElement` through operator overloading.
6. **Efficient Exponentiation**: Implements exponentiation by squaring for elements in \(k\) and \(l\), optimizing computational efficiency.
7. **Order Calculation**: Adds functionality to determine the multiplicative order of elements in \(l^\times\).
8. **Generator Identification**: Includes a method in `FiniteField` to find generators for the group \(l^\times\), which is essential for constructing cyclic groups.
9. **BSGS Algorithm**: Features the Baby-Step Giant-Step algorithm for addressing the discrete logarithm problem in \(l\), enhancing the security analysis.

## Runnig the project
* To run the tests for different sections of the project, you can use the `tests.py` script: <br>
`python tests.py`
This script will automatically run all predefined unit tests, verifying the correctness of each module.
* For a practical demonstration of the project, particularly the BSGS algorithm, use the `orchestrator.py` script: <br>
  `python orchestrator.py`
  #### Note: Ensure that Python is installed on your system and that all dependencies specified in the project's requirements.txt are installed: <br> `pip install -r requirements.txt`



## Additional Information
- Each element in \(l\) is represented as a vector [a0, a1, ..., an−1], simplifying polynomial operations.

