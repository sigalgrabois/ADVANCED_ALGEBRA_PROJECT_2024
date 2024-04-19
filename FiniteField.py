import itertools
from typing import Set

import numpy as np
from galois import is_prime

from FiniteFieldElement import FiniteFieldElement
from PrimeFieldElement import PrimeFieldElement


def is_irreducible(p, f_x):
    """
     This function checks if a polynomial f(x) is irreducible over a prime field GF(p).
    For polynomials of degree > 3, it is assumed to be irreducible.
    Otherwise, it checks each element in the prime field to verify irreducibility.
    :param p: The prime number that defines the prime field GF(p).
    :param f_x: The polynomial f(x) to check for irreducibility.
    """

    if len(f_x) > 3:
        return True  # Assume irreducibility for polynomials of degree > 3

    # Check each element in the prime field to verify irreducibility
    for i in range(p):
        if sum(coeff * (i ** idx) for idx, coeff in enumerate(f_x)) % p == 0:
            return False  # Polynomial is reducible if a root is found

    return True  # Polynomial is irreducible if no roots are found


def check_params(p, fx_list):
    """
    This function checks the parameters for the FiniteField class.
    :param p: The prime number that defines the prime field GF(p).
    :param fx_list: The polynomial f(x) to check for irreducibility.
    """
    if not is_prime(p):
        raise ValueError(f"p ({p}) must be prime")
    if fx_list[0] == 0 or 2 >= len(fx_list):
        # for irreduciability, deg(f(x)) must be at least 2 and free coefficient must be non-zero
        raise ValueError(f"The value of the last coefficient cannot be zero: {fx_list}")
    if not is_irreducible(p, fx_list):
        raise ValueError(f"The polynomial {fx_list} is not irreducible for prime {p}")


def to_monic(p, f_x):
    """
    This function converts a polynomial f(x) to its monic representation in GF(p).
    :param p:
    :param f_x:
    :return:
    """
    leading_coeff = f_x[-1]
    if leading_coeff == 1:  # Already monic polynomial
        return f_x
    else:
        element = PrimeFieldElement(leading_coeff, p)
        inverse = element.inverse()
        inverse_coeff = inverse.coeff
        return [(coef * inverse_coeff) % p for coef in f_x]


class FiniteField:
    """
    This class represents a finite field (l) formed as an extension of a prime field, GF(p) by an irreducible
    polynomial f(x).
    The class assumes p is indeed prime and verifies irreducibility for polynomials whose degree is
    up to 3.
    It is assumed polynomial coefficients are arranged from free element (leftest coefficient) to the
    highest degree
    """

    def __init__(self, p, f_x):
        check_params(p, f_x)
        self.p = p  # prime whose corresponding field is the kernel of the described finite field
        self.f_x_original = f_x  # given irreducible polynomial whose highest degree coefficient may be larger than 1
        self.f_x_degree = len(
            f_x) - 1  # the degree of the given irreduciable polynomial which also equivalent to the field extension dimension
        self.f_x_monic = to_monic(p, f_x)  # monic representation of the given irreducible polynomial
        self.field_size = p ** self.f_x_degree  # number of elements above the described finite field

        # Calculate congregate equivalency;
        # the equivalence polynomial representation of x^(f_x_degree) deduced from the irreducible polynomial.
        lower_power_coeff_vec = np.array(
            self.f_x_monic[:-1])  # vector of f(x)_monic coefficients without the coefficient of X^(deg(f(x))
        # congurate equivalency: given polynomial
        # f(x)=a_0+a_1X+...+a_(n-1)X^(n-1)+X^n. perfoming modulu f(x)
        # means X^n=-a_0-a_1X-...-a_(n-1)X^(n-1)
        self.congruate_equivalency = (-lower_power_coeff_vec) % self.p

    def elements(self):
        """
        Generate all elements in the finite field, starting from the most significant coefficient.
        :return: A generator yielding all elements in the finite field.
        """
        # Generate all possible coefficients
        coefficients = range(self.p)

        # Generate all combinations of coefficients for polynomials up to degree f_x_degree - 1
        # We reverse the order of coefficients in each combination to prioritize higher degree terms
        for coeffs in itertools.product(coefficients, repeat=self.f_x_degree):
            # we used to reversed for the generator that in this case is left most significant
            reversed_coeffs = tuple(reversed(coeffs))  # Reverse the tuple of coefficients
            # Yield the corresponding FiniteFieldElement
            yield FiniteFieldElement(self, reversed_coeffs)

    def find_generator(self):
        """
        Find a generator of the multiplicative group of the finite field.
        :return: A generator of the multiplicative group.
        """
        e0_element = FiniteFieldElement(self, [0] * self.f_x_degree)
        e1_element = FiniteFieldElement(self, [1] + [0] * (self.f_x_degree - 1))
        elements = self.elements()
        for alpha in elements:
            if alpha == e0_element or alpha == e1_element:
                continue  # Skip the zero element
            order_alpha = alpha.multiplicative_order()
            if order_alpha == self.field_size - 1:
                return alpha
        raise ValueError("No generator found in the finite field.")

    def __eq__(self, other):
        return isinstance(other, FiniteField) and self.p == other.p and self.f_x_monic == other.f_x_monic

    def __repr__(self):
        return f'FiniteField({self.p}, {self.f_x_original})'

    def __str__(self):
        return f"F_{self.p}({self.f_x_original})"

    def __hash__(self):
        # Convert coefficients tuple to a hashable type (e.g., frozenset)
        return hash((self.p, tuple(self.f_x_original)))
