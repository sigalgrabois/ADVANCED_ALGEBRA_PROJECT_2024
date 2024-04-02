import numpy as np
from galois import is_prime

from PrimeFieldElement import PrimeFieldElement


def is_irreducible(p, f_x):
    # This function checks if a polynomial f(x) is irreducible over a prime field GF(p).
    # For polynomials of degree > 3, it is assumed to be irreducible.
    # Otherwise, it checks each element in the prime field to verify irreducibility.

    if len(f_x) > 3:
        return True  # Assume irreducibility for polynomials of degree > 3

    # Check each element in the prime field to verify irreducibility
    for i in range(p):
        if sum(coef * (i ** idx) for idx, coef in enumerate(reversed(f_x))) % p == 0:
            return False  # Polynomial is reducible if a root is found

    return True  # Polynomial is irreducible if no roots are found


def check_params(p, fx_list):
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
    This class represents a finite field formed as an extension of a prime field, GF(p) by an irreduciable polynomial f(x).
    The class assumes p is indeed prime and verifies irreduciability for
    polynomials whose degree is up to 3. It is assumed polynomial
    coefficients are arranged from free element (leftest coefficient) to the highest degree
    """

    def __init__(self, p, f_x):
        check_params(p, f_x)
        self.p = p  # prime whose corresponding field is the kernel of the described finite field
        self.f_x_original = f_x  # given irreduciable polynomial whose highest degree coefficient may be larger than 1
        self.f_x_degree = len(
            f_x) - 1  # the degree of the given irreduciable polynomial which also equivalent to the field extension dimension
        self.f_x_monic = to_monic(p,
                                  f_x)  # monic representation of the given irreduciable polynomial
        self.field_size = p ** self.f_x_degree  # number of elements above the described finite field

        # Calculate congruate equivalency
        lower_power_coeff_vec = np.array(self.f_x_monic[:-1])  # coefficients of x^(f_x_degree) in the monic polynomial
        # self.congruate_equivalency = np.mod(-lower_power_coeff_vec, self.p)
        self.congruate_equivalency = (-lower_power_coeff_vec) % self.p

    def __eq__(self, other):
        return isinstance(other, FiniteField) and self.p == other.p and self.f_x_monic == other.f_x_monic

    def __repr__(self):
        return f'FiniteField({self.p}, {self.f_x_original})'

    def __str__(self):
        return f"F_{self.p}({self.f_x_original})"
