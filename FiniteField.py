from galois import is_prime


def is_irreducible(p, f_x):
    if len(f_x) > 3:
        return True  # Assume irreducibility for polynomials of degree > 3
    for i in range(p):
        if sum(coef * (i ** idx) for idx, coef in enumerate(reversed(f_x))) % p == 0:
            return False
    return True


def check_params(p, fx_list):
    if not is_prime(p):
        raise ValueError(f"p ({p}) must be prime")
    if fx_list[0] == 0:
        raise ValueError(f"The value of the last coefficient cannot be zero: {fx_list}")
    if not is_irreducible(p, fx_list):
        raise ValueError(f"The polynomial {fx_list} is not irreducible for prime {p}")


def to_monic(p, f_x):
    leading_coeff = f_x[-1]
    if leading_coeff == 1:
        return f_x
    inverse_leading = pow(leading_coeff, -1, p)
    return [coef * inverse_leading % p for coef in f_x]


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
        self.f_x_degree = len(f_x) - 1  # monic representation of the given irreduciable polynomial
        self.f_x_monic = to_monic(p,
                                  f_x)  # the degree of the given irreduciable polynomial which also equivalent to the field extension dimension
        self.field_size = p ** self.f_x_degree  # number of elements above the described finite field

        # Calculate congruate equivalency
        self.congruate_equivalency = [-self.f_x_monic[i] % p for i in range(
            self.f_x_degree)]  # the equivalence polynomial representation of x^(f_x_degree) deduced from the irreduciable polynomial

    def __eq__(self, other):
        return isinstance(other, FiniteField) and self.p == other.p and self.f_x_monic == other.f_x_monic

    def __str__(self):
        return f"FiniteField GF({self.p}^{self.f_x_degree}) with polynomial {self.f_x_monic}"
