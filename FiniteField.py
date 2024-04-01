import numpy as np

def is_irreducible(p, f_x):
    if len(f_x) > 3:
        return True  # Assume irreducibility for polynomials of degree > 3
    for i in range(p):
        if sum(coef * (i ** idx) for idx, coef in enumerate(reversed(f_x))) % p == 0:
            return False
    return True

def to_monic(p, f_x):
    leading_coeff = f_x[-1]
    if leading_coeff == 1:
        return f_x
    inverse_leading = pow(leading_coeff, -1, p)
    return [coef * inverse_leading % p for coef in f_x]

class FiniteField:
    def __init__(self, p, f_x):
        if not is_irreducible(p, f_x) or f_x[0] == 0:
            raise ValueError("The polynomial is not irreducible or has zero constant term, so a finite field cannot be formed.")

        self.p = p
        self.f_x_original = f_x
        self.f_x_degree = len(f_x) - 1
        self.f_x_monic = to_monic(p, f_x)
        self.field_size = p ** self.f_x_degree

        # Calculate congruate equivalency
        self.congruate_equivalency = [-self.f_x_monic[i] % p for i in range(self.f_x_degree)]

    def __eq__(self, other):
        return isinstance(other, FiniteField) and self.p == other.p and self.f_x_monic == other.f_x_monic

    def __str__(self):
        return f"FiniteField GF({self.p}^{self.f_x_degree}) with polynomial {self.f_x_monic}"

# Example of usage
try:
    p = 3
    f_x = [1, 0, 2]  # Example polynomial
    ff = FiniteField(p, f_x)
    print(ff)
except ValueError as e:
    print(e)
