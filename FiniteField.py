from galois import is_prime


def is_irreducible(p, fx_list):
    if len(fx_list) <= 3:  # Check only for degree 2 or 3
        for i in range(1, p):  # Start from 1 because 0 will always satisfy the polynomial
            total = sum(coef * (i ** idx) for idx, coef in enumerate(reversed(fx_list)))
            if total % p == 0:
                return False
    return True


def check_params(p, fx_list):
    if not is_prime(p):
        raise ValueError(f"p ({p}) must be prime")
    if fx_list[-1] == 0:
        raise ValueError(f"The value of the last coefficient cannot be zero: {fx_list}")
    if not is_irreducible(p, fx_list):
        raise ValueError(f"The polynomial {fx_list} is not irreducible for prime {p}")


class FiniteField:
    def __init__(self, p, fx_list):
        check_params(p, fx_list)
        self.p = p
        self.fx_list_coefficients = fx_list
        self.degree = len(fx_list) - 1
        self.elements = None  # To be generated lazily

    def __eq__(self, other):
        return self.p == other.p and self.fx_list_coefficients == other.fx_list_coefficients

    def __str__(self):
        return f"F[{self.p}]({self.fx_list_coefficients})"
