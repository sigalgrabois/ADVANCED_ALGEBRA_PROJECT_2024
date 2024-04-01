import numpy as np
from utilities import calc_matrix_representation_above_finite_field

def polynomial_mod(poly, mod_poly, p):
    """
    Reduce a polynomial by another polynomial in a finite field of order p.

    :param poly: The polynomial to be reduced, given as a list of coefficients.
                 Coefficients are ordered from the lowest degree to the highest.
    :param mod_poly: The modulus polynomial, given as a list of coefficients.
    :param p: The prime order of the finite field.
    :return: The remainder of the division of poly by mod_poly in GF(p).
    """
    # Convert lists to NumPy polynomial objects
    poly = np.poly1d(list(reversed(poly)))
    mod_poly = np.poly1d(list(reversed(mod_poly)))

    # Perform polynomial division in the field
    quotient, remainder = np.polydiv(poly, mod_poly)

    # The coefficients should be reduced modulo p
    remainder_coeffs = np.mod(remainder.coeffs, p)

    # Convert back to the standard list format (lowest degree to highest)
    remainder_list = list(reversed(remainder_coeffs))

    # Ensure the list has the correct length
    deg_diff = len(mod_poly.coeffs) - 1 - len(remainder_list)
    if deg_diff > 0:
        remainder_list.extend([0] * deg_diff)

    return remainder_list


class FiniteFieldElement:
    def __init__(self, finite_field, coeffs_vector):
        if len(coeffs_vector) > finite_field.f_x_degree:
            raise ValueError(f"Element's degree above the given field should not exceed {finite_field.f_x_degree}")

        self.l = finite_field
        self.received_a = coeffs_vector
        self.a = [coeff % self.l.p for coeff in coeffs_vector]

        self.is_0 = all(coeff == 0 for coeff in self.a)

        if self.is_0:
            self.matrix_representation = np.zeros((self.l.f_x_degree, self.l.f_x_degree), dtype=int)
        else:
            self.matrix_representation = self.calc_matrix_representation()

    def calc_matrix_representation(self):
        # Assuming calcMatrixRepresentationaboveFiniteField is already defined
        return calc_matrix_representation_above_finite_field(self.a, self.l)

    def pretty_printing(self):
        print("Polynomial Representation:")
        poly_repr = '+'.join(f"{coeff}x^{i}" if i > 1 else (f"{coeff}x" if i == 1 else str(coeff))
                             for i, coeff in enumerate(self.a) if coeff != 0)
        print(poly_repr or '0')

        print("Vector Representation:")
        print(self.a)

        print("Matrix Representation:")
        print(self.matrix_representation)

    def __add__(self, other):
        if self.l != other.l:
            raise ValueError("Both elements must be above the same field")
        sum_coeffs = [(x + y) % self.l.p for x, y in zip(self.a, other.a)]
        return FiniteFieldElement(self.l, sum_coeffs)

    def __sub__(self, other):
        if self.l != other.l:
            raise ValueError("Both elements must be above the same field")
        sub_coeffs = [(x - y) % self.l.p for x, y in zip(self.a, other.a)]
        return FiniteFieldElement(self.l, sub_coeffs)

    def __mul__(self, other):
        if self.l != other.l:
            raise ValueError("Both elements must be above the same field")

        # Direct polynomial multiplication followed by reduction (for simplicity here)
        # In practice, you need to ensure this respects finite field rules
        result_coeffs = np.polynomial.polynomial.polymul(self.a, other.a)
        result_coeffs = np.mod(result_coeffs, self.l.p)  # Reduce each coefficient modulo p

        # Reduce the result by the field's polynomial if necessary
        # (Assuming you have a method to do polynomial division and find the remainder)
        result_coeffs = polynomial_mod(result_coeffs, self.l.f_x_original, self.l.p)

        return FiniteFieldElement(self.l, result_coeffs)

    def __truediv__(self, other):
        if self.l != other.l:
            raise ValueError("Both elements must be above the same field")
        if other.is_0:
            raise ZeroDivisionError("Division by zero element in finite field not defined")

        # Division implementation here should consider field's properties
        # Simplified version for demonstration, assuming matrix inverse is defined elsewhere
        inverse_matrix = np.linalg.inv(other.matrix_representation)  # This is not generally correct in finite fields
        div_matrix = self.matrix_representation @ inverse_matrix
        return FiniteFieldElement(self.l, (div_matrix % self.l.p)[0])
    def __str__(self):
        return " + ".join(f"{coeff} (mod {self.l.p})*x^{i}" for i, coeff in enumerate(self.a) if coeff != 0)

