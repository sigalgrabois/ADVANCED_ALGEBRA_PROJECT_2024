import numpy as np

from FiniteField import FiniteField
from PrimeFieldElement import PrimeFieldElement
import galois


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
    def __init__(self, l, a):
        """
        Initialize a finite field element.
        :param a: the given vector
        :param l: the finite field above which a is considered
        """
        self.l = l
        self.received_a = a
        # apply the modulu operation on the coefficients
        self.a = [coeff % l.p for coeff in a]
        if len(self.a) > l.f_x_degree:
            raise ValueError(f"Element's degree above the given field should not exceed {l.f_x_degree}")

        self.is_0 = all(
            coeff == 0 for coeff in self.a)  # all coefficients are non negative, their sum euqals 0 if they are all 0

        if self.is_0:
            self.matrix_representation = np.zeros((self.l.f_x_degree, self.l.f_x_degree), dtype=int)
        else:
            self.matrix_representation = self.calc_matrix_representation()

    def calc_matrix_representation(self):
        a, l = self.a, self.l
        n = l.f_x_degree
        p = l.p
        element_matrix_representation = np.zeros((n, n), dtype=int)

        # Initialization: first column equals the element itself
        temp_polynomial = np.zeros(n, dtype=int)
        temp_polynomial[:len(a)] = a  # Copy coefficients into the correct positions
        element_matrix_representation[:, 0] = temp_polynomial

        for i in range(1, n):
            # Multiplication in x
            shift_right_polynomial = np.roll(temp_polynomial, shift=1)
            if temp_polynomial[-1] != 0:
                congruent_equivalent = [(temp_polynomial[-1] * c) % p for c in l.congruate_equivalency]
                temp_polynomial = (congruent_equivalent + shift_right_polynomial) % p
            else:
                temp_polynomial = shift_right_polynomial

            element_matrix_representation[:, i] = temp_polynomial
        return element_matrix_representation

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
        if self.is_0:
            return FiniteFieldElement(self.l, [0] * self.l.f_x_degree)  # Return zero element

        # Calculate the inverse of the determinant of the divisor's matrix
        determinant_of_divider = np.linalg.det(other.matrix_representation)

        inverse_determinant = PrimeFieldElement(int(round(determinant_of_divider)), other.l.p).inverse()

        # Calculate the adjugate matrix and multiply it by the inverse determinant
        adjugate_matrix = np.round(np.linalg.inv(other.matrix_representation) * determinant_of_divider)

        inverse_matrix = np.mod(adjugate_matrix * inverse_determinant.a, other.l.p)

        # Perform matrix division and get the first column as the result
        division_matrix = self.matrix_representation @ inverse_matrix

        formed_element_matrix = np.mod(division_matrix, other.l.p)
        result = formed_element_matrix[:, 0]

        return FiniteFieldElement(other.l, result.tolist())

    def __str__(self):
        return " + ".join(f"{coeff} (mod {self.l.p})*x^{i}" for i, coeff in enumerate(self.a) if coeff != 0)

    def __repr__(self):
        coeffs_repr = ", ".join(f"{coeff} (mod {self.l.p})" for coeff in self.a)
        return f"FiniteFieldElement(FiniteField({self.l.p}, {self.l.f_x_original}), [{coeffs_repr}])"

    def __eq__(self, other):
        if isinstance(other, FiniteFieldElement):
            return self.a == other.a and self.l == other.l
        return False


if __name__ == '__main__':
    def run_section_4():
        print(f"====================================")
        print(f"section (4) - finite field element and matrix representation")
        print(f"====================================")

        # Define a Finite Field:
        p = 47  # prime number to set the field
        fx_coeff = [42, 3, 0, 1]  # a irreducible poly' coeff': for a_n*x^n+...+a_1*x+a_0 -> [a_0, a_1, ...]
        l = FiniteField(p, fx_coeff)  # the finite field object

        # Define a poly' by its coeff'
        a_coeff = [1, 2, 3]
        b_coeff = [1, 1, 1]
        a = FiniteFieldElement(l, a_coeff)  # an object of finite field element
        b = FiniteFieldElement(l, b_coeff)  # an object of finite field element

        print(f"polynomial a coeff' are:\n{a_coeff}")
        print(f"polynomial a in matrix representation:\n{a.matrix_representation}\n")  # [[1, 2, 3], [15, 39, 2], [10, 9, 39]]

        print(f"polynomial b coeff' are:\n{b_coeff}")
        print(f"polynomial b in matrix representation:\n{b.matrix_representation}")  # [[1, 1, 1], [5, 45  1], [5, 2, 45]]



