import numpy as np
import galois

import FiniteField
import PrimeFieldElement
from sympy import Matrix


# TODO:
# 1. delete the extra printings in the code below
# 3. check consistency of the code - it brings back an element.

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
            shift_right_polynomial = np.concatenate(([0], temp_polynomial[:n - 1]))  # np.roll(temp_polynomial, shift=1)
            if temp_polynomial[-1] != 0:
                congruent_equivalent = np.mod(temp_polynomial[-1] * l.congruate_equivalency,
                                              p)  # [(temp_polynomial[-1] * c) % p for c in l.congruate_equivalency]
                temp_polynomial = (congruent_equivalent + shift_right_polynomial) % p
            else:
                temp_polynomial = shift_right_polynomial

            element_matrix_representation[:, i] = temp_polynomial
        return element_matrix_representation

    def poly_print(self):
        # Generate polynomial string representation
        poly_repr = ' + '.join(
            f"{coeff}x^{i}" if coeff != 0 and i > 0 else f"{coeff}"
            for i, coeff in enumerate(self.a)
            if coeff != 0
        ) or '0'
        return f"Polynomial Representation: {poly_repr}"

    def vec_print(self):
        print(self.a)

    def matrix_print(self):
        for row in self.matrix_representation:
            print(' '.join(str(x) for x in row))

    def pretty_printing(self):
        self.poly_print()
        self.vec_print()
        self.matrix_print()

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
        result_matrix = np.matmul(self.matrix_representation,other.matrix_representation)
        coeffs = result_matrix[:, 0].tolist()
        coeffs_res = [coeff % self.l.p for coeff in coeffs]
        return FiniteFieldElement(self.l, coeffs_res)

    def __truediv__(self, other):
        if self.l != other.l:
            raise ValueError("Both elements must be above the same field")
        if other.is_0:
            raise ZeroDivisionError("Division by zero element in finite field not defined")
        if self.is_0:
            return FiniteFieldElement(self.l, [0] * self.l.f_x_degree)  # Return zero element

        # Calculate the inverse of the determinant of the divisor's matrix
        determinant_of_divider = np.linalg.det(other.matrix_representation)

        inverse_determinant = PrimeFieldElement.PrimeFieldElement(int(round(determinant_of_divider)),
                                                                  other.l.p).inverse()

        # Calculate the adjugate matrix and multiply it by the inverse determinant
        # change the matrix to sympy matrix
        adjugate_matrix = Matrix(other.matrix_representation).adjugate()

        inverse_matrix = np.mod(adjugate_matrix * inverse_determinant.a, other.l.p)

        # Perform matrix division and get the first column as the result
        division_matrix = self.matrix_representation @ inverse_matrix

        formed_element_matrix = np.mod(division_matrix, other.l.p)
        result = formed_element_matrix[:, 0]

        return FiniteFieldElement(other.l, result.tolist())

    def __pow__(self, exponent):
        e1_element = FiniteFieldElement(self.l, [1] + [0] * (self.l.f_x_degree - 1))
        base = self
        if exponent == 0:
            # Anything raised to the power of 0 is 1
            return e1_element
        elif exponent == 1:
            # Return the element itself
            return base
        elif exponent < 0:
            # Compute the inverse and exponentiate with the absolute value of the exponent
            base = e1_element / base
            exponent = -1 * exponent
        # Apply exponentiation by squaring
        result = e1_element
        while exponent > 0:
            if exponent % 2 == 1:
                result *= base  # Multiply result by the base if the current bit of the exponent is 1
            base *= base  # Square the base
            exponent //= 2  # Shift exponent to the right by 1 bit
        return result

    def multiplicative_order(self):
        """
        Compute the multiplicative order of the element.
        :return: The multiplicative order of the element.
        """
        if self.is_0:
            raise ValueError("Multiplicative order is not defined for zero element")

        power = 1
        curr_element = self
        e1_element = FiniteFieldElement(self.l, [1] + [0] * (self.l.f_x_degree - 1))
        while curr_element != e1_element:
            curr_element *= self  # Multiply by self
            power += 1
        return power

    def __str__(self):
        return self.poly_print()

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
        l = FiniteField.FiniteField(p, fx_coeff)  # the finite field object

        # Define a poly' by its coeff'
        a_coeff = [1, 2, 3]
        b_coeff = [1, 1, 1]
        a = FiniteFieldElement(l, a_coeff)  # an object of finite field element
        b = FiniteFieldElement(l, b_coeff)  # an object of finite field element

        print(f"polynomial a coeff' are:\n{a_coeff}")
        print(
            f"polynomial a in matrix representation:\n{a.matrix_representation}\n")  # [[1, 2, 3], [15, 39, 2], [10, 9, 39]]

        print(f"polynomial b coeff' are:\n{b_coeff}")
        print(
            f"polynomial b in matrix representation:\n{b.matrix_representation}")  # [[1, 1, 1], [5, 45  1], [5, 2, 45]]

        p = 5
        f_x = [2, 1, 1]
        field = FiniteField.FiniteField(p, f_x)
        a_coeff = [3, 2]
        b_coeff = [1, 3]
        a = FiniteFieldElement(field, a_coeff)
        b = FiniteFieldElement(field, b_coeff)

        print(a / b)
        print(a * b)


    run_section_4()
