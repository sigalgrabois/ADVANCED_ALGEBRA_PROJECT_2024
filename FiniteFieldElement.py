import numpy as np
import FiniteField
import PrimeFieldElement
from sympy import Matrix


class FiniteFieldElement:
    """
    This class represents a finite field element. Assuming the finite field, l, is an 'n' dimension extension field of a prime field GF(p),
    the constructor is given l and a vector of length n which represents an element on the field
    NOTES:
    1. It is assumed given vector is of the form [a_0,a_1,...,a_n-1] which equivalent to the polynomial
    a_0+a_1X+,...,a_n-1X^(n-1)
    2. Even if the given coefficients {a_0,...,a_n-1} exceeds the range [0,p-1], they would be
    translated to this range via mod p operation
    """

    def __init__(self, l, a):
        """
        Initialize a finite field element.
        :param a: the given vector in form of [a_0, a_1, ..., a_n] where a_i is the coefficient of x^i
        :param l: the finite field above which a is considered
        """
        self.l = l
        self.received_a = a
        # apply the modulu operation on the coefficients
        self.a = [coeff % l.p for coeff in a]
        if len(self.a) > l.f_x_degree:
            raise ValueError(f"Element's degree above the given field should not exceed {l.f_x_degree}")

        # a boolean variable indicating whether the given a is the 0 element of the field
        self.is_0 = all(
            coeff == 0 for coeff in self.a)  # all coefficients are non-negative, their sum equals 0 if they are all 0

        # representation of the given element a as a matrix above GLn(GF(p))
        if self.is_0:
            self.matrix_representation = np.zeros((self.l.f_x_degree, self.l.f_x_degree), dtype=int)
        else:

            self.matrix_representation = self.calc_matrix_representation()

    def calc_matrix_representation(self):
        """
        The function calculates a matrix representation of the given element, a,
        which is formed by the elementary basis {1,x,...,x^(n-1)} and the linear transformation x^i->a*x^i s.t.
        element_matrix_representation = [a, aX,..., aX^(n-1)]
        :return: matrix representation of the element
        """
        a, l = self.a, self.l
        n = l.f_x_degree  # extension dimension of the finite field
        p = l.p
        element_matrix_representation = np.zeros((n, n), dtype=int)  # would hold the matrix representation of a

        # Initialization: first column equals the element itself
        temp_polynomial = np.zeros(n, dtype=int)
        temp_polynomial[:len(a)] = a  # Copy coefficients into the correct positions
        element_matrix_representation[:, 0] = temp_polynomial

        for i in range(1, n):
            # Multiplication in x: if the highest degree of 'temp_polynomial' (the last element) is zero,
            # then a cyclic shift one place right of it, is sufficient to represent multiplication.
            # Otherwise, the result should be the sum of the congruate equivalent of the highest degree
            # and the shift right of the remaining elements of 'temp_polynomial'

            shift_right_polynomial = np.concatenate(([0], temp_polynomial[:n - 1]))
            if temp_polynomial[-1] != 0:  # then the congruent equivalent of the highest degree should be calculated
                congruent_equivalent = np.mod(temp_polynomial[-1] * l.congruate_equivalency,
                                              p)
                temp_polynomial = (congruent_equivalent + shift_right_polynomial) % p  # mod p is applied to ensure
                # coefficients are above GF(p)
            else:
                # cyclic shift right is sufficient
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
        return f"{poly_repr}"

    def vec_print(self):
        print(self.a)

    def matrix_print(self):
        for row in self.matrix_representation:
            print(' '.join(str(x) for x in row))

    def pretty_printing(self):
        print("polynomial representation:")
        print(self.poly_print())
        print("vector representation:")
        self.vec_print()
        print("matrix representation:")
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
        result_matrix = np.matmul(self.matrix_representation, other.matrix_representation)
        coeffs = result_matrix[:, 0].tolist()
        coeffs_res = [coeff % self.l.p for coeff in coeffs]
        return FiniteFieldElement(self.l, coeffs_res)

    def __truediv__(self, other):
        """
        Perform division of two elements in the finite field.
        :param other: The element to divide by.
        :return: The result of the division.

        it is desired to calculate multiplication by 'obj.matrix_representation*inverse(
        other_obj.matrix_representation)'. Only that matrix inversion may result in a fractional representation which
        a simple modulo wouldn't resolve. To combat this we rely on the formula: A^(-1)=det(A)^(-1)*Adjugate(A) The
        result of Adjugate(A) is assembled with integers (except perhaps a negligible deviation resulting from
        calculation is performed, hence round is used)

        """
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
        # print the object in a readable format

        return f"FiniteFieldElement(FiniteField(p:{self.l.p}, fx:{self.l.f_x_original}), a:{self.a})"

    def __repr__(self):
        coeffs_repr = ", ".join(f"{coeff} (mod {self.l.p})" for coeff in self.a)
        return f"FiniteFieldElement(FiniteField({self.l.p}, {self.l.f_x_original}), [{coeffs_repr}])"

    def __eq__(self, other):
        if isinstance(other, FiniteFieldElement):
            return self.a == other.a and self.l == other.l
        return False
