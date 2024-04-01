import numpy as np
from galois import is_prime


def gcd(a, b):
    """
    This function calculates the greatest common divisor of two integers a and b

    :param a:
    :param b:
    :return:
    """
    while b != 0:
        a, b = b, a % b
    return a


def xgcd(a, b):
    """
    This function performs the extended Euclides algorithm which receives two
    complete numbers a and b and returns their g.c.d, (a,b)=d and two complete
    numbers s,t which uphold: d=s*a+t*b
    """

    # Sanity checks
    if not a.is_integer() or not b.is_integer():
        raise ValueError("Both inputs must be integers.")

    if (a % 1) != 0 or (b % 1) != 0:
        raise ValueError("both inputs must be complete numbers.")

    # Edge cases
    if a == 0 and b == 0:
        raise ValueError("At least one input must not be equal to 0.")

    # The GCD calculation is done for the absolute values of a, b
    mag_a, mag_b = abs(a), abs(b)

    if a == 0 or b == 0:
        d = max(mag_a, mag_b)
        s = 1
        t = 1
    else:
        # Iterative process
        earlier_residue, early_residue = mag_a, mag_b
        temp_residue = 1  # Default value to ensure entering the loop
        s, t = 0, 0
        last_s, last_t = 1, 0
        current_s, current_t = 0, 1

        while early_residue != 0:
            quotient = earlier_residue // early_residue
            earlier_residue, early_residue = early_residue, earlier_residue - quotient * early_residue
            last_s, current_s = current_s, last_s - quotient * current_s
            last_t, current_t = current_t, last_t - quotient * current_t

        d = earlier_residue
        s = last_s
        t = last_t

    # Adjust the signs of s and t
    s = s if a >= 0 else -s
    t = t if b >= 0 else -t

    # Sanity check
    if d <= 0 or d != s * a + t * b:
        raise ValueError("Something went wrong with the calculations.")

    return d, s, t


def calc_matrix_representation_above_finite_field(a, l):
    n = l.f_x_degree  # extension dimension of the finite field
    p = l.p  # the prime used as the kernel of the extended field
    element_matrix_representation = np.zeros((n, n), dtype=int)

    # Initialization: first column equals the element itself
    temp_polynomial = np.array(a + [0] * (n - len(a)))  # Extend a to match the field's degree if necessary
    element_matrix_representation[:, 0] = temp_polynomial

    for i in range(1, n):
        # Multiplication in x
        shift_right_polynomial = np.roll(temp_polynomial, shift=1)
        shift_right_polynomial[0] = 0  # Set the first element to 0 after the shift

        # If the highest degree term is nonzero, calculate the congruent equivalent
        if temp_polynomial[-1] != 0:
            congruent_equivalent = [(temp_polynomial[-1] * c) % p for c in l.congruate_equivalency]
            temp_polynomial = (congruent_equivalent + shift_right_polynomial) % p
        else:
            temp_polynomial = shift_right_polynomial

        element_matrix_representation[:, i] = temp_polynomial

    return element_matrix_representation
