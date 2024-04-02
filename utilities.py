import numpy as np
from galois import is_prime

import PrimeFieldElement


def xgcd(a, b):
    """
    Performs the Extended Euclidean Algorithm which calculates the greatest common divisor (gcd)
    of two integers a and b, along with integers s and t such that d = s * a + t * b.

    Args:
        a (int): First input integer.
        b (int): Second input integer.

    Returns:
        d (int): The greatest common divisor of a and b.
        s (int): Coefficient such that d = s * a + t * b.
        t (int): Coefficient such that d = s * a + t * b.

    Raises:
        ValueError: If either a or b is not an integer or both are zero.

    Usage example :
    d, s, t = xgcd(234, 61)
    Output: 1 6 -23
    """

    # Sanity checks to ensure that inputs are integers
    if not a.is_integer() or not b.is_integer():
        raise ValueError("Both inputs must be integers.")

    # Additional check to ensure inputs are not fractional parts (should be redundant due to is_integer check)
    if (a % 1) != 0 or (b % 1) != 0:
        raise ValueError("Both inputs must be complete numbers (integers).")

    # Handle the edge case where both inputs are zero
    if a == 0 and b == 0:
        raise ValueError("At least one input must not be equal to 0.")

    # Calculate the GCD using the absolute values of a and b to ensure positive results
    mag_a, mag_b = abs(a), abs(b)

    # Handle the case where one of the inputs is zero
    if a == 0 or b == 0:
        d = max(mag_a, mag_b)
        s = 1
        t = 1
    else:
        # Iterative process using the Extended Euclidean Algorithm
        earlier_residue, early_residue = mag_a, mag_b
        last_s, last_t = 1, 0  # Initializing previous coefficients
        current_s, current_t = 0, 1  # Initializing current coefficients

        # Loop until the remainder becomes zero
        while early_residue != 0:
            quotient = earlier_residue // early_residue
            earlier_residue, early_residue = early_residue, earlier_residue - quotient * early_residue
            last_s, current_s = current_s, last_s - quotient * current_s
            last_t, current_t = current_t, last_t - quotient * current_t

        # The last non-zero remainder is the gcd
        d = earlier_residue
        s = last_s
        t = last_t

    # Adjust the signs of s and t if a or b were negative
    s = s if a >= 0 else -s
    t = t if b >= 0 else -t

    # Perform a final sanity check to ensure the correctness of the algorithm
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
