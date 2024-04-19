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
    if not isinstance(a, int) or not isinstance(b, int):
        raise ValueError("Both inputs must be integers.")

    # Additional check to ensure inputs are not fractional parts (should be redundant due to is_integer check)
    if (a % 1) != 0 or (b % 1) != 0:
        raise ValueError("Both inputs must be complete numbers (integers).")

    # Handle the edge case where both inputs are zero
    if a == 0 and b == 0:
        raise ValueError("At least one input must not be equal to 0.")

    # G.C.D calculation is done on the absolute values of the inputs and the signs are added back later on
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
    # return the gcd and the coefficients s and t
    return d, s, t
