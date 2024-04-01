from galois import is_prime
from utilities import *


class PrimeFieldElement:
    """
        This class represents a prime field element. For construction, it
        receives two integers (a,p), where p is prime number and a is an element in
        GF(p). NOTE: GF(p):={0,1,...p-1}, hence any received 'a' would be translated to this set of values by appling modulu p. The class provides implementation of
        basic arithmetic operations above GF(p).''''''
    """

    # constructor
    def __init__(self, a: int, p: int):
        if not is_prime(p):
            raise ValueError(f"p ({p}) is not prime")
        self.a = a % p
        self.p = p

    def __add__(self, other):
        if isinstance(other, PrimeFieldElement) and self.p == other.p:
            return PrimeFieldElement((self.a + other.a), self.p)
        else:
            raise ValueError("Cannot perform addition of elements with different prime fields")

    def __sub__(self, other):
        if isinstance(other, PrimeFieldElement) and self.p == other.p:
            return PrimeFieldElement((self.a - other.a), self.p)
        else:
            raise ValueError("Cannot perform subtraction with different prime fields")

    def __mul__(self, other):
        if isinstance(other, PrimeFieldElement) and self.p == other.p:
            return PrimeFieldElement((self.a * other.a), self.p)
        else:
            raise ValueError("Cannot perform multiplication with different prime fields")

    def __truediv__(self, other):
        if isinstance(other, PrimeFieldElement) and self.p == other.p:
            if other.a == 0:
                raise ValueError("Cannot divide by zero")
            else:
                return self * other.inverse()
        else:
            raise ValueError("Cannot perform division with different prime fields")

    def __eq__(self, other):
        if isinstance(other, PrimeFieldElement):
            return self.a == other.a and self.p == other.p
        return False

    def __pow__(self, power, modulo=None):
        return PrimeFieldElement(pow(self.a, power, self.p), self.p)

    def inverse(self):
        a, p = self.a, self.p
        if a == 0:
            raise ZeroDivisionError("Cannot compute inverse of zero")
        elif gcd(a, p) != 1:
            raise ValueError("Element does not have an inverse in the prime field")
        else:
            d, s, t = xgcd(a, p)
            return PrimeFieldElement(s, p)
