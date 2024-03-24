from galois import is_prime


class PrimeFieldElement:
    # constructor
    def __init__(self, a, p):
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
            _, s, _ = extended_gcd(a, p)
            return PrimeFieldElement(s % p, p)


def gcd(a, b):
    while b != 0:
        a, b = b, a % b
    return a


def extended_gcd(a, b):
    if b == 0:
        return a, 1, 0
    else:
        d, s, t = extended_gcd(b, a % b)
        return d, t, s - (a // b) * t
