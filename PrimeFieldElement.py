
class PrimeFieldElement:
    def __init__(self, a, p):
        self.a = a % p
        self.p = p

    def __add__(self, other):
        if isinstance(other, PrimeFieldElement) and self.p == other.p:
            return PrimeFieldElement((self.a + other.a) % self.p, self.p)
        else:
            raise ValueError("Cannot perform addition with different prime fields")

    def __sub__(self, other):
        if isinstance(other, PrimeFieldElement) and self.p == other.p:
            return PrimeFieldElement((self.a - other.a) % self.p, self.p)
        else:
            raise ValueError("Cannot perform subtraction with different prime fields")

    def __mul__(self, other):
        if isinstance(other, PrimeFieldElement) and self.p == other.p:
            return PrimeFieldElement((self.a * other.a) % self.p, self.p)
        else:
            raise ValueError("Cannot perform multiplication with different prime fields")

    def __truediv__(self, other):
        if isinstance(other, PrimeFieldElement) and self.p == other.p:
            if other.a == 0:
                raise ZeroDivisionError("Cannot divide by zero")
            else:
                return self * other.inverse()
        else:
            raise ValueError("Cannot perform division with different prime fields")

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