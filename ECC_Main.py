class EllipticCurve:
    """An elliptic curve over a prime field.
    The field is specified by the parameter 'p'.
    The curve coefficients are 'a' and 'b'.
    The base point of the cyclic subgroup is 'g'.
    The order of the subgroup is 'n'.
    """

    def __init__(self, p, a, b, g, n):
        self.p = p
        self.a = a
        self.b = b
        self.g = g
        self.n = n

        assert pow(2, p - 1, p) == 1
        assert (4 * a * a * a + 27 * b * b) % p != 0
        assert self.lying_on_curve(g)
        assert self.multiply(n, g) is None

    def lying_on_curve(self, point):
        if point is None:
            return True

        x, y = point
        return (y * y - x * x * x - self.a * x - self.b) % self.p == 0

    def multiply(self, n, point):
        if n % self.n == 0 or point is None:
            return None

        if n < 0:
            return self.negative_point(self.multiply(-n, point))

        result = None
        addend = point

        while n:
            if n & 1:
                result = self.add(result, addend)
            addend = self.double_point(addend)
            n >>= 1

        return result

    def negative_point(self, point):
        if point is None:
            return None

        x, y = point
        result = x, -y % self.p
        assert self.lying_on_curve(result)
        return result

    def add(self, point1, point2):

        assert self.lying_on_curve(point1)
        assert self.lying_on_curve(point2)

        if point1 is None:
            return point2
        if point2 is None:
            return point1

        x1, y1 = point1
        x2, y2 = point2

        if x1 == x2 and y1 != y2:
            return None

        if x1 == x2:
            s = (3 * x1 * x1 + self.a) * inverse_mod(2 * y1, self.p)
        else:
            s = (y1 - y2) * inverse_mod(x1 - x2, self.p)

        x3 = s * s - x1 - x2
        y3 = s * (x3 - x1) + y1
        result = (x3 % self.p, -y3 % self.p)

        assert self.lying_on_curve(result)
        return result

    def double_point(self, point):
        return self.add(point, point)

    def __str__(self):
        a = abs(self.a)
        b = abs(self.b)
        a_sign = '-' if self.a < 0 else '+'
        b_sign = '-' if self.b < 0 else '+'

        return 'y^2 = (x^3 {} {}x {} {}) mod {}'.format(a_sign, a, b_sign, b, self.p)


def inverse_mod(n, p):
    if n == 0:
        raise ZeroDivisionError('Division by zero')
    if n < 0:
        return p - inverse_mod(-n, p)

    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, n

    while r != 0:
        quot = old_r // r
        old_r, r = r, old_r - r * quot
        old_s, s = s, old_s - s * quot
        old_t, t = t, old_s - t * quot

    gcd, x, y = old_r, old_s, old_t
    assert gcd == 1
    assert (n * x) % p == 1

    return x % p


Elliptic_curve = EllipticCurve(
    p=10177,
    a=1,
    b=-1,
    g=(1, 1),
    n=10331,
)
