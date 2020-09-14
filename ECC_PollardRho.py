import random

from ECC_Main import inverse_mod, Elliptic_curve as E_curve


class PollardRho:
    def __init__(self, point1, point2):
        self.point1 = point1
        self.point2 = point2
        self.a1_add = random.randrange(1, E_curve.n)
        self.b1_add = random.randrange(1, E_curve.n)
        self.x1_add = E_curve.add(E_curve.multiply(self.a1_add, point1), E_curve.multiply(self.b1_add, point2))

        self.a2_add = random.randrange(1, E_curve.n)
        self.b2_add = random.randrange(1, E_curve.n)
        self.x2_add = E_curve.add(E_curve.multiply(self.a2_add, point1), E_curve.multiply(self.b2_add, point2))

    def __iter__(self):
        window = E_curve.n // 3 + 1
        x = None
        a = 0
        b = 0

        while True:
            if x is None:
                i = 0
            else:
                i = x[0] // window

            if i == 0:
                a += self.a1_add
                b += self.b1_add
                x = E_curve.add(x, self.x1_add)
            elif i == 1:
                a *= 2
                b *= 2
                x = E_curve.double_point(x)
            elif i == 2:
                a += self.a2_add
                b += self.b2_add
                x = E_curve.add(x, self.x2_add)
            else:
                raise AssertionError(i)

            a = a % E_curve.n
            b = b % E_curve.n

            yield x, a, b


def log(p, q, counter=None):
    assert E_curve.lying_on_curve(p)
    assert E_curve.lying_on_curve(q)

    # Pollard's Rho may fail sometimes: it may find a1 == a2 and b1 == b2,
    # leading to a division by zero error. Because PollardRhoSequence uses
    # random coefficients, we have more chances of finding the logarithm
    # if we try again, without affecting the asymptotic time complexity.
    # We try at most three times before giving up.
    for i in range(3):
        sequence = PollardRho(p, q)

        tortoise = iter(sequence)
        hare = iter(sequence)

        # The range is from 0 to curve.n - 1, but actually the algorithm will
        # stop much sooner (either finding the logarithm, or failing with a
        # division by zero).
        for j in range(E_curve.n):
            x1, a1, b1 = next(tortoise)

            x2, a2, b2 = next(hare)
            x2, a2, b2 = next(hare)

            if x1 == x2:
                if b1 == b2:
                    # This would lead to a division by zero. Try with
                    # another random sequence.
                    break

                x = (a1 - a2) * inverse_mod(b2 - b1, E_curve.n)
                logarithm = x % E_curve.n
                steps = i * E_curve.n + j + 1
                return logarithm, steps

    raise AssertionError('logarithm not found')


def main():
    x = random.randrange(1, E_curve.n)
    p = E_curve.g
    q = E_curve.multiply(x, p)

    print('Curve: {}'.format(E_curve))
    print('Curve order: {}'.format(E_curve.n))
    print('p = (0x{:x}, 0x{:x})'.format(*p))
    print('q = (0x{:x}, 0x{:x})'.format(*q))
    print(x, '* p = q')

    y, steps = log(p, q)
    print('log(p, q) =', y)
    print('Took', steps, 'steps')

    assert x == y


if __name__ == '__main__':
    main()
