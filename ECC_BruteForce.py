import random

from ECC_Main import Elliptic_curve as E_curve


def log(p, q):
    assert E_curve.lying_on_curve(p)
    assert E_curve.lying_on_curve(q)

    start = random.randrange(E_curve.n)
    r = E_curve.multiply(start, p)

    for x in range(E_curve.n):
        if q == r:
            logarithm = (start + x) % E_curve.n
            steps = x + 1
            return logarithm, steps
        r = E_curve.add(r, p)

    raise AssertionError('logarithm not found')


def main():
    x = random.randrange(E_curve.n)
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