import math
import random

from ECC_Main import Elliptic_curve as E_curve


def log(p, q):
    assert E_curve.lying_on_curve(p)
    assert E_curve.lying_on_curve(q)

    m = int(math.sqrt(E_curve.n)) + 1

    # Compute the baby steps and store them in the 'precomputed' hash table.
    r = None
    stored_hash = {None: 0}

    for a in range(1, m):
        r = E_curve.add(r, p)
        stored_hash[r] = a

    # Now compute the giant steps and check the hash table for any
    # matching point.
    r = q
    s = E_curve.multiply(m, E_curve.negative_point(p))

    for b in range(m):
        try:
            a = stored_hash[r]
        except KeyError:
            pass
        else:
            steps = m + b
            logarithm = a + m * b
            return logarithm, steps

        r = E_curve.add(r, s)

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
