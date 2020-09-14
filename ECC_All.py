import functools
import multiprocessing
import sys
import time

import ECC_BruteForce
import ECC_BabyGiantStep
import ECC_PollardRho

from ECC_Main import Elliptic_curve as E_curve


def compute_one(func, x):
    p = E_curve.g
    q = E_curve.multiply(x, p)

    try:
        y, steps = func(p, q)
    except Exception as exc:
        return x, str(exc)

    return x, y, steps


def compute_all(func):
    total_steps = 0
    compute_func = functools.partial(compute_one, func)

    with multiprocessing.Pool() as pool:
        results = pool.imap_unordered(compute_func, range(E_curve.n))

        for i, (x, y, steps) in enumerate(results):
            total_steps += steps

            if x != y:
                print('\nERROR: expected {}, got: {}'.format(x, y))

            if i % 100 == 0:
                print('\rComputing all logarithms: {:.2f}% done'
                      .format(100 * i / (E_curve.n - 1)), end='')
                sys.stdout.flush()

    print('\rComputing all logarithms: 100.00% done')

    return total_steps / E_curve.n


def main():
    all_mods = [
        ECC_BruteForce,
        ECC_BabyGiantStep,
        ECC_PollardRho,
    ]

    print('Curve order: {}'.format(E_curve.n))

    for mod in all_mods:
        print('Using {}'.format(mod.__name__))

        start = time.monotonic()
        average_steps = compute_all(mod.log)
        stop = time.monotonic()

        total_seconds = stop - start
        minutes = int(total_seconds // 60)
        seconds = round(total_seconds - 60 * minutes)

        print('Took {}m {}s ({} steps on average)'
              .format(minutes, seconds, round(average_steps)))


if __name__ == '__main__':
    main()