#!/usr/bin/env python3
"""
combine modular counts with the chinese remainder theorem

reads one residue per line.  the k-th residue is taken modulo the k-th prime
of the table below, which mirrors prime[] in na-query.c / orig.c -- so residues
must be listed in prime order starting from prime[0].  blank lines and #
comments are ignored.

example (these two residues combine to 736983568, the count for the 4x4 square):

  25278      # count mod prime[0] = 29443
  10602      # count mod prime[1] = 29453
"""

import argparse
import math
import sys
from pathlib import Path

PRIMES = [
    29443, 29453, 29473, 29483, 29501, 29527, 29531, 29537, 29567, 29569,
    29573, 29581, 29587, 29599, 29611, 29629, 29633, 29641, 29663, 29669,
    29671, 29683, 29717, 29723, 29741, 29753, 29759, 29761, 29789, 29803,
    29819, 29833, 29837, 29851, 29863, 29867, 29873, 29879, 29881, 29917,
    29921, 29927, 29947, 29959, 29983, 29989, 30011, 30013, 30029, 30047,
    30059, 30071, 30089, 30091, 30097, 30103, 30109, 30113, 30119, 30133,
    30137, 30139, 30161, 30169, 30181, 30187, 30197, 30203, 30211, 30223,
    30241, 30253, 30259, 30269, 30271, 30293, 30307, 30313, 30319, 30323,
    30341, 30347, 30367, 30389, 30391, 30403, 30427, 30431, 30449, 30467,
    30469, 30491, 30493, 30497, 30509, 30517, 30529, 30539, 30553, 30557,
    30559, 30577, 30593, 30631, 30637, 30643, 30649, 30661, 30671, 30677,
    30689, 30697, 30703, 30707, 30713, 30727, 30757, 30763, 30773, 30781,
    30803, 30809, 30817, 30829, 30839, 30841, 30851, 30853, 30859, 30869,
    30871, 30881, 30893, 30911, 30931, 30937, 30941, 30949, 30971, 30977,
    30983, 31013, 31019, 31033, 31039, 31051, 31063, 31069, 31079, 31081,
    31091, 31121, 31123, 31139, 31147, 31151, 31153, 31159, 31177, 31181,
    31183, 31189, 31193, 31219, 31223, 31231, 31237, 31247, 31249, 31253,
    31259, 31267, 31271, 31277, 31307, 31319, 31321, 31327, 31333, 31337,
    31357, 31379, 31387, 31391, 31393, 31397, 31469, 31477, 31481, 31489,
    31511, 31513, 31517, 31531, 31541, 31543, 31547, 31567, 31573, 31583,
]


def parse_residues(lines):
    residues = []
    for line in lines:
        line = line.split('#', 1)[0].strip()
        if line:
            residues.append(int(line))
    return residues


def combine_pair(a, m, b, n):
    g = math.gcd(m, n)
    if (b - a) % g != 0:
        raise ValueError(f'inconsistent congruences mod {m} and mod {n}')

    m_reduced = m // g
    n_reduced = n // g
    step = ((b - a) // g * pow(m_reduced, -1, n_reduced)) % n_reduced
    modulus = m * n_reduced
    value = (a + m * step) % modulus
    return value, modulus


def combine(congruences):
    value = 0
    modulus = 1
    for residue, next_modulus in congruences:
        value, modulus = combine_pair(value, modulus, residue, next_modulus)
    return value, modulus


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('files', nargs='*', help='files to read, or stdin')
    args = parser.parse_args()

    if args.files:
        lines = []
        for filename in args.files:
            lines.extend(Path(filename).read_text().splitlines())
    else:
        lines = sys.stdin.read().splitlines()

    residues = parse_residues(lines)
    if not residues:
        raise SystemExit('no residues provided')
    if len(residues) > len(PRIMES):
        raise SystemExit(
            f'{len(residues)} residues given but only '
            f'{len(PRIMES)} primes are tabulated')

    congruences = [(r % PRIMES[i], PRIMES[i])
                   for i, r in enumerate(residues)]
    value, modulus = combine(congruences)
    print(f'value {value}')
    print(f'modulus {modulus}')


if __name__ == '__main__':
    main()
