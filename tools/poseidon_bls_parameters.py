#!/usr/bin/env python3
"""Generate Poseidon round constants and sequential Cauchy matrices for BLS12-381.

The LFSR and rejection sampling follow the previously supplied constants.py.
This generates parameters only; it is not a permutation implementation.
"""
import argparse
import json

P = 0x73EDA753299D7D483339D80809A1D80553BDA402FFFE5BFEFFFFFFFF00000001
# Poseidon paper Table 2; poseidon-X is an additional requested configuration.
INSTANCES = {
    'poseidon-128-t3': (3, 8, 57),
    'poseidon-128-t5': (5, 8, 60),
    'poseidon-80-t3': (3, 8, 33),
    'poseidon-80-t5': (5, 8, 35),
    'poseidon-256-t6': (6, 8, 120),
    'poseidon-256-t10': (10, 8, 120),
    'poseidon-x': (9, 8, 57),
}


def generate_constants(t, rf, rp):
    state = [int(c) for c in (f'{1:02b}{1:04b}{255:012b}{t:012b}'
                              f'{rf:010b}{rp:010b}' + '1' * 30)]

    def advance():
        bit = state[0] ^ state[13] ^ state[23] ^ state[38] ^ state[51] ^ state[62]
        state.pop(0)
        state.append(bit)
        return bit

    for _ in range(160):
        advance()
    result = []
    while len(result) < (rf + rp) * t:
        value = 0
        count = 0
        while count < 255:
            select, bit = advance(), advance()
            if select:
                value = (value << 1) | bit
                count += 1
        if value < P:
            result.append(value)
    return [result[i:i+t] for i in range(0, len(result), t)]


def matrix(t):
    return [[pow(t + i + j, -1, P) for j in range(t)] for i in range(t)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('instance', choices=INSTANCES)
    args = parser.parse_args()
    t, rf, rp = INSTANCES[args.instance]
    data = dict(instance=args.instance, t=t, full_rounds=rf, partial_rounds=rp,
                exponent=5, modulus=hex(P), matrix=matrix(t),
                round_constants=generate_constants(t, rf, rp))
    print(json.dumps(data, indent=2))


if __name__ == '__main__':
    main()
