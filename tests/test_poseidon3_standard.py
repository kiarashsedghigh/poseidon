"""Check original constants, round states, and compatibility without SageMath.

Build in sources first, then run from the repository root:
    python3 tests/test_poseidon3_standard.py
    python3 tests/test_poseidon3_standard.py sources/lib_pos_asm.so
Only load the assembly library on an x86_64 CPU with ADX support.
"""

import ast
import ctypes
from pathlib import Path
import random
import sys
import unittest


ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / "sources/lib_pos.so"
P = 2**251 + 17 * 2**192 + 1
R = 2**256
Felt = ctypes.c_uint64 * 4
State = Felt * 3
PARAMS = {}
for node in ast.parse((ROOT / "sage/poseidon3_data.sage").read_text()).body:
    if isinstance(node, ast.Assign):
        PARAMS[node.targets[0].id] = ast.literal_eval(node.value)
KEYS = PARAMS["RoundKeys3"]
MATRIX = PARAMS["Matrix3"]


def encode(values, montgomery=False):
    values = [(v * R if montgomery else v) % P for v in values]
    return State(*(Felt(*((v >> (64 * i)) & (2**64 - 1)
                         for i in range(4))) for v in values))


def decode(state, montgomery=False):
    factor = pow(R, -1, P) if montgomery else 1
    return [(sum(int(v[i]) << (64 * i) for i in range(4)) * factor) % P
            for v in state]


def reference_round(state, round_idx):
    state = [(x + c) % P for x, c in zip(state, KEYS[round_idx])]
    full = round_idx < 4 or round_idx >= 87
    for i in range(3) if full else (2,):
        state[i] = pow(state[i], 3, P)
    return [sum(m * x for m, x in zip(row, state)) % P for row in MATRIX]


class StandardPermutationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib = ctypes.CDLL(str(LIBRARY))
        for name in ("permutation_3", "permutation_3_standard"):
            fn = getattr(cls.lib, name)
            fn.argtypes = [ctypes.POINTER(Felt)]
            fn.restype = None
        cls.lib.round_3_standard.argtypes = [ctypes.POINTER(Felt), ctypes.c_int,
                                            ctypes.c_uint8]
        cls.lib.round_3_standard.restype = None
        rng = random.Random(314159)
        cls.cases = [[0, 0, 0], [1, 1, 1], [P - 1] * 3, [0, 1, P - 1]]
        cls.cases += [[rng.randrange(P) for _ in range(3)] for _ in range(100)]

    def test_original_montgomery_constants(self):
        self.assertEqual(len(KEYS), 91)
        constants = (State * 91).in_dll(self.lib, "ROUND_CONST_MONT_P3")
        for r, key in enumerate(KEYS):
            self.assertEqual(len(key), 3)
            self.assertEqual(decode(constants[r], montgomery=True), key)

    def test_every_round_against_reference(self):
        for case in self.cases:
            state = encode(case, montgomery=True)
            expected = case[:]
            for r in range(91):
                with self.subTest(input=case, round=r):
                    mode = 0xF if r < 4 or r >= 87 else 0x1
                    self.lib.round_3_standard(state, r, mode)
                    expected = reference_round(expected, r)
                    self.assertEqual(decode(state, montgomery=True), expected)

    def test_complete_permutation(self):
        for case in self.cases:
            with self.subTest(input=case):
                original, standard = encode(case), encode(case)
                self.lib.permutation_3(original)
                self.lib.permutation_3_standard(standard)
                # Compare exact output limbs, including canonical reduction.
                self.assertEqual(bytes(standard), bytes(original))
                expected = case[:]
                for r in range(91):
                    expected = reference_round(expected, r)
                self.assertEqual(decode(standard), expected)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        LIBRARY = Path(sys.argv[1]).resolve()
    unittest.main(argv=[sys.argv[0]], verbosity=2)
