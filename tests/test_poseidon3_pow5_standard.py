"""Non-SD baseline: Python parity and comparison with SD after decoding.

Run from the repository root after building:
  python3 tests/test_poseidon3_pow5_standard.py
  python3 tests/test_poseidon3_pow5_standard.py sources/lib_pos_pow5_standard_asm.so sources/lib_pos_sd_asm.so
"""
import ctypes
from pathlib import Path
import random
import sys
import unittest
import contextlib
import io
import test_poseidon3_sd as sd

LIBRARY = sd.ROOT / 'sources/lib_pos_pow5_standard.so'
SD_LIBRARY = sd.ROOT / 'sources/lib_pos_sd.so'


def reference_round_standard(state, r):
    state = [(x+c) % sd.P for x, c in zip(state, sd.KEYS[r])]
    for i in range(3) if r < 4 or r >= 87 else (0,):
        state[i] = pow(state[i], 5, sd.P)
    return [sum(m*x for m, x in zip(row, state)) % sd.P for row in sd.MATRIX]


class FifthPowerStandardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib = ctypes.CDLL(str(LIBRARY))
        cls.sd_lib = ctypes.CDLL(str(SD_LIBRARY))
        field_pointer = ctypes.POINTER(sd.Felt)
        for name in ('state_to_montgomery_3_standard', 'state_from_montgomery_3_standard',
                     'mix_layer_3_standard', 'permutation_3_pow5_standard'):
            getattr(cls.lib, name).argtypes = [field_pointer]
            getattr(cls.lib, name).restype = None
        cls.lib.round_3_pow5_standard.argtypes = [field_pointer, ctypes.c_int, ctypes.c_uint8]
        cls.lib.round_3_pow5_standard.restype = ctypes.c_int
        cls.lib.poseidon_hash_3_pow5_standard.argtypes = [ctypes.POINTER(ctypes.c_uint64), field_pointer]
        cls.lib.poseidon_hash_3_pow5_standard.restype = None
        cls.sd_lib.poseidon3_init_sd.argtypes = [ctypes.POINTER(sd.Context), field_pointer, ctypes.c_bool]
        cls.sd_lib.poseidon3_init_sd.restype = ctypes.c_int
        for name in ('permutation_3_sd', 'permutation_3_restored_sd'):
            getattr(cls.sd_lib, name).argtypes = [field_pointer, ctypes.POINTER(sd.Context)]
            getattr(cls.sd_lib, name).restype = None
        rng = random.Random(271828)
        cls.cases = [[0]*3, [3, 1, 2], [sd.P-1]*3, [sd.R-1, sd.P, 2*sd.P]]
        cls.cases += [[rng.randrange(sd.R) for _ in range(3)] for _ in range(100)]

    def test_every_round_and_complete_output(self):
        for case in self.cases:
            actual = sd.state(case)
            self.lib.state_to_montgomery_3_standard(actual)
            expected = [x % sd.P for x in case]
            self.assertEqual(sd.values(actual, True), expected)
            for r in range(91):
                with self.subTest(case=case, round=r):
                    mode = 15 if r < 4 or r >= 87 else 1
                    self.assertEqual(self.lib.round_3_pow5_standard(actual, r, mode), 0)
                    expected = reference_round_standard(expected, r)
                    self.assertEqual(sd.values(actual, True), expected)
            self.lib.state_from_montgomery_3_standard(actual)
            self.assertEqual([sd.integer(x) for x in actual], expected)
            complete = sd.state(case)
            self.lib.permutation_3_pow5_standard(complete)
            self.assertEqual(bytes(complete), bytes(actual))

    def test_sd_comparison(self):
        for a in ((1,1,1), (1,1,2), (3,5,7), (sd.P-1,sd.P-2,2), (sd.P+1,2*sd.P+3,sd.R-1)):
            for permute in (False, True):
                ctx = sd.Context()
                self.assertEqual(self.sd_lib.poseidon3_init_sd(ctypes.byref(ctx), sd.state(a), permute), 0)
                order = [2,0,1] if permute else [0,1,2]
                for case in self.cases[:12]:
                    with self.subTest(a=a, permute=permute, case=case):
                        baseline, encoded, restored = sd.state(case), sd.state(case), sd.state(case)
                        self.lib.permutation_3_pow5_standard(baseline)
                        self.sd_lib.permutation_3_sd(encoded, ctypes.byref(ctx))
                        self.sd_lib.permutation_3_restored_sd(restored, ctypes.byref(ctx))
                        self.assertEqual(bytes(restored), bytes(baseline))
                        expected = [sd.integer(baseline[j])*a[j] % sd.P for j in order]
                        self.assertEqual([sd.integer(x) for x in encoded], expected)
                        if a == (1,1,1) and not permute:
                            self.assertEqual(bytes(encoded), bytes(baseline))

    def test_hash_against_uploaded_python_and_aliasing(self):
        _, reference = sd.reference((1,1,1), False)
        for preimage in ([1,2], [0,0], [sd.P-1,sd.P+1], [sd.R-1,sd.R-1]):
            inputs = (sd.Felt*2)(*(sd.felt(x) for x in preimage))
            result = sd.Felt()
            self.lib.poseidon_hash_3_pow5_standard(result, inputs)
            with contextlib.redirect_stdout(io.StringIO()):
                expected = reference.poseidon_hash(preimage)
            self.assertEqual(sd.integer(result), expected)
            self.lib.poseidon_hash_3_pow5_standard(inputs[0], inputs)
            self.assertEqual(sd.integer(inputs[0]), expected)

    def test_mix_layer(self):
        for case in self.cases:
            actual = sd.state(case)
            self.lib.state_to_montgomery_3_standard(actual)
            self.lib.mix_layer_3_standard(actual)
            expected = [sum(m*x for m,x in zip(row,case)) % sd.P for row in sd.MATRIX]
            self.assertEqual(sd.values(actual, True), expected)

    def test_invalid_round_without_mutation(self):
        for r, mode in ((-1,15), (91,1), (0,0), (0,2)):
            actual = sd.state([1,2,3])
            before = bytes(actual)
            self.assertEqual(self.lib.round_3_pow5_standard(actual,r,mode), -1)
            self.assertEqual(bytes(actual),before)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        LIBRARY = Path(sys.argv[1]).resolve()
    if len(sys.argv) > 2:
        SD_LIBRARY = Path(sys.argv[2]).resolve()
    unittest.main(argv=[sys.argv[0]], verbosity=2)
