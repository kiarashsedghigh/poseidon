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
import test_poseidon3_sd as sd

LIBRARY = sd.ROOT / 'sources/lib_pos_pow5_standard.so'
SD_LIBRARY = sd.ROOT / 'sources/lib_pos_sd.so'


class StandardContext(ctypes.Structure):
    _fields_ = [('modulus', sd.Felt), ('round_constants', sd.TABLE),
                ('mds_full', sd.Matrix), ('mds_partial', sd.Matrix)]


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
        context_pointer = ctypes.POINTER(StandardContext)
        cls.lib.poseidon3_init_standard.argtypes = [context_pointer]
        cls.lib.poseidon3_init_standard.restype = None
        cls.ctx = StandardContext()
        cls.lib.poseidon3_init_standard(ctypes.byref(cls.ctx))
        cls.lib.permutation_3_pow5_standard.argtypes = [field_pointer, context_pointer]
        cls.lib.permutation_3_pow5_standard.restype = None
        cls.lib.mix_layer_3_standard.argtypes = [field_pointer, ctypes.c_uint8, context_pointer]
        cls.lib.mix_layer_3_standard.restype = ctypes.c_int
        cls.lib.round_3_pow5_standard.argtypes = [field_pointer, ctypes.c_int, ctypes.c_uint8, context_pointer]
        cls.lib.round_3_pow5_standard.restype = ctypes.c_int
        cls.sd_lib.poseidon3_init_sd.argtypes = [ctypes.POINTER(sd.Context), field_pointer]
        cls.sd_lib.poseidon3_init_sd.restype = ctypes.c_int
        for name in ('permutation_3_sd',):
            getattr(cls.sd_lib, name).argtypes = [field_pointer, ctypes.POINTER(sd.Context)]
            getattr(cls.sd_lib, name).restype = None
        rng = random.Random(271828)
        cls.cases = [[0]*3, [3, 1, 2], [sd.P-1]*3, [sd.R-1, sd.P, 2*sd.P]]
        cls.cases += [[rng.randrange(sd.R) for _ in range(3)] for _ in range(100)]

    def test_shared_modulus(self):
        plain_prime = sd.Felt.in_dll(self.lib, 'F251_MODULUS')
        sd_prime = sd.Felt.in_dll(self.sd_lib, 'F251_MODULUS')
        ctx = sd.Context()
        self.assertEqual(self.sd_lib.poseidon3_init_sd(
            ctypes.byref(ctx), sd.state([1, 1, 2])), 0)
        self.assertEqual(list(plain_prime), list(sd.felt(sd.P)))
        self.assertEqual(list(plain_prime), list(sd_prime))
        self.assertEqual(list(plain_prime), list(ctx.modulus))

    def test_every_round_and_complete_output(self):
        for case in self.cases:
            actual = sd.state([x*sd.R % sd.P for x in case])
            expected = [x % sd.P for x in case]
            self.assertEqual(sd.values(actual, True), expected)
            for r in range(91):
                with self.subTest(case=case, round=r):
                    mode = 15 if r < 4 or r >= 87 else 1
                    self.assertEqual(self.lib.round_3_pow5_standard(actual, r, mode, ctypes.byref(self.ctx)), 0)
                    expected = reference_round_standard(expected, r)
                    self.assertEqual(sd.values(actual, True), expected)
            complete = sd.state(case)
            self.lib.permutation_3_pow5_standard(complete, ctypes.byref(self.ctx))
            self.assertEqual([sd.integer(x) for x in complete], expected)

    def test_sd_comparison(self):
        permute = bool(ctypes.c_int.in_dll(self.sd_lib, 'poseidon3_permute_sd').value)
        scale = bool(ctypes.c_int.in_dll(self.sd_lib, 'poseidon3_scale_sd').value)
        for supplied in ((1,1,1), (1,1,2), (3,5,7), (sd.P-1,sd.P-2,2), (sd.P+1,2*sd.P+3,sd.R-1)):
            ctx = sd.Context()
            self.assertEqual(self.sd_lib.poseidon3_init_sd(ctypes.byref(ctx), sd.state(supplied)), 0)
            a = supplied if scale else (1,1,1)
            order = [2,0,1] if permute else [0,1,2]
            for case in self.cases[:12]:
                with self.subTest(a=a, permute=permute, case=case):
                    baseline, encoded = sd.state(case), sd.state(case)
                    self.lib.permutation_3_pow5_standard(baseline, ctypes.byref(self.ctx))
                    self.sd_lib.permutation_3_sd(encoded, ctypes.byref(ctx))
                    # Decode only in the test; the C permutation returns SD state.
                    decoded = [sd.integer(encoded[(i+1)%3 if permute else i])*pow(a[i], -1, sd.P) % sd.P for i in range(3)]
                    self.assertEqual(decoded, [sd.integer(x) for x in baseline])
                    expected = [sd.integer(baseline[j])*a[j] % sd.P for j in order]
                    self.assertEqual([sd.integer(x) for x in encoded], expected)
                    if a == (1,1,1) and not permute:
                        self.assertEqual(bytes(encoded), bytes(baseline))


    def test_mix_layer(self):
        for case in self.cases:
            for mode in (15, 1):
                actual = sd.state(case)
                values = sd.values(actual, True)
                self.assertEqual(self.lib.mix_layer_3_standard(actual, mode, ctypes.byref(self.ctx)), 0)
                expected = [sum(m*x for m,x in zip(row,values)) % sd.P for row in sd.MATRIX]
                self.assertEqual(sd.values(actual, True), expected)

    def test_initialized_context(self):
        self.assertEqual(list(self.ctx.modulus), list(sd.felt(sd.P)))
        for r in range(91):
            self.assertEqual(sd.values(self.ctx.round_constants[r], True), list(sd.KEYS[r]))
        for matrix in (self.ctx.mds_full, self.ctx.mds_partial):
            self.assertTrue(all(sd.integer(x) < sd.P for row in matrix for x in row))
            self.assertEqual([sd.values(row, True) for row in matrix],
                             [[x % sd.P for x in row] for row in sd.MATRIX])
        before = bytes(self.ctx)
        for case in self.cases[:12]:
            actual = sd.state(case)
            self.lib.permutation_3_pow5_standard(actual, ctypes.byref(self.ctx))
        self.assertEqual(bytes(self.ctx), before)
        actual = sd.state([1, 2, 3])
        before = bytes(actual)
        self.assertEqual(self.lib.mix_layer_3_standard(actual, 0, ctypes.byref(self.ctx)), -1)
        self.assertEqual(bytes(actual), before)

    def test_invalid_round_without_mutation(self):
        for r, mode in ((-1,15), (91,1), (0,0), (0,2)):
            actual = sd.state([1,2,3])
            before = bytes(actual)
            self.assertEqual(self.lib.round_3_pow5_standard(actual,r,mode,ctypes.byref(self.ctx)), -1)
            self.assertEqual(bytes(actual),before)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        LIBRARY = Path(sys.argv[1]).resolve()
    if len(sys.argv) > 2:
        SD_LIBRARY = Path(sys.argv[2]).resolve()
    unittest.main(argv=[sys.argv[0]], verbosity=2)
