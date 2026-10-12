"""Functional SD parity with uploaded Python, using the requested C parameters.

Run: python3 tests/test_poseidon3_sd.py [sources/lib_pos_sd_asm.so]
The reference fixtures are the user's original uploads. Only field/round-count/
base-matrix/constant inputs are substituted, not their SD transformations.
"""
import ast
import contextlib
import ctypes
import io
from pathlib import Path
import random
import sys
import types
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
LIBRARY = ROOT / 'sources/lib_pos_sd.so'
P = 2**251 + 17*2**192 + 1
R = 2**256
Felt = ctypes.c_uint64 * 4
State = Felt * 3
Matrix = State * 3
TABLE = State * 91
PARAMS = {}
for node in ast.parse((ROOT / 'sage/poseidon3_data.sage').read_text()).body:
    if isinstance(node, ast.Assign):
        PARAMS[node.targets[0].id] = ast.literal_eval(node.value)
KEYS = PARAMS['RoundKeys3']
MATRIX = PARAMS['Matrix3']


class Context(ctypes.Structure):
    _fields_ = [('multipliers', State), ('round_constants', TABLE),
                ('mds_full', Matrix), ('mds_partial', Matrix),
                ('permute', ctypes.c_bool), ('exponent', ctypes.c_uint8),
                ('partial_index', ctypes.c_uint8)]


def felt(x):
    return Felt(*((x >> (64*i)) & (2**64-1) for i in range(4)))


def state(xs):
    return State(*(felt(x) for x in xs))


def integer(x):
    return sum(int(x[i]) << (64*i) for i in range(4))


def values(xs, mont=False):
    factor = pow(R, -1, P) if mont else 1
    return [integer(x)*factor % P for x in xs]


def reference(a, permute):
    """Load uploads with C parameters and preserve the actual SD Python code."""
    params = types.ModuleType('params')
    basic = types.ModuleType('basic')
    constants = types.ModuleType('constants')
    poseidon = types.ModuleType('poseidon')
    directory = ROOT / 'tests/reference_sd'
    exec(compile((directory/'params.py').read_text(), 'params.py', 'exec'), params.__dict__)
    params.A0_CONST, params.A1_CONST, params.A2_CONST = a
    params.PERMUTE = permute
    with patch.dict(sys.modules, params=params, basic=basic, constants=constants):
        exec(compile((directory/'basic.py').read_text(), 'basic.py', 'exec'), basic.__dict__)
        basic.P = P
        basic.ROUNDPARTIAL = {3: 83}
        ops = basic.PrimeFieldOps
        full, partial = ops.get_mds_matrix, ops.get_mds_matrix_partial
        inverse = ops.get_inverse

        def c_base_matrix(builder):
            # In the upload the first nine inverse calls construct the base
            # Cauchy matrix. Substitute the C base entries at those calls;
            # subsequent inverse calls (the SD formulas) remain untouched.
            entries = iter([x % P for row in MATRIX for x in row])
            remaining = [9]
            def substitute(x):
                if remaining[0]:
                    remaining[0] -= 1
                    return next(entries)
                return inverse(x)
            with patch.object(ops, 'get_inverse', substitute):
                return builder(3)

        full_matrix = c_base_matrix(full)
        partial_matrix = c_base_matrix(partial)
        # Cache the public matrices; this changes no arithmetic.
        ops.get_mds_matrix = lambda t: full_matrix
        ops.get_mds_matrix_partial = lambda t: partial_matrix
        def generate(t, rf, rp):
            assert (t, rf, rp) == (3, 8, 83)
            return [x for row in KEYS for x in row]
        constants.generate_constants = generate
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile((directory/'poseidon.py').read_text(), 'poseidon.py', 'exec'), poseidon.__dict__)
    return basic, poseidon


class SDTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib = ctypes.CDLL(str(LIBRARY))
        cls.original = ctypes.CDLL(str(ROOT/'sources/lib_pos.so'))
        cls.original.permutation_3.argtypes = [ctypes.POINTER(Felt)]
        cls.original.permutation_3.restype = None
        pointer = ctypes.POINTER(Context)
        signatures = {
            'poseidon3_init_sd': ([pointer, ctypes.POINTER(Felt), ctypes.c_bool], ctypes.c_int),
            'poseidon3_init_custom_sd': ([pointer, ctypes.POINTER(Felt), ctypes.c_bool, ctypes.c_uint8, ctypes.c_uint8], ctypes.c_int),
            'prepare_state_3_sd': ([ctypes.POINTER(Felt), pointer], None),
            'restore_state_3_sd': ([ctypes.POINTER(Felt), pointer], None),
            'permutation_3_sd': ([ctypes.POINTER(Felt), pointer], None),
            'permutation_3_restored_sd': ([ctypes.POINTER(Felt), pointer], None),
            'round_3_sd': ([ctypes.POINTER(Felt), ctypes.c_int, ctypes.c_uint8, pointer], ctypes.c_int),
            'mix_layer_3_sd': ([ctypes.POINTER(Felt), ctypes.c_uint8, pointer], ctypes.c_int),
            'poseidon_hash_3_sd': ([ctypes.POINTER(ctypes.c_uint64), ctypes.POINTER(Felt), pointer], None),
        }
        for name, (args, result) in signatures.items():
            getattr(cls.lib, name).argtypes = args
            getattr(cls.lib, name).restype = result
        cls.configs = []
        for a in ((1, 1, 2), (1, 1, 1), (3, 5, 7), (P-1, P-2, 2), (P+1, 2*P+3, R-1)):
            for permute in (False, True):
                ctx = Context()
                assert cls.lib.poseidon3_init_sd(ctypes.byref(ctx), state(a), permute) == 0
                basic, poseidon = reference(a, permute)
                cls.configs.append((a, permute, ctx, basic, poseidon))
        rng = random.Random(12345)
        cls.cases = [[0, 0, 0], [3, 1, 2], [P-1]*3, [R-1, P, 2*P]]
        cls.cases += [[rng.randrange(P) for _ in range(3)] for _ in range(8)]

    def test_constants_and_matrices(self):
        for a, permute, ctx, basic, _ in self.configs:
            with self.subTest(a=a, permute=permute):
                self.assertEqual(ctx.exponent, 5)
                self.assertEqual(ctx.partial_index, 1 if permute else 0)
                order = [2, 0, 1] if permute else [0, 1, 2]
                for r in range(91):
                    expected = [KEYS[r][j]*a[j] % P for j in order]
                    self.assertEqual(values(ctx.round_constants[r], True), expected)
                for mode, actual in (('full', ctx.mds_full), ('partial', ctx.mds_partial)):
                    # The original mixing method also performs pi*M*pi^-1.
                    for j in range(3):
                        basis = [int(i == j) for i in range(3)]
                        column = basic.PrimeFieldOps.mds_mixing(basis, mode)
                        self.assertEqual([values(row, True)[j] for row in actual], column)

    def test_every_round_and_final_outputs(self):
        for a, permute, ctx, basic, _ in self.configs:
            ops = basic.PrimeFieldOps
            order = [2, 0, 1] if permute else [0, 1, 2]
            for case in self.cases:
                actual = state(case)
                self.lib.prepare_state_3_sd(actual, ctypes.byref(ctx))
                expected = [case[j]*a[j] % P for j in order]
                self.assertEqual(values(actual, True), expected)
                for r in range(91):
                    full = r < 4 or r >= 87
                    c = [KEYS[r][j]*a[j] % P for j in order]
                    expected = ops.add_round_constants(expected, c)
                    if full:
                        expected = ops.s_boxes(expected)
                    else:
                        active = 1 if permute else 0
                        expected[active] = ops.s_box(expected[active])
                    expected = ops.mds_mixing(expected, 'full' if full else 'partial')
                    with self.subTest(a=a, permute=permute, case=case, round=r):
                        self.assertEqual(self.lib.round_3_sd(actual, r, 15 if full else 1, ctypes.byref(ctx)), 0)
                        self.assertEqual(values(actual, True), expected)
                complete = state(case)
                self.lib.permutation_3_sd(complete, ctypes.byref(ctx))
                self.assertEqual([integer(x) for x in complete], expected)
                restored = state(case)
                self.lib.permutation_3_restored_sd(restored, ctypes.byref(ctx))
                decoded = [expected[(i+1)%3 if permute else i]*pow(a[i], -1, P) % P for i in range(3)]
                self.assertEqual([integer(x) for x in restored], decoded)

    def test_uploaded_hash_function(self):
        for a, permute, ctx, _, python in self.configs:
            for preimage in ([1, 2], [0, 0], [P-1, P+1]):
                with self.subTest(a=a, permute=permute, preimage=preimage):
                    inputs = (Felt*2)(*(felt(x) for x in preimage))
                    result = Felt()
                    self.lib.poseidon_hash_3_sd(result, inputs, ctypes.byref(ctx))
                    with contextlib.redirect_stdout(io.StringIO()):
                        expected = python.poseidon_hash(preimage)
                    self.assertEqual(integer(result), expected)

    def test_cube_mode_matches_existing_c_after_restoration(self):
        for a, permute, _, _, _ in self.configs:
            ctx = Context()
            self.assertEqual(self.lib.poseidon3_init_custom_sd(ctypes.byref(ctx), state(a), permute, 3, 2), 0)
            for case in self.cases:
                original, actual = state(case), state(case)
                self.original.permutation_3(original)
                self.lib.permutation_3_restored_sd(actual, ctypes.byref(ctx))
                self.assertEqual(bytes(actual), bytes(original))

    def test_reject_invalid_configuration_without_mutation(self):
        for zero in (0, P, 2*P):
            for slot in range(3):
                a = [1, 2, 3]
                a[slot] = zero
                ctx = Context()
                ctypes.memset(ctypes.byref(ctx), 0xA5, ctypes.sizeof(ctx))
                before = bytes(ctx)
                self.assertEqual(self.lib.poseidon3_init_sd(ctypes.byref(ctx), state(a), True), -1)
                self.assertEqual(bytes(ctx), before)
        for exponent, active in ((0, 0), (4, 0), (5, 3), (3, 255)):
            ctx = Context()
            before = bytes(ctx)
            self.assertEqual(self.lib.poseidon3_init_custom_sd(ctypes.byref(ctx), state([1, 1, 2]), True, exponent, active), -1)
            self.assertEqual(bytes(ctx), before)
        ctx = self.configs[0][2]
        for r, mode in ((-1, 15), (91, 1), (0, 0), (0, 2)):
            actual = state([1, 2, 3])
            before = bytes(actual)
            self.assertEqual(self.lib.round_3_sd(actual, r, mode, ctypes.byref(ctx)), -1)
            self.assertEqual(bytes(actual), before)
        actual = state([1, 2, 3])
        before = bytes(actual)
        self.assertEqual(self.lib.mix_layer_3_sd(actual, 2, ctypes.byref(ctx)), -1)
        self.assertEqual(bytes(actual), before)

    def test_prepare_restore_and_context_reuse(self):
        for _, _, ctx, _, _ in self.configs:
            before = bytes(ctx)
            for case in self.cases:
                actual = state(case)
                self.lib.prepare_state_3_sd(actual, ctypes.byref(ctx))
                self.lib.restore_state_3_sd(actual, ctypes.byref(ctx))
                self.assertEqual([integer(x) for x in actual], [x % P for x in case])
                self.lib.permutation_3_sd(actual, ctypes.byref(ctx))
            self.assertEqual(bytes(ctx), before)


if __name__ == '__main__':
    if len(sys.argv) > 1:
        LIBRARY = Path(sys.argv[1]).resolve()
    unittest.main(argv=[sys.argv[0]], verbosity=2)
