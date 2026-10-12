"""Run after building: python3 tests/test_f251_pow5.py [path/to/library.so]."""

import ctypes
from pathlib import Path
import random
import sys
import unittest


LIBRARY = Path(__file__).resolve().parents[1] / "sources/lib_pos.so"
P = 2**251 + 17 * 2**192 + 1
R = 2**256
Felt = ctypes.c_uint64 * 4


def encode(value):
    return Felt(*((value >> (64 * i)) & (2**64 - 1) for i in range(4)))


def decode(value):
    return sum(int(value[i]) << (64 * i) for i in range(4))


class FifthPowerTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.lib = ctypes.CDLL(str(LIBRARY))
        for name in ("f251_to_montgomery", "f251_from_montgomery",
                     "f251_montgomery_pow5"):
            fn = getattr(cls.lib, name)
            fn.argtypes = [ctypes.POINTER(ctypes.c_uint64)] * 2
            fn.restype = None
        rng = random.Random(5)
        cls.cases = [0, 1, 2, P - 1, P, P + 1, 2 * P, 16 * P, R - 1]
        for bit in (64, 128, 192, 255):
            cls.cases.extend((2**bit - 1, 2**bit, 2**bit + 1))
        cls.cases += [rng.randrange(R) for _ in range(1000)]

    def test_fifth_power_and_in_place_operation(self):
        for value in self.cases:
            with self.subTest(value=value):
                source, mx, result, ordinary = encode(value), Felt(), Felt(), Felt()
                self.lib.f251_to_montgomery(mx, source)
                saved = bytes(mx)
                self.lib.f251_montgomery_pow5(result, mx)
                self.assertEqual(bytes(mx), saved)
                expected = pow(value, 5, P)
                self.assertEqual(decode(result) % P, expected * R % P)
                self.lib.f251_from_montgomery(ordinary, result)
                self.assertEqual(decode(ordinary), expected)
                self.lib.f251_montgomery_pow5(mx, mx)
                self.lib.f251_from_montgomery(mx, mx)
                self.assertEqual(decode(mx), expected)

    def test_unreduced_montgomery_inputs(self):
        # Round additions may leave any 256-bit representative, not just < p.
        inverse_r = pow(R, -1, P)
        for value in self.cases:
            with self.subTest(value=value):
                mx, result = encode(value), Felt()
                expected = pow(value * inverse_r % P, 5, P)
                self.lib.f251_montgomery_pow5(result, mx)
                self.assertEqual(decode(result) % P, expected * R % P)
                self.lib.f251_montgomery_pow5(mx, mx)
                self.lib.f251_from_montgomery(mx, mx)
                self.assertEqual(decode(mx), expected)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        LIBRARY = Path(sys.argv[1]).resolve()
    unittest.main(argv=[sys.argv[0]], verbosity=2)
