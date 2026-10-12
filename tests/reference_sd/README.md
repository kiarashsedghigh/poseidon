These three Python files are the user's SD reference uploads, with only line
endings and trailing whitespace normalized.
`test_poseidon3_sd.py` loads them into isolated modules and substitutes the C
field, 83 partial rounds, original `RoundKeys3`, and base `Matrix3`, as requested.
The original Python scaling, matrix adjustments, permutation, fifth-power S-box,
and hash input/output conventions are exercised directly.

The nine inverse calls constructing the Python Cauchy base matrix are replaced
with C's nine base matrix entries; the inverse calls used for the SD adjustment
formulas remain unchanged. Public transformed matrices are cached for speed.
The constants generator is replaced with the repository's original constants;
the Filecoin generator is intentionally not used for this C-parameter variant.

The fixtures' BLS12-381 parameters do not describe the new C path. See
`../../docs/poseidon3_sd.md` for its actual parameters and API.
