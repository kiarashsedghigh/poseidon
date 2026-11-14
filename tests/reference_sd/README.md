These two Python files are the user's `basic.py` and `params.py` reference
uploads, with only line endings and trailing whitespace normalized.
`test_poseidon3_sd.py` loads their arithmetic and SD matrix transformations
and substitutes the C prime and base matrix. Tests use the C original
`RoundKeys3` and 4/83/4 schedule, as requested.

The first nine inverse calls constructing the Python Cauchy matrix are replaced
with the C base matrix entries; inverses in the SD formulas remain unchanged.
Public transformed matrices are cached for speed. Tests operate on all three
state values and never invoke a hashing wrapper. The unused hashing-only
`poseidon.py` fixture has been removed.

See `../../docs/poseidon3_sd.md` for the actual C parameters and API.
