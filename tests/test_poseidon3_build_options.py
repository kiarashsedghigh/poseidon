"""Build and verify all SD compile-time combinations on both C backends.

Run: python3 tests/test_poseidon3_build_options.py
Uses CC (default gcc); temporary builds leave the checkout untouched.
"""
import os
from pathlib import Path
import shlex
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parents[1]


def run(args):
    subprocess.run(args, cwd=ROOT, check=True)


def main():
    compiler = shlex.split(os.environ.get('CC', 'gcc'))
    common = ['-shared', '-fPIC', '-O3', '-Wall', '-Wextra', '-Werror', '-Isources']
    fields = ['sources/f251.c', 'sources/f251_int128.c', 'sources/poseidon_rc.c']
    with tempfile.TemporaryDirectory(prefix='poseidon3-options-') as directory:
        for iso in (0, 1):
            backend = ['-DISO_C=1'] if iso else []
            plain = str(Path(directory) / f'plain-{iso}.so')
            run(compiler + common + backend + fields +
                ['sources/poseidon3_p5_plain.c', '-o', plain])
            for rotate in (0, 1):
                for scale in (0, 1):
                    print(f'ISO_C={iso}, SD_PERMUTE={rotate}, SD_SCALE={scale}', flush=True)
                    library = str(Path(directory) / f'sd-{iso}-{rotate}-{scale}.so')
                    flags = [f'-DSD_PERMUTE={rotate}', f'-DSD_SCALE={scale}']
                    run(compiler + common + backend + flags + fields +
                        ['sources/posedion3_p5_sd.c', '-o', library])
                    run([sys.executable, 'tests/test_poseidon3_sd.py', library])
                    run([sys.executable, 'tests/test_poseidon3_pow5_standard.py', plain, library])
    print('All 8 builds and 96 tests passed.')


if __name__ == '__main__':
    main()
