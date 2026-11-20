"""Fail-closed input pins for the re-runnable model-less recovery consumers.

Every file a consumer reads must have a row in the sibling RECOVERY_PINS.tsv
(column path, relative to the project root, or 00_REVIEW/modelless_recovery/<name>
for retained recovery tables) and must match it before it is opened for use.
All required keys are checked for presence before any file is hashed, and a
consumer refuses any read that was not verified first.

A pin is a full SHA256, or "fp..." where the only historical record is the
provenance-ledger fingerprint of a file over 512 MiB (SHA256 of the decimal
size, the first 16 MiB and the last 16 MiB, first two hex digits replaced by
"fp"; tools/provenance/prov.py file_hash). RECOVERY_PINS.tsv names the authority
of every pin.
"""
import csv, hashlib
from pathlib import Path

PINS_FILE = Path(__file__).resolve().with_name('RECOVERY_PINS.tsv')
RECOVERY_PREFIX = '00_REVIEW/modelless_recovery/'
FP_BYTES = 16 << 20


class PinError(RuntimeError):
    pass


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for b in iter(lambda: f.read(1 << 22), b''):
            h.update(b)
    return h.hexdigest()


def fingerprint(path):
    size = Path(path).stat().st_size
    h = hashlib.sha256(str(size).encode())
    with open(path, 'rb') as f:
        h.update(f.read(FP_BYTES))
        f.seek(max(size - FP_BYTES, 0))
        h.update(f.read())
    return 'fp' + h.hexdigest()[2:]


def load(path=PINS_FILE):
    pins = {}
    with open(path) as f:
        for r in csv.DictReader(f, delimiter='\t'):
            old = pins.setdefault(r['path'], r['sha256'])
            if old != r['sha256']:
                raise PinError(f"PIN_CONFLICT {r['path']}")
    return pins


class Verifier:
    def __init__(self, pins_file=PINS_FILE):
        self.pins = load(pins_file)
        self.verified = {}

    def require(self, items):
        """items: iterable of (pin key, path). Presence first, then content."""
        items = list(items)
        missing = sorted({k for k, _ in items if k not in self.pins})
        if missing:
            raise PinError('PIN_MISSING ' + ';'.join(missing))
        for key, path in items:
            want = self.pins[key]
            got = fingerprint(path) if want.startswith('fp') else sha256(path)
            if got != want:
                raise PinError(f'PIN_MISMATCH {key} expected={want} observed={got}')
            self.verified[str(path)] = key

    def check(self, path):
        if str(path) not in self.verified:
            raise PinError(f'UNVERIFIED_READ {path}')
        return path
