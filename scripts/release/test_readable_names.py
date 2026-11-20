"""Negative fixture for the deposit filename gate."""
import tempfile
from pathlib import Path
from check_readable_names import check
with tempfile.TemporaryDirectory() as d:
    root=Path(d); bad=root/('a'*64+'.faa');bad.write_text('>x\nA\n')
    try: check(Path(__file__).resolve().parents[2],root)
    except ValueError as e: assert 'HASH_NAME_FAIL' in str(e)
    else: raise AssertionError('hash filename accepted')
print('READABLE_NAMES_NEGATIVE_PASS')
