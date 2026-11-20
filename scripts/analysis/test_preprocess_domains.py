#!/usr/bin/env python3
"""Behavioral tests for the domain-preparation command's output boundary."""
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import importlib.util
import itertools

SCRIPT = Path(__file__).with_name('preprocess_domains.py')
FAMILIES = ('GH5_2', 'PL1_8', 'PL4_4')
spec = importlib.util.spec_from_file_location('preprocess', SCRIPT)
preprocess = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preprocess)


def hit(start, end, evalue=1e-10, hmm_start=1):
    return dict(ali_from=start, ali_to=end, i_evalue=evalue,
                hmm_from=hmm_start, hmm_to=90, hmm_len=100, tlen=1000)


class OverlapTests(unittest.TestCase):
    def test_best_later_hit(self):
        a, b = hit(1, 120, 1e-6), hit(2, 120, 1e-20)
        self.assertEqual(preprocess.filter_overlapping([a, b]), [b])
        self.assertEqual(preprocess.filter_overlapping([a, b], 'legacy-diagnostic'), [a])

    def test_bridge_checks_every_retained_interval(self):
        a, b, bridge = hit(1, 20, 1e-40), hit(30, 50, 1e-30), hit(15, 35, 1e-20)
        self.assertEqual(preprocess.filter_overlapping([bridge, b, a]), [a, b])
        bridge['i_evalue'] = 1e-50
        self.assertEqual(preprocess.filter_overlapping([a, b, bridge]), [bridge])

    def test_all_permutations_chain_nested_ties(self):
        a, b, c = hit(1, 20, 1e-20), hit(20, 40, 1e-30), hit(41, 60, 1e-20)
        nested, tied, disjoint = hit(23, 30, 1e-10), hit(20, 45, 1e-30), hit(80, 90, 1e-5)
        values = [a, b, c, nested, tied, disjoint]
        original = [dict(h) for h in values]
        for order in itertools.permutations(values):
            self.assertEqual(preprocess.filter_overlapping(list(order)), [b, c, disjoint])
        self.assertEqual(values, original)

    def test_coordinate_boundaries_and_single_residue(self):
        a, shared, adjacent = hit(1, 10, 1e-30), hit(10, 10), hit(11, 20)
        self.assertEqual(preprocess.filter_overlapping([a, shared, adjacent]), [a, adjacent])
        self.assertEqual(preprocess.filter_overlapping([shared, dict(shared)]), [shared])

    def test_complete_tie_order_and_coordinate_output(self):
        a, b = hit(20, 40, 1e-20, 2), hit(20, 40, 1e-20, 1)
        early = hit(1, 10, 1e-6)
        self.assertEqual(preprocess.filter_overlapping([a, b, early]), [early, b])


class ParserTests(unittest.TestCase):
    def parse(self, line):
        with tempfile.TemporaryDirectory(prefix='pcwde_domtbl_test_') as td:
            p = Path(td) / 'hits.txt'
            p.write_text(line)
            return preprocess.parse_domtbl(p, 'GH5_2.hmm')

    def row(self):
        return 'protein - 200 GH5_2.hmm - 100 1e-20 80 0 1 1 1e-20 1e-20 80 0 1 90 11 100 10 102 0.98'.split()

    def test_hmmsearch_orientation_without_description(self):
        h = self.parse(' '.join(self.row()))['protein'][0]
        self.assertEqual((h['ali_from'], h['ali_to'], h['hmm_len'], h['tlen']), (11, 100, 100, 200))

    def test_hmmscan_rejected(self):
        with self.assertRaises(ValueError):
            self.parse('# Program: hmmscan\n')
        row = self.row()
        row[0], row[3] = row[3], row[0]
        with self.assertRaises(ValueError):
            self.parse(' '.join(row))

    def test_invalid_coordinates_evalues_and_short_rows(self):
        for field, value in ((17, '0'), (18, '201'), (15, '0'), (16, '101'), (12, 'nan'), (12, '-1')):
            row = self.row()
            row[field] = value
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                self.parse(' '.join(row))
        with self.assertRaises(ValueError):
            self.parse('truncated row\n')


class OutputBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='pcwde_preprocess_test_')
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.inputs = self.root / 'inputs'
        self.inputs.mkdir()
        for family in FAMILIES:
            (self.inputs / (family + '.faa')).write_text('>example\n' + 'A' * 200 + '\n')
            (self.inputs / (family + '_domtbl.txt')).write_text('# Program: hmmsearch\n')
        self.before = {p.name: p.read_bytes() for p in self.inputs.iterdir()}

    def run_cli(self, *args):
        env = os.environ.copy()
        env['PCWDE_PROTEIN_DIR'] = str(self.inputs)
        return subprocess.run([sys.executable, '-B', str(SCRIPT), *map(str, args)],
                              cwd=self.root, env=env, capture_output=True, text=True)

    def assert_inputs_preserved(self):
        self.assertEqual(self.before, {p.name: p.read_bytes() for p in self.inputs.iterdir()})

    def assert_outputs(self, out):
        self.assertEqual({p.name for p in out.iterdir()},
                         {'preprocess_report.tsv', *(f + '_trimmed.faa' for f in FAMILIES)})
        for family in FAMILIES:
            self.assertEqual((out / (family + '_trimmed.faa')).read_bytes(), self.before[family + '.faa'])
        self.assert_inputs_preserved()

    def test_default(self):
        result = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assert_outputs(self.root / 'build/proteins')

    def test_explicit(self):
        out = self.root / 'chosen/output'
        result = self.run_cli('--protein-dir', self.inputs, '--output-dir', out)
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assert_outputs(out)

    def test_missing_each_required_input_creates_no_output(self):
        for name, payload in self.before.items():
            with self.subTest(name=name):
                p = self.inputs / name
                p.unlink()
                result = self.run_cli()
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('required input', result.stderr)
                self.assertFalse((self.root / 'build').exists())
                p.write_bytes(payload)
        self.assert_inputs_preserved()

    def test_missing_input_root(self):
        result = self.run_cli('--protein-dir', self.root / 'missing')
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse((self.root / 'build').exists())

    def test_existing_outputs_untouched(self):
        out = self.root / 'existing'
        out.mkdir()
        for empty in (True, False):
            if not empty:
                (out / 'sentinel').write_text('preserve me')
            result = self.run_cli('--output-dir', out)
            self.assertNotEqual(result.returncode, 0)
            self.assertEqual(list(out.iterdir()), [] if empty else [out / 'sentinel'])
        self.assertEqual((out / 'sentinel').read_text(), 'preserve me')
        self.assert_inputs_preserved()

    def test_input_output_conflicts(self):
        for path in (self.inputs, self.inputs / 'nested', self.root):
            with self.subTest(path=path):
                result = self.run_cli('--output-dir', path)
                self.assertNotEqual(result.returncode, 0)
                self.assertIn('must be separate', result.stderr)
        alias = self.root / 'alias'
        alias.symlink_to(self.inputs, target_is_directory=True)
        self.assertNotEqual(self.run_cli('--output-dir', alias).returncode, 0)
        self.assert_inputs_preserved()

    def test_file_output_conflict(self):
        out = self.root / 'file'
        out.write_text('preserve file')
        self.assertNotEqual(self.run_cli('--output-dir', out).returncode, 0)
        self.assertEqual(out.read_text(), 'preserve file')
        self.assert_inputs_preserved()


if __name__ == '__main__':
    unittest.main()
