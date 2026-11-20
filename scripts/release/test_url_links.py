"""Regression checks for prose URL targets and sentence punctuation."""
import ast
from pathlib import Path
import re


def load(path):
    node = next(n for n in ast.parse(Path(path).read_text()).body
                if isinstance(n, ast.FunctionDef) and n.name == 'tex_escape')
    scope = {'re': re}
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), 'exec'), scope)
    return scope['tex_escape']


def main():
    escape = load(Path(__file__).with_name('build_manuscript.py'))
    fixtures = {
        'See https://github.com/Megachile/cynipoid-pcwde-evolution.':
            r'See \url{https://github.com/Megachile/cynipoid-pcwde-evolution}.',
        r'https\://github.com/broadinstitute/picard':
            r'\url{https://github.com/broadinstitute/picard}',
        '(http://CAZy.org),': r'(\url{http://CAZy.org}),',
        'http://cazy.org);': r'\url{http://cazy.org});',
        'https://doi.org/10.1016/0300-9629(74)90850-0':
            r'\url{https://doi.org/10.1016/0300-9629(74)90850-0}',
        '[label](https://example.org/a-b)': r'\href{https://example.org/a-b}{label}',
        '[https://example.org/a-b](https://example.org/a-b)':
            r'\href{https://example.org/a-b}{https://example.org/a-b}',
        r'http\://www\.repeatmasker.org.': r'\url{http://www.repeatmasker.org}.',
        'A & B_1; 10^-3; α⁻²; *gene*; **bold**':
            r'A \& B\_1; \(10^{-3}\); α\textsuperscript{\ensuremath{-2}}; \textit{gene}; \textbf{bold}',
    }
    for source, expected in fixtures.items():
        assert escape(source) == expected, (source, escape(source), expected)
    print(f'URL_LINK_FIXTURES_PASS cases={len(fixtures)}')


if __name__ == '__main__':
    main()
