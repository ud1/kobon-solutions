import json
import unittest
from decimal import Decimal
from fractions import Fraction
from pathlib import Path

from verification.quick_check import parse_lines, replay_word
from verification.rationalize import (
    CERTIFICATE_LINE_EQUATION,
    compare_implicit_realization,
    fraction_text,
    implicit_line,
    make_certificate,
    first_sources_without_parallel_pairs,
    rationalize_for_certificate,
    rationalize_lines,
    verify_certificate_payload,
)


HARD_CONFIGURATIONS = (
    "10-1cn3vuo6pqt4l",
    "10-1dunp0e9zfkyb",
    "10-1e4haqnnyjfsm",
    "10-2gqgkycvqv4rm",
    "10-2oa6vxr4rtzpg",
    "10-3hqalf90ikufz",
    "10-3sedisclbqqx0",
    "10-axal6boj31rh",
    "10-d5eiieobbm8m",
    "10-upn8v3d0qxaf",
)


class RationalizeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.repo_root = Path(__file__).resolve().parents[1]

    def test_implicit_verifier_supports_vertical_lines(self):
        structure = replay_word("0 1 0", 3)
        lines = (
            (Fraction(1), Fraction(0), Fraction(0)),
            (Fraction(0), Fraction(1), Fraction(0)),
            (Fraction(1), Fraction(1), Fraction(1)),
        )
        self.assertEqual(compare_implicit_realization(lines, structure), (True, ""))

    def test_all_current_hard_configurations_have_exact_realizations(self):
        for name in HARD_CONFIGURATIONS:
            with self.subTest(name=name):
                source = self.repo_root / "gallery" / "data" / "10" / f"{name}.json"
                with source.open(encoding="utf-8") as stream:
                    data = json.load(stream, parse_float=Decimal)
                structure = replay_word(data["gens"], 10)
                result = rationalize_lines(
                    parse_lines(data, 10), structure, max_seed_denominator=10
                )
                payload = make_certificate(
                    source, self.repo_root, result, data, 10, structure
                )
                self.assertEqual(
                    payload["line_equation"], CERTIFICATE_LINE_EQUATION
                )
                self.assertEqual(
                    payload["lines_frac"],
                    [
                        [fraction_text(value) for value in implicit_line(line)]
                        for line in result.lines
                    ],
                )
                verified = verify_certificate_payload(payload, source, self.repo_root)
                self.assertEqual(verified, tuple(map(implicit_line, result.lines)))

    def test_first_example_of_every_plain_series_has_an_exact_realization(self):
        data_root = self.repo_root / "gallery" / "data"
        sources = first_sources_without_parallel_pairs(data_root)
        for source in sources:
            with self.subTest(source=source.relative_to(data_root).as_posix()):
                n = int(source.parent.name)
                with source.open(encoding="utf-8") as stream:
                    data = json.load(stream, parse_float=Decimal)
                structure = replay_word(data["gens"], n)
                self.assertFalse(structure.parallel_pairs)
                result = rationalize_for_certificate(
                    parse_lines(data, n),
                    structure,
                    max_seed_denominator=10,
                )
                payload = make_certificate(
                    source, self.repo_root, result, data, n, structure
                )
                verified = verify_certificate_payload(payload, source, self.repo_root)
                self.assertEqual(verified, tuple(map(implicit_line, result.lines)))


if __name__ == "__main__":
    unittest.main()
