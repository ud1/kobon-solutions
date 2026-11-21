import json
import unittest
from decimal import Decimal
from fractions import Fraction
from pathlib import Path
from random import Random
from tempfile import TemporaryDirectory

from verification.quick_check import (
    CheckError,
    WordStructure,
    check_one_file,
    compare_realization,
    count_triangles,
    geometry_rows,
    parse_lines,
    quick_rationalize,
    replay_word,
    validate_catalogs,
)


class QuickCheckTests(unittest.TestCase):
    def test_suboptimal_series_is_discovered_and_checked(self):
        with TemporaryDirectory() as temporary_directory:
            root = Path(temporary_directory)
            data_root = root / "data"
            series = data_root / "3-s"
            series.mkdir(parents=True)
            source = series / "example.json"
            source.write_text(
                json.dumps({"gens": "0 1 0", "lines": [[2, 0], [1, -1], [0, 0]]}),
                encoding="utf-8",
            )
            (data_root / "index.json").write_text('[{"id": "3-s"}]', encoding="utf-8")
            (data_root / "3-s.json").write_text('["3-s/example"]', encoding="utf-8")

            self.assertEqual(validate_catalogs(data_root), ([source], []))
            result = check_one_file(str(source), str(data_root), str(root / "certificates"))
            self.assertEqual((result.n, result.triangle_count, result.status), (3, 1, "exact"))

    def test_integer_comparison_agrees_with_sorted_fraction_rows(self):
        random = Random(1729)
        for n in range(3, 9):
            # A reduced word for the complete reversal, with no triples.
            gens = " ".join(str(i) for width in range(n - 1, 0, -1) for i in range(width))
            structure = replay_word(gens, n)
            for _ in range(40):
                lines = tuple(
                    (Fraction(random.randrange(-20, 21), random.randrange(1, 11)),
                     Fraction(random.randrange(-20, 21), random.randrange(1, 11)))
                    for _ in range(n)
                )
                try:
                    rows, parallel = geometry_rows(lines)
                    matches = parallel == structure.parallel_pairs and all(
                        actual == wanted or actual == tuple(reversed(wanted))
                        for wanted, actual in zip(structure.rows, rows)
                    )
                except CheckError:
                    matches = False
                else:
                    # Also exercise successful comparisons, including random
                    # parallel classes and multiple intersection points.
                    realized = WordStructure(rows, parallel, ())
                    self.assertEqual(compare_realization(lines, realized), (True, ""))
                self.assertEqual(compare_realization(lines, structure)[0], matches)

    def test_reversed_line_orientation_is_accepted(self):
        structure = replay_word("0 1 0", 3)
        lines = ((Fraction(-2), Fraction(0)),
                 (Fraction(-1), Fraction(-1)),
                 (Fraction(0), Fraction(0)))
        self.assertEqual(compare_realization(lines, structure), (True, ""))

    def test_unexpected_concurrency_and_coincident_lines_are_rejected(self):
        structure = replay_word("0 1 0", 3)
        concurrent = ((Fraction(2), Fraction(0)),
                      (Fraction(1), Fraction(0)),
                      (Fraction(0), Fraction(0)))
        self.assertFalse(compare_realization(concurrent, structure)[0])
        coincident = (concurrent[0], concurrent[0], concurrent[2])
        self.assertFalse(compare_realization(coincident, structure)[0])

    def test_exact_simple_arrangement(self):
        structure = replay_word("0 1 0", 3)
        lines = parse_lines(
            {"lines": [[Decimal("2"), Decimal("0")], [Decimal("1"), Decimal("-1")], [Decimal("0"), Decimal("0")]]},
            3,
        )
        self.assertEqual(compare_realization(lines, structure), (True, ""))
        self.assertEqual(count_triangles(structure.rows), 1)

    def test_quick_rationalization_repairs_a_triple(self):
        structure = replay_word("0*", 3)
        lines = parse_lines(
            {"lines": [[Decimal("2"), Decimal("0")], [Decimal("1"), Decimal("0.001")], [Decimal("0"), Decimal("0")]]},
            3,
        )
        self.assertFalse(compare_realization(lines, structure)[0])
        repaired = quick_rationalize(lines, structure)
        self.assertEqual(compare_realization(repaired, structure), (True, ""))
        self.assertEqual(count_triangles(structure.rows), 0)

    def test_exact_parallel_pair(self):
        structure = replay_word("1 0", 3)
        lines = parse_lines(
            {"lines": [[Decimal("1"), Decimal("0")], [Decimal("1"), Decimal("1")], [Decimal("0"), Decimal("0")]]},
            3,
        )
        self.assertEqual(structure.parallel_pairs, frozenset({(0, 1)}))
        self.assertEqual(compare_realization(lines, structure), (True, ""))

    def test_wrong_crossing_order_is_not_accepted(self):
        structure = replay_word("0 1 0 2 1 0", 4)
        lines = parse_lines(
            {
                "lines": [
                    [Decimal("3"), Decimal("0")],
                    [Decimal("2"), Decimal("1")],
                    [Decimal("1"), Decimal("6")],
                    [Decimal("0"), Decimal("10")],
                ]
            },
            4,
        )
        self.assertFalse(compare_realization(lines, structure)[0])
        repaired = quick_rationalize(lines, structure)
        self.assertFalse(compare_realization(repaired, structure)[0])


if __name__ == "__main__":
    unittest.main()
