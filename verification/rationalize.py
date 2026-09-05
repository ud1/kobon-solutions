#!/usr/bin/env python3
"""Construct exact rational realizations for the hard gallery configurations.

The construction works in the dual projective plane.  A source line
``y = m*x + b`` becomes the dual point ``(m:b:1)`` and every declared triple
intersection becomes a line through three dual points.  Starting with a small
set of rational seed points, the remaining points are constructed using only
cross products.  If the incidence graph leaves one constraint, the script
solves it exactly when it is linear in one seed coordinate.  Certificates
store the result as implicit equations ``a*x + b*y = c``.

No third-party Python packages are required.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from decimal import Decimal
from fractions import Fraction
from itertools import combinations
from pathlib import Path
from typing import Sequence

if __package__ in {None, ""}:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from verification.quick_check import (
    CheckError,
    EventRows,
    Line,
    SERIES_RE,
    WordStructure,
    certificate_path,
    check_one_file,
    compare_realization,
    count_triangles,
    first_row_difference,
    parse_lines,
    quick_rationalize,
    replay_word,
    validate_catalogs,
)


class RationalizationError(CheckError):
    """The supported exact construction could not realize a configuration."""


class CertificateError(CheckError):
    """An exact realization certificate is invalid."""


class Polynomial:
    """Small univariate polynomial over the rationals."""

    def __init__(self, coefficients: object = 0):
        if isinstance(coefficients, Polynomial):
            values = coefficients.coefficients
        elif isinstance(coefficients, (int, Fraction)):
            values = (Fraction(coefficients),)
        else:
            try:
                values = tuple(Fraction(value) for value in coefficients)  # type: ignore[union-attr]
            except TypeError as exc:
                raise TypeError(f"invalid polynomial coefficients: {coefficients!r}") from exc
        self.coefficients = list(values)
        while len(self.coefficients) > 1 and not self.coefficients[-1]:
            self.coefficients.pop()

    def __add__(self, other: object) -> Polynomial:
        right = Polynomial(other)
        return Polynomial(
            (self.coefficients[index] if index < len(self.coefficients) else 0)
            + (right.coefficients[index] if index < len(right.coefficients) else 0)
            for index in range(max(len(self.coefficients), len(right.coefficients)))
        )

    def __neg__(self) -> Polynomial:
        return Polynomial(-value for value in self.coefficients)

    def __sub__(self, other: object) -> Polynomial:
        return self + -Polynomial(other)

    def __mul__(self, other: object) -> Polynomial:
        right = Polynomial(other)
        values = [Fraction(0)] * (
            len(self.coefficients) + len(right.coefficients) - 1
        )
        for left_index, left in enumerate(self.coefficients):
            for right_index, value in enumerate(right.coefficients):
                values[left_index + right_index] += left * value
        return Polynomial(values)

    @property
    def degree(self) -> int:
        return len(self.coefficients) - 1

    @property
    def is_zero(self) -> bool:
        return not any(self.coefficients)


Scalar = Fraction | Polynomial
Homogeneous = tuple[Scalar, Scalar, Scalar]
FractionPoint = tuple[Fraction, Fraction, Fraction]
ImplicitLine = tuple[Fraction, Fraction, Fraction]  # a*x + b*y = c
CERTIFICATE_LINE_EQUATION = "a*x + b*y = c"


@dataclass(frozen=True)
class ConstructionStep:
    kind: str  # "incidence-line" or "point"
    target: int
    left: int
    right: int


@dataclass(frozen=True)
class ConstructionPlan:
    seeds: tuple[int, ...]
    steps: tuple[ConstructionStep, ...]


@dataclass(frozen=True)
class RationalizationResult:
    lines: tuple[Line, ...]
    max_seed_denominator: int
    seeds: tuple[int, ...]
    adjusted_line: int | None
    adjusted_coordinate: str | None
    method: str = "dual-incidence-construction"


def cross(left: Homogeneous, right: Homogeneous) -> Homogeneous:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def dot(left: Homogeneous, right: Homogeneous) -> Scalar:
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2]


def construction_plans(
    n: int, triple_points: Sequence[frozenset[int]]
) -> tuple[ConstructionPlan, ...]:
    """Return every plan with the smallest usable set of seed dual points."""

    for seed_count in range(1, n + 1):
        plans: list[ConstructionPlan] = []
        for seeds in combinations(range(n), seed_count):
            seed_set = set(seeds)
            # A triple made entirely of arbitrary seeds would add an
            # unsupported constraint before construction has started.
            if any(point <= seed_set for point in triple_points):
                continue

            known_points = set(seeds)
            known_incidence_lines: set[int] = set()
            steps: list[ConstructionStep] = []
            changed = True
            while changed:
                changed = False
                for index, point in enumerate(triple_points):
                    if index in known_incidence_lines:
                        continue
                    known = sorted(point & known_points)
                    if len(known) >= 2:
                        known_incidence_lines.add(index)
                        steps.append(
                            ConstructionStep(
                                "incidence-line", index, known[0], known[1]
                            )
                        )
                        changed = True

                for point in range(n):
                    if point in known_points:
                        continue
                    incident = [
                        index
                        for index, lines in enumerate(triple_points)
                        if point in lines and index in known_incidence_lines
                    ]
                    if len(incident) >= 2:
                        known_points.add(point)
                        steps.append(
                            ConstructionStep("point", point, incident[0], incident[1])
                        )
                        changed = True

            if len(known_points) == n:
                plans.append(ConstructionPlan(seeds, tuple(steps)))

        if plans:
            return tuple(plans)
    return ()


def build_points(
    seed_points: dict[int, Homogeneous], plan: ConstructionPlan
) -> dict[int, Homogeneous]:
    points = dict(seed_points)
    incidence_lines: dict[int, Homogeneous] = {}
    for step in plan.steps:
        if step.kind == "incidence-line":
            incidence_lines[step.target] = cross(
                points[step.left], points[step.right]
            )
        else:
            points[step.target] = cross(
                incidence_lines[step.left], incidence_lines[step.right]
            )
    return points


def incidence_residuals(
    points: dict[int, Homogeneous],
    triple_points: Sequence[frozenset[int]],
) -> tuple[Scalar, ...]:
    result: list[Scalar] = []
    for point in triple_points:
        first, second, third = sorted(point)
        result.append(dot(points[third], cross(points[first], points[second])))
    return tuple(result)


def fraction_seed_points(
    source_lines: Sequence[Line], plan: ConstructionPlan, max_denominator: int
) -> dict[int, Homogeneous]:
    return {
        line: (
            source_lines[line][0].limit_denominator(max_denominator),
            source_lines[line][1].limit_denominator(max_denominator),
            Fraction(1),
        )
        for line in plan.seeds
    }


def candidate_seed_points(
    source_lines: Sequence[Line],
    triple_points: Sequence[frozenset[int]],
    plan: ConstructionPlan,
    max_denominator: int,
) -> list[tuple[dict[int, Homogeneous], int | None, str | None]]:
    base = fraction_seed_points(source_lines, plan, max_denominator)
    points = build_points(base, plan)
    if not any(incidence_residuals(points, triple_points)):
        return [(base, None, None)]

    candidates: list[tuple[dict[int, Homogeneous], int | None, str | None]] = []
    seen: set[tuple[int, int, Fraction]] = set()
    for line in plan.seeds:
        for coordinate, name in ((0, "m"), (1, "b")):
            symbolic: dict[int, Homogeneous] = {
                index: tuple(Polynomial(value) for value in point)  # type: ignore[assignment]
                for index, point in base.items()
            }
            variable_point = list(symbolic[line])
            variable_point[coordinate] = Polynomial((0, 1))
            symbolic[line] = tuple(variable_point)  # type: ignore[assignment]

            symbolic_points = build_points(symbolic, plan)
            nonzero = [
                residual
                for residual in incidence_residuals(symbolic_points, triple_points)
                if isinstance(residual, Polynomial) and not residual.is_zero
            ]
            if len(nonzero) != 1 or nonzero[0].degree != 1:
                continue
            polynomial = nonzero[0]
            coefficient = polynomial.coefficients[1]
            if not coefficient:
                continue
            root = -polynomial.coefficients[0] / coefficient
            key = (line, coordinate, root)
            if key in seen:
                continue
            seen.add(key)

            adjusted = dict(base)
            point = list(adjusted[line])
            point[coordinate] = root
            adjusted[line] = tuple(point)  # type: ignore[assignment]
            candidates.append((adjusted, line, name))
    return candidates


def normalize_points(points: dict[int, Homogeneous], n: int) -> tuple[Line, ...]:
    lines: list[Line] = []
    for line in range(n):
        point = points[line]
        if any(isinstance(value, Polynomial) for value in point):
            raise RationalizationError("symbolic value remained after construction")
        x, y, scale = point
        if not scale:
            raise RationalizationError(f"dual point {line} lies at infinity")
        lines.append((x / scale, y / scale))  # type: ignore[operator]
    return tuple(lines)


def fraction_text(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def implicit_line(line: Line) -> ImplicitLine:
    """Convert ``y = m*x + b`` to ``a*x + b*y = c`` exactly."""

    slope, intercept = line
    return slope, Fraction(-1), -intercept


def implicit_geometry_rows(
    lines: Sequence[ImplicitLine],
) -> tuple[EventRows, frozenset[tuple[int, int]]]:
    """Build exact event rows directly from implicit line equations."""

    n = len(lines)
    buckets: list[dict[Fraction, set[int]]] = [dict() for _ in range(n)]
    parallel_pairs: set[tuple[int, int]] = set()

    for i, j in combinations(range(n), 2):
        a_i, b_i, c_i = lines[i]
        a_j, b_j, c_j = lines[j]
        determinant = a_i * b_j - a_j * b_i
        if not determinant:
            if a_i * c_j == a_j * c_i and b_i * c_j == b_j * c_i:
                raise CertificateError(f"coincident lines {i} and {j}")
            parallel_pairs.add((i, j))
            continue

        x = (c_i * b_j - b_i * c_j) / determinant
        y = (a_i * c_j - c_i * a_j) / determinant
        parameter_i = b_i * x - a_i * y
        parameter_j = b_j * x - a_j * y
        buckets[i].setdefault(parameter_i, set()).add(j)
        buckets[j].setdefault(parameter_j, set()).add(i)

    rows = tuple(
        tuple(frozenset(by_parameter[value]) for value in sorted(by_parameter))
        for by_parameter in buckets
    )
    return rows, frozenset(parallel_pairs)


def compare_implicit_realization(
    lines: Sequence[ImplicitLine], expected: WordStructure
) -> tuple[bool, str]:
    try:
        actual_rows, actual_parallel = implicit_geometry_rows(lines)
    except CertificateError as exc:
        return False, str(exc)

    if actual_parallel != expected.parallel_pairs:
        missing = sorted(expected.parallel_pairs - actual_parallel)
        extra = sorted(actual_parallel - expected.parallel_pairs)
        return False, f"parallel pairs differ: missing={missing[:3]} extra={extra[:3]}"
    if any(
        actual != wanted and actual != tuple(reversed(wanted))
        for wanted, actual in zip(expected.rows, actual_rows)
    ):
        return False, first_row_difference(expected.rows, actual_rows)
    return True, ""


def display_path(path: Path, repo_root: Path) -> Path:
    try:
        return path.resolve().relative_to(repo_root.resolve())
    except ValueError:
        return path


def realization_complexity(lines: Sequence[Line]) -> tuple[int, int]:
    values = [value for line in lines for value in line]
    max_bits = max(
        max(value.numerator.bit_length(), value.denominator.bit_length())
        for value in values
    )
    encoded_size = sum(len(fraction_text(value)) for value in values)
    return max_bits, encoded_size


def denominator_candidates(max_denominator: int) -> tuple[int, ...]:
    """Prefer tiny fractions, then grow geometrically for sensitive examples."""

    dense_limit = 512
    candidates = list(range(1, min(max_denominator, dense_limit) + 1))
    value = 768
    while value < max_denominator:
        candidates.append(value)
        value = value * 3 // 2
    if max_denominator > dense_limit:
        candidates.append(max_denominator)
    return tuple(dict.fromkeys(candidates))


def simple_rationalize_lines(
    source_lines: Sequence[Line],
    structure: WordStructure,
    *,
    max_denominator: int = 5000,
) -> RationalizationResult:
    """Round free coefficients and impose the linear incidence constraints."""

    for denominator in denominator_candidates(max_denominator):
        rounded = tuple(
            (
                slope.limit_denominator(denominator),
                intercept.limit_denominator(denominator),
            )
            for slope, intercept in source_lines
        )
        try:
            lines = quick_rationalize(rounded, structure)
        except CheckError:
            continue
        matches, _ = compare_realization(lines, structure)
        if matches:
            method = (
                "coefficient-rounding-and-linear-projection"
                if structure.multiple_points or structure.parallel_pairs
                else "coefficient-rounding"
            )
            return RationalizationResult(
                lines=lines,
                max_seed_denominator=denominator,
                seeds=(),
                adjusted_line=None,
                adjusted_coordinate=None,
                method=method,
            )
    raise RationalizationError(
        f"no simple exact realization found with coefficient denominators up to "
        f"{max_denominator}"
    )


def rationalize_lines(
    source_lines: Sequence[Line],
    structure: WordStructure,
    *,
    max_seed_denominator: int = 200,
) -> RationalizationResult:
    n = len(source_lines)
    if structure.parallel_pairs:
        raise RationalizationError(
            "full dual construction currently expects no parallel pairs"
        )
    if not structure.multiple_points:
        raise RationalizationError("configuration has no multiple points")
    if any(len(point) != 3 for point in structure.multiple_points):
        raise RationalizationError("only triple points are supported")

    plans = construction_plans(n, structure.multiple_points)
    if not plans:
        raise RationalizationError("incidence graph has no supported construction plan")

    for denominator in range(1, max_seed_denominator + 1):
        valid: list[tuple[tuple[int, int], RationalizationResult]] = []
        for plan in plans:
            try:
                candidates = candidate_seed_points(
                    source_lines,
                    structure.multiple_points,
                    plan,
                    denominator,
                )
            except (KeyError, ZeroDivisionError):
                continue
            for seeds, adjusted_line, adjusted_coordinate in candidates:
                try:
                    points = build_points(seeds, plan)
                    if any(incidence_residuals(points, structure.multiple_points)):
                        continue
                    lines = normalize_points(points, n)
                except (KeyError, ZeroDivisionError, RationalizationError):
                    continue
                matches, _ = compare_realization(lines, structure)
                if not matches:
                    continue
                result = RationalizationResult(
                    lines=lines,
                    max_seed_denominator=denominator,
                    seeds=plan.seeds,
                    adjusted_line=adjusted_line,
                    adjusted_coordinate=adjusted_coordinate,
                )
                valid.append((realization_complexity(lines), result))
        if valid:
            return min(valid, key=lambda item: item[0])[1]

    raise RationalizationError(
        f"no exact realization found with seed denominators up to "
        f"{max_seed_denominator}"
    )


def rationalize_for_certificate(
    source_lines: Sequence[Line],
    structure: WordStructure,
    *,
    max_simple_denominator: int = 5000,
    max_seed_denominator: int = 200,
) -> RationalizationResult:
    """Use the simple method when applicable, otherwise the dual construction."""

    try:
        source_repaired = quick_rationalize(source_lines, structure)
        simple_applicable, _ = compare_realization(source_repaired, structure)
    except CheckError:
        simple_applicable = False
    if simple_applicable:
        return simple_rationalize_lines(
            source_lines, structure, max_denominator=max_simple_denominator
        )
    return rationalize_lines(
        source_lines,
        structure,
        max_seed_denominator=max_seed_denominator,
    )


def load_source(source: Path) -> tuple[dict[str, object], int, WordStructure, tuple[Line, ...]]:
    series_match = SERIES_RE.fullmatch(source.parent.name)
    if series_match is None:
        raise RationalizationError(f"cannot derive n from {source.parent.name!r}")
    n = int(series_match.group(1))
    try:
        with source.open(encoding="utf-8") as stream:
            data = json.load(stream, parse_float=Decimal)
    except (OSError, json.JSONDecodeError) as exc:
        raise RationalizationError(f"cannot read {source}: {exc}") from exc
    if not isinstance(data, dict):
        raise RationalizationError("source must contain a JSON object")
    structure = replay_word(data.get("gens"), n)
    lines = parse_lines(data, n)
    return data, n, structure, lines


def make_certificate(
    source: Path,
    repo_root: Path,
    result: RationalizationResult,
    data: dict[str, object],
    n: int,
    structure: WordStructure,
) -> dict[str, object]:
    source_relative = source.resolve().relative_to(repo_root.resolve()).as_posix()
    return {
        "version": 1,
        "source": source_relative,
        "source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "gens": data["gens"],
        "n": n,
        "triangle_count": count_triangles(structure.rows),
        "line_equation": CERTIFICATE_LINE_EQUATION,
        "lines_frac": [
            [fraction_text(a), fraction_text(b), fraction_text(c)]
            for a, b, c in map(implicit_line, result.lines)
        ],
        "declared_triple_points": [
            sorted(point) for point in structure.multiple_points
        ],
        "declared_parallel_pairs": [
            list(pair) for pair in sorted(structure.parallel_pairs)
        ],
    }


def parse_fraction(value: object, field: str) -> Fraction:
    if not isinstance(value, str) or not value:
        raise CertificateError(f"{field}: expected a fraction string")
    try:
        result = Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise CertificateError(f"{field}: invalid fraction {value!r}") from exc
    if fraction_text(result) != value:
        raise CertificateError(f"{field}: fraction is not in canonical form")
    return result


def verify_certificate_payload(
    payload: object, source: Path, repo_root: Path
) -> tuple[ImplicitLine, ...]:
    if not isinstance(payload, dict):
        raise CertificateError("certificate must contain a JSON object")
    if payload.get("version") != 1:
        raise CertificateError("unsupported certificate version")

    expected_source = source.resolve().relative_to(repo_root.resolve()).as_posix()
    if payload.get("source") != expected_source:
        raise CertificateError(
            f"source mismatch: expected {expected_source!r}, got {payload.get('source')!r}"
        )
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if payload.get("source_sha256") != digest:
        raise CertificateError("source_sha256 does not match the source file")

    data, n, structure, _ = load_source(source)
    if payload.get("gens") != data.get("gens"):
        raise CertificateError("gens does not match the source file")
    if payload.get("n") != n:
        raise CertificateError(f"n: expected {n}, got {payload.get('n')!r}")
    triangle_count = count_triangles(structure.rows)
    if payload.get("triangle_count") != triangle_count:
        raise CertificateError(
            f"triangle_count: expected {triangle_count}, "
            f"got {payload.get('triangle_count')!r}"
        )
    if payload.get("line_equation") != CERTIFICATE_LINE_EQUATION:
        raise CertificateError("unsupported line_equation")

    raw_lines = payload.get("lines_frac")
    if not isinstance(raw_lines, list) or len(raw_lines) != n:
        raise CertificateError(f"lines_frac: expected {n} lines")
    lines: list[ImplicitLine] = []
    for index, raw_line in enumerate(raw_lines):
        if not isinstance(raw_line, list) or len(raw_line) != 3:
            raise CertificateError(f"lines_frac[{index}]: expected [a, b, c]")
        line = (
            parse_fraction(raw_line[0], f"lines_frac[{index}][0]"),
            parse_fraction(raw_line[1], f"lines_frac[{index}][1]"),
            parse_fraction(raw_line[2], f"lines_frac[{index}][2]"),
        )
        if not line[0] and not line[1]:
            raise CertificateError(
                f"lines_frac[{index}]: a and b cannot both be zero"
            )
        lines.append(line)

    expected_triples = [sorted(point) for point in structure.multiple_points]
    if payload.get("declared_triple_points") != expected_triples:
        raise CertificateError("declared_triple_points does not match gens")
    expected_parallel = [list(pair) for pair in sorted(structure.parallel_pairs)]
    if payload.get("declared_parallel_pairs") != expected_parallel:
        raise CertificateError("declared_parallel_pairs does not match gens")

    matches, reason = compare_implicit_realization(tuple(lines), structure)
    if not matches:
        raise CertificateError(f"exact realization does not match gens: {reason}")
    return tuple(lines)


def verify_certificate(certificate: Path, certificate_root: Path, repo_root: Path) -> Path:
    try:
        relative = certificate.resolve().relative_to(certificate_root.resolve())
    except ValueError as exc:
        raise CertificateError(f"certificate is outside {certificate_root}") from exc
    if relative.suffix != ".json":
        raise CertificateError("certificate name must end in .json")
    source_relative = relative
    source = repo_root / "gallery" / "data" / source_relative
    if not source.is_file():
        raise CertificateError(f"mirrored source does not exist: {source}")
    try:
        with certificate.open(encoding="utf-8") as stream:
            payload = json.load(stream)
    except (OSError, json.JSONDecodeError) as exc:
        raise CertificateError(f"cannot read certificate: {exc}") from exc
    verify_certificate_payload(payload, source, repo_root)
    return source


def classify_file(arguments: tuple[str, str, str]):
    return check_one_file(*arguments)


def required_sources(
    data_root: Path, certificate_root: Path, jobs: int
) -> list[Path]:
    sources, catalog_errors = validate_catalogs(data_root)
    if catalog_errors:
        raise RationalizationError("; ".join(catalog_errors))
    arguments = [
        (str(source), str(data_root), str(certificate_root)) for source in sources
    ]
    if jobs == 1:
        results = [classify_file(item) for item in arguments]
    else:
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            results = list(executor.map(classify_file, arguments, chunksize=8))
    wanted = {
        result.path
        for result in results
        if result.status in {"certificate-required", "certificate-present"}
    }
    return [source for source in sources if source.relative_to(data_root).as_posix() in wanted]


def first_sources_without_parallel_pairs(data_root: Path) -> list[Path]:
    """Select the first (lowest-ratio) entry of every plain N catalog."""

    directories = sorted(
        (
            path
            for path in data_root.iterdir()
            if path.is_dir() and path.name.isdigit()
        ),
        key=lambda path: int(path.name),
    )
    sources: list[Path] = []
    for directory in directories:
        catalog_path = data_root / f"{directory.name}.json"
        try:
            with catalog_path.open(encoding="utf-8") as stream:
                catalog = json.load(stream)
        except (OSError, json.JSONDecodeError) as exc:
            raise RationalizationError(f"cannot read {catalog_path}: {exc}") from exc
        if not isinstance(catalog, list) or not catalog or not isinstance(catalog[0], str):
            raise RationalizationError(f"{catalog_path}: expected a non-empty string list")
        relative = Path(catalog[0] + ".json")
        if relative.parent.name != directory.name:
            raise RationalizationError(
                f"{catalog_path}: first entry belongs to {relative.parent}"
            )
        source = data_root / relative
        if not source.is_file():
            raise RationalizationError(f"catalog source does not exist: {source}")
        sources.append(source)
    return sources


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sources", nargs="*", type=Path, help="source JSON files")
    parser.add_argument(
        "--all-required",
        action="store_true",
        help="find and rationalize every quick-check certificate fallback",
    )
    parser.add_argument(
        "--first-per-n",
        action="store_true",
        help="rationalize the first catalog entry for every series without -K",
    )
    parser.add_argument(
        "--max-seed-denominator",
        type=int,
        default=200,
        help="maximum seed approximation denominator (default: 200)",
    )
    parser.add_argument(
        "--max-simple-denominator",
        type=int,
        default=5000,
        help="maximum rounded coefficient denominator (default: 5000)",
    )
    parser.add_argument(
        "--certificates",
        type=Path,
        default=repo_root / "gallery" / "certificates",
        help="certificate output root",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(8, os.cpu_count() or 1),
        help="workers used by --all-required",
    )
    parser.add_argument(
        "--dry-run", action="store_true", help="construct and verify without writing"
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    repo_root = Path(__file__).resolve().parents[1]
    data_root = repo_root / "gallery" / "data"
    certificate_root = args.certificates.resolve()
    if (
        args.max_seed_denominator < 1
        or args.max_simple_denominator < 1
        or args.jobs < 1
    ):
        print("denominator and jobs must be positive", file=sys.stderr)
        return 2
    selection_count = bool(args.sources) + args.all_required + args.first_per_n
    if selection_count > 1:
        print("pass sources, --all-required, or --first-per-n", file=sys.stderr)
        return 2
    if selection_count == 0:
        print("pass sources, --all-required, or --first-per-n", file=sys.stderr)
        return 2

    try:
        if args.all_required:
            sources = required_sources(data_root, certificate_root, args.jobs)
        elif args.first_per_n:
            sources = first_sources_without_parallel_pairs(data_root)
        else:
            sources = [source.resolve() for source in args.sources]
    except RationalizationError as exc:
        print(f"[error] {exc}", file=sys.stderr)
        return 1

    failures = 0
    for source in sources:
        try:
            data, n, structure, source_lines = load_source(source)
            result = rationalize_for_certificate(
                source_lines,
                structure,
                max_simple_denominator=args.max_simple_denominator,
                max_seed_denominator=args.max_seed_denominator,
            )
            payload = make_certificate(
                source, repo_root, result, data, n, structure
            )
            verify_certificate_payload(payload, source, repo_root)
            destination = certificate_path(source, data_root, certificate_root)
            if not args.dry_run:
                destination.parent.mkdir(parents=True, exist_ok=True)
                destination.write_text(
                    json.dumps(payload, ensure_ascii=False, indent=2) + "\n",
                    encoding="utf-8",
                )
            max_bits, encoded_size = realization_complexity(result.lines)
            action = "verified" if args.dry_run else "wrote"
            print(
                f"[{action}] {display_path(destination, repo_root)}: "
                f"method={result.method}, "
                f"denominator={result.max_seed_denominator}, "
                f"max_bits={max_bits}, fraction_chars={encoded_size}"
            )
        except (CheckError, OSError, ValueError) as exc:
            failures += 1
            print(f"[failed] {source}: {exc}", file=sys.stderr)

    print(f"Rationalized {len(sources) - failures}/{len(sources)} configurations")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
