#!/usr/bin/env python3
"""Fast, exact smoke-check for every source configuration in the gallery.

The checker treats JSON decimal numbers as exact rationals.  It first compares
the stored realization with the event rows reconstructed from ``gens``.  If
that fails because the word declares parallel or triple intersections, a
small deterministic repair is attempted:

* slopes in every expected parallel class are replaced by their exact mean;
* with those slopes fixed, intercepts are projected onto the exact linear
  concurrency constraints of the triple points.

The repaired realization is accepted only if every exact event row agrees
with ``gens``.  No repaired coordinates are written.  A failed repair means
that a persistent certificate is required at the mirrored path under
``gallery/certificates``.

Only Python's standard library is required.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import asdict, dataclass
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from itertools import combinations
from math import lcm
from pathlib import Path
from typing import Iterable, Sequence


GENERATOR_RE = re.compile(r"^(\d+)(\*)?$")
SERIES_RE = re.compile(r"^(\d+)(?:-(\d+))?$")


class CheckError(ValueError):
    """A source configuration is structurally invalid."""


Event = frozenset[int]
EventRows = tuple[tuple[Event, ...], ...]
Line = tuple[Fraction, Fraction]  # y = slope*x + intercept


@dataclass(frozen=True)
class WordStructure:
    rows: EventRows
    parallel_pairs: frozenset[tuple[int, int]]
    multiple_points: tuple[frozenset[int], ...]


@dataclass(frozen=True)
class CheckResult:
    path: str
    series: str
    n: int
    triangle_count: int | None
    status: str
    reason: str = ""
    certificate: str = ""
    elapsed_ms: int = 0


def parse_decimal(value: object, *, field: str) -> Fraction:
    if isinstance(value, bool):
        raise CheckError(f"{field}: boolean is not a number")
    try:
        decimal = value if isinstance(value, Decimal) else Decimal(str(value))
    except (InvalidOperation, ValueError) as exc:
        raise CheckError(f"{field}: invalid number {value!r}") from exc
    if not decimal.is_finite():
        raise CheckError(f"{field}: number must be finite")
    return Fraction(decimal)


def parse_lines(data: dict[str, object], n: int) -> tuple[Line, ...]:
    raw_lines = data.get("lines")
    if not isinstance(raw_lines, list) or len(raw_lines) != n:
        size = len(raw_lines) if isinstance(raw_lines, list) else "not a list"
        raise CheckError(f"lines: expected {n} entries, got {size}")

    lines: list[Line] = []
    for index, raw_line in enumerate(raw_lines):
        if not isinstance(raw_line, list) or len(raw_line) != 2:
            raise CheckError(f"lines[{index}]: expected [slope, intercept]")
        slope = parse_decimal(raw_line[0], field=f"lines[{index}][0]")
        intercept = parse_decimal(raw_line[1], field=f"lines[{index}][1]")
        lines.append((slope, intercept))
    return tuple(lines)


def replay_word(gens: object, n: int) -> WordStructure:
    if not isinstance(gens, str) or not gens.strip():
        raise CheckError("gens: expected a non-empty string")

    wires = list(range(n))
    rows: list[list[Event]] = [[] for _ in range(n)]
    met_pairs: set[tuple[int, int]] = set()
    multiple_points: set[frozenset[int]] = set()

    for position, token in enumerate(gens.split(), start=1):
        match = GENERATOR_RE.fullmatch(token)
        if match is None:
            raise CheckError(f"gens token {position}: unsupported token {token!r}")

        generator = int(match.group(1))
        is_triple = match.group(2) is not None
        width = 3 if is_triple else 2
        if generator < 0 or generator + width > n:
            raise CheckError(
                f"gens token {position}: {token!r} is outside 0..{n - width}"
            )

        block = wires[generator : generator + width]
        if block != sorted(block):
            raise CheckError(
                f"gens token {position}: lines {block} do not meet for the first time"
            )

        event_lines = frozenset(block)
        pairs = {(min(a, b), max(a, b)) for a, b in combinations(block, 2)}
        repeated = pairs & met_pairs
        if repeated:
            raise CheckError(
                f"gens token {position}: repeated line pair {sorted(repeated)[0]}"
            )
        met_pairs.update(pairs)

        for line in block:
            rows[line].append(frozenset(event_lines - {line}))

        if is_triple:
            multiple_points.add(event_lines)
            wires[generator], wires[generator + 2] = (
                wires[generator + 2],
                wires[generator],
            )
        else:
            wires[generator], wires[generator + 1] = (
                wires[generator + 1],
                wires[generator],
            )

    all_pairs = {(i, j) for i in range(n) for j in range(i + 1, n)}
    parallel_pairs = frozenset(all_pairs - met_pairs)
    return WordStructure(
        rows=tuple(tuple(row) for row in rows),
        parallel_pairs=parallel_pairs,
        multiple_points=tuple(sorted(multiple_points, key=lambda point: tuple(sorted(point)))),
    )


def geometry_rows(lines: Sequence[Line]) -> tuple[EventRows, frozenset[tuple[int, int]]]:
    n = len(lines)
    buckets: list[dict[Fraction, set[int]]] = [defaultdict(set) for _ in range(n)]
    parallel_pairs: set[tuple[int, int]] = set()

    for i, j in combinations(range(n), 2):
        slope_i, intercept_i = lines[i]
        slope_j, intercept_j = lines[j]
        if slope_i == slope_j:
            if intercept_i == intercept_j:
                raise CheckError(f"coincident lines {i} and {j}")
            parallel_pairs.add((i, j))
            continue

        x = (intercept_j - intercept_i) / (slope_i - slope_j)
        buckets[i][x].add(j)
        buckets[j][x].add(i)

    rows = tuple(
        tuple(frozenset(by_x[x]) for x in sorted(by_x))
        for by_x in buckets
    )
    return rows, frozenset(parallel_pairs)


def first_row_difference(expected: EventRows, actual: EventRows) -> str:
    for line, (expected_row, actual_row) in enumerate(zip(expected, actual)):
        if actual_row != expected_row and actual_row != tuple(reversed(expected_row)):
            def show(row: Sequence[Event]) -> str:
                return " ".join(
                    "{" + ",".join(str(v) for v in sorted(event)) + "}"
                    for event in row
                )

            return (
                f"line {line}: expected [{show(expected_row)}], "
                f"got [{show(actual_row)}]"
            )
    return "event rows differ"


def compare_realization(lines: Sequence[Line], expected: WordStructure) -> tuple[bool, str]:
    # Clear denominators once per line. Intersection x-coordinates can then
    # be compared by integer cross products, without reducing and sorting
    # hundreds of Fractions for each realization.
    coefficients = []
    for slope, intercept in lines:
        denominator = lcm(slope.denominator, intercept.denominator)
        coefficients.append((
            slope.numerator * (denominator // slope.denominator),
            intercept.numerator * (denominator // intercept.denominator),
            denominator,
        ))

    intersections: list[dict[int, tuple[int, int]]] = [{} for _ in lines]
    actual_parallel = set()
    for i, j in combinations(range(len(lines)), 2):
        ai, bi, di = coefficients[i]
        aj, bj, dj = coefficients[j]
        numerator = bj * di - bi * dj
        denominator = ai * dj - aj * di
        if denominator == 0:
            if numerator == 0:
                return False, f"coincident lines {i} and {j}"
            actual_parallel.add((i, j))
            continue
        if denominator < 0:
            numerator, denominator = -numerator, -denominator
        intersections[i][j] = intersections[j][i] = (numerator, denominator)

    if actual_parallel != expected.parallel_pairs:
        missing = sorted(expected.parallel_pairs - actual_parallel)
        extra = sorted(actual_parallel - expected.parallel_pairs)
        return False, f"parallel pairs differ: missing={missing[:3]} extra={extra[:3]}"
    # A pseudoline has no preferred orientation.  Depending on which side of
    # the affine chart its direction lies, the stored x-order is either the
    # generator row or its reverse.
    for line, row in enumerate(expected.rows):
        previous = None
        direction = 0
        for event in row:
            coordinate = intersections[line][next(iter(event))]
            numerator, denominator = coordinate
            if any(
                numerator * intersections[line][other][1]
                != intersections[line][other][0] * denominator
                for other in event
            ):
                return False, f"line {line}: declared event {sorted(event)} is split"
            if previous is not None:
                difference = numerator * previous[1] - previous[0] * denominator
                sign = (difference > 0) - (difference < 0)
                if sign == 0:
                    return False, f"line {line}: distinct events coincide"
                if direction and direction != sign:
                    return False, f"line {line}: event order differs"
                direction = sign
            previous = coordinate
    return True, ""


def solve_linear_system(
    matrix: Sequence[Sequence[Fraction]], rhs: Sequence[Fraction]
) -> list[Fraction]:
    """Solve a consistent, possibly singular system; free variables become zero."""

    row_count = len(matrix)
    if row_count == 0:
        return []
    column_count = len(matrix[0])
    augmented = [list(row) + [rhs[i]] for i, row in enumerate(matrix)]
    pivot_columns: list[int] = []
    pivot_row = 0

    for column in range(column_count):
        found = next(
            (row for row in range(pivot_row, row_count) if augmented[row][column]),
            None,
        )
        if found is None:
            continue
        augmented[pivot_row], augmented[found] = augmented[found], augmented[pivot_row]
        divisor = augmented[pivot_row][column]
        augmented[pivot_row] = [value / divisor for value in augmented[pivot_row]]

        for row in range(row_count):
            if row == pivot_row:
                continue
            factor = augmented[row][column]
            if factor:
                augmented[row] = [
                    value - factor * pivot
                    for value, pivot in zip(augmented[row], augmented[pivot_row])
                ]

        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == row_count:
            break

    for row in range(pivot_row, row_count):
        if not any(augmented[row][column] for column in range(column_count)):
            if augmented[row][-1]:
                raise CheckError("inconsistent rationalization system")

    solution = [Fraction(0) for _ in range(column_count)]
    for row, column in enumerate(pivot_columns):
        solution[column] = augmented[row][-1]
    return solution


def project_to_constraints(
    values: Sequence[Fraction], constraints: Sequence[Sequence[Fraction]]
) -> list[Fraction]:
    """Exact orthogonal projection of VALUES onto A*x=0."""

    if not constraints:
        return list(values)

    gram = [
        [sum(a * b for a, b in zip(left, right)) for right in constraints]
        for left in constraints
    ]
    residual = [sum(a * value for a, value in zip(row, values)) for row in constraints]
    multipliers = solve_linear_system(gram, residual)
    return [
        value
        - sum(constraints[row][column] * multipliers[row] for row in range(len(constraints)))
        for column, value in enumerate(values)
    ]


def parallel_components(
    n: int, parallel_pairs: Iterable[tuple[int, int]]
) -> tuple[tuple[int, ...], ...]:
    parents = list(range(n))

    def find(value: int) -> int:
        while parents[value] != value:
            parents[value] = parents[parents[value]]
            value = parents[value]
        return value

    def union(left: int, right: int) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    pairs = set(parallel_pairs)
    for left, right in pairs:
        union(left, right)

    groups: dict[int, list[int]] = defaultdict(list)
    for line in range(n):
        groups[find(line)].append(line)

    components = tuple(tuple(group) for group in groups.values() if len(group) > 1)
    for component in components:
        expected_pairs = {
            (min(left, right), max(left, right))
            for left, right in combinations(component, 2)
        }
        if not expected_pairs <= pairs:
            raise CheckError(f"non-transitive parallel relation in {component}")
    return components


def quick_rationalize(lines: Sequence[Line], expected: WordStructure) -> tuple[Line, ...]:
    n = len(lines)
    slopes = [line[0] for line in lines]
    intercepts = [line[1] for line in lines]

    for component in parallel_components(n, expected.parallel_pairs):
        shared_slope = sum((slopes[line] for line in component), Fraction(0)) / len(component)
        for line in component:
            slopes[line] = shared_slope

    for left, right in combinations(range(n), 2):
        if slopes[left] == slopes[right] and (left, right) not in expected.parallel_pairs:
            raise CheckError(f"rationalization made non-parallel lines {left} and {right} parallel")

    constraints: list[list[Fraction]] = []
    for point in expected.multiple_points:
        if len(point) != 3:
            raise CheckError(f"quick rationalizer supports triples only, got {sorted(point)}")
        first, second, third = sorted(point)
        row = [Fraction(0) for _ in range(n)]
        row[first] = slopes[second] - slopes[third]
        row[second] = slopes[third] - slopes[first]
        row[third] = slopes[first] - slopes[second]
        if not any(row):
            raise CheckError(f"triple point {sorted(point)} has one direction")
        constraints.append(row)

    repaired_intercepts = project_to_constraints(intercepts, constraints)
    return tuple(zip(slopes, repaired_intercepts))


def count_triangles(rows: EventRows) -> int:
    n = len(rows)
    positions: list[dict[int, int]] = [dict() for _ in range(n)]
    for line, row in enumerate(rows):
        for position, event in enumerate(row):
            for other in event:
                positions[line][other] = position

    count = 0
    for first, second, third in combinations(range(n), 3):
        if (
            second not in positions[first]
            or third not in positions[first]
            or third not in positions[second]
        ):
            continue
        if positions[first][second] == positions[first][third]:
            continue
        if (
            abs(positions[first][second] - positions[first][third]) == 1
            and abs(positions[second][first] - positions[second][third]) == 1
            and abs(positions[third][first] - positions[third][second]) == 1
        ):
            count += 1
    return count


def certificate_path(source: Path, data_root: Path, certificate_root: Path) -> Path:
    relative = source.relative_to(data_root)
    return certificate_root / relative


def check_one_file(source_name: str, data_root_name: str, certificate_root_name: str) -> CheckResult:
    started = time.perf_counter()
    source = Path(source_name)
    data_root = Path(data_root_name)
    certificate_root = Path(certificate_root_name)
    relative = source.relative_to(data_root)
    series = relative.parent.as_posix()
    n = 0

    def result(status: str, *, triangle_count: int | None = None, reason: str = "") -> CheckResult:
        cert = certificate_path(source, data_root, certificate_root)
        return CheckResult(
            path=relative.as_posix(),
            series=series,
            n=n,
            triangle_count=triangle_count,
            status=status,
            reason=reason,
            certificate=cert.relative_to(certificate_root.parent).as_posix(),
            elapsed_ms=round((time.perf_counter() - started) * 1000),
        )

    try:
        series_match = SERIES_RE.fullmatch(relative.parent.name)
        if series_match is None:
            raise CheckError(f"invalid series directory {relative.parent.name!r}")
        n = int(series_match.group(1))
        expected_parallel_count = int(series_match.group(2) or 0)

        with source.open(encoding="utf-8") as stream:
            data = json.load(stream, parse_float=Decimal)
        if not isinstance(data, dict):
            raise CheckError("top-level JSON value must be an object")

        structure = replay_word(data.get("gens"), n)
        if len(structure.parallel_pairs) != expected_parallel_count:
            raise CheckError(
                f"series declares {expected_parallel_count} parallel pairs, "
                f"gens imply {len(structure.parallel_pairs)}"
            )
        lines = parse_lines(data, n)
        triangle_count = count_triangles(structure.rows)

        matches, reason = compare_realization(lines, structure)
        if matches:
            return result("exact", triangle_count=triangle_count)

        try:
            repaired = quick_rationalize(lines, structure)
            repaired_matches, repaired_reason = compare_realization(repaired, structure)
        except CheckError as exc:
            repaired_matches, repaired_reason = False, str(exc)

        if repaired_matches:
            return result("quick-rationalized", triangle_count=triangle_count)

        cert = certificate_path(source, data_root, certificate_root)
        if cert.is_file():
            return result(
                "certificate-present",
                triangle_count=triangle_count,
                reason=repaired_reason or reason,
            )
        return result(
            "certificate-required",
            triangle_count=triangle_count,
            reason=repaired_reason or reason,
        )
    except (CheckError, json.JSONDecodeError, OSError) as exc:
        return result("invalid", reason=str(exc))


def validate_catalogs(data_root: Path) -> tuple[list[Path], list[str]]:
    errors: list[str] = []
    directories = sorted(
        path for path in data_root.iterdir() if path.is_dir() and SERIES_RE.fullmatch(path.name)
    )
    directory_names = {path.name for path in directories}

    try:
        with (data_root / "index.json").open(encoding="utf-8") as stream:
            index = json.load(stream)
        index_names = {
            item.get("id") for item in index if isinstance(item, dict) and isinstance(item.get("id"), str)
        }
        if index_names != directory_names:
            errors.append(
                "index.json series mismatch: "
                f"missing={sorted(directory_names - index_names)} "
                f"extra={sorted(index_names - directory_names)}"
            )
    except (OSError, json.JSONDecodeError, TypeError) as exc:
        errors.append(f"cannot read index.json: {exc}")

    sources: list[Path] = []
    for directory in directories:
        actual = {
            path.relative_to(data_root).with_suffix("").as_posix()
            for path in directory.glob("*.json")
        }
        list_path = data_root / f"{directory.name}.json"
        try:
            with list_path.open(encoding="utf-8") as stream:
                listed_value = json.load(stream)
            if not isinstance(listed_value, list) or not all(
                isinstance(item, str) for item in listed_value
            ):
                raise CheckError("expected an array of strings")
            listed = set(listed_value)
            if len(listed) != len(listed_value):
                errors.append(f"{list_path.name}: duplicate entries")
            if listed != actual:
                errors.append(
                    f"{list_path.name}: file list mismatch: "
                    f"missing={sorted(actual - listed)[:3]} extra={sorted(listed - actual)[:3]}"
                )
        except (OSError, json.JSONDecodeError, CheckError) as exc:
            errors.append(f"{list_path.name}: {exc}")
        sources.extend(sorted(directory.glob("*.json")))

    return sources, errors


def summarize(results: Sequence[CheckResult]) -> dict[str, int]:
    return dict(sorted(Counter(result.status for result in results).items()))


def write_report(
    report_path: Path,
    data_root: Path,
    certificate_root: Path,
    results: Sequence[CheckResult],
    catalog_errors: Sequence[str],
    elapsed_seconds: float,
) -> None:
    payload = {
        "version": 1,
        "data_root": data_root.as_posix(),
        "certificate_root": certificate_root.as_posix(),
        "elapsed_seconds": round(elapsed_seconds, 3),
        "summary": summarize(results),
        "catalog_errors": list(catalog_errors),
        "files": [asdict(result) for result in results],
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data",
        type=Path,
        default=repo_root / "gallery" / "data",
        help="gallery data directory",
    )
    parser.add_argument(
        "--certificates",
        type=Path,
        default=repo_root / "gallery" / "certificates",
        help="mirrored certificate directory",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=min(8, os.cpu_count() or 1),
        help="worker process count (default: min(8, CPU count))",
    )
    parser.add_argument("--report", type=Path, help="write a detailed JSON report")
    parser.add_argument(
        "--allow-certificate-required",
        action="store_true",
        help="inventory mode: do not fail solely because certificates are missing",
    )
    parser.add_argument(
        "--show",
        choices=("all", "problems", "none"),
        default="problems",
        help="which per-file results to print",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    data_root = args.data.resolve()
    certificate_root = args.certificates.resolve()
    if args.jobs < 1:
        print("--jobs must be at least 1", file=sys.stderr)
        return 2

    started = time.perf_counter()
    sources, catalog_errors = validate_catalogs(data_root)
    arguments = [
        (str(source), str(data_root), str(certificate_root))
        for source in sources
    ]

    if args.jobs == 1:
        results = [check_one_file(*item) for item in arguments]
    else:
        with ProcessPoolExecutor(max_workers=args.jobs) as executor:
            results = list(executor.map(_check_tuple, arguments, chunksize=8))

    results.sort(key=lambda result: result.path)
    elapsed = time.perf_counter() - started

    visible_statuses = {"invalid", "certificate-required", "certificate-present"}
    for result in results:
        if args.show == "none":
            break
        if args.show == "all" or result.status in visible_statuses:
            suffix = f": {result.reason}" if result.reason else ""
            print(f"[{result.status}] {result.path}{suffix}")

    for error in catalog_errors:
        print(f"[catalog-invalid] {error}")

    summary = summarize(results)
    print(
        f"Checked {len(results)} configurations in {elapsed:.2f}s: "
        + ", ".join(f"{status}={count}" for status, count in summary.items())
    )
    if catalog_errors:
        print(f"Catalog errors: {len(catalog_errors)}")

    if args.report:
        write_report(
            args.report,
            data_root,
            certificate_root,
            results,
            catalog_errors,
            elapsed,
        )
        print(f"Report: {args.report}")

    if catalog_errors or summary.get("invalid", 0):
        return 1
    if summary.get("certificate-required", 0) and not args.allow_certificate_required:
        return 1
    return 0


def _check_tuple(arguments: tuple[str, str, str]) -> CheckResult:
    return check_one_file(*arguments)


if __name__ == "__main__":
    raise SystemExit(main())
