# Gallery quick verification

`quick_check.py` performs a fast, exact smoke-check of every source JSON in
`gallery/data` using only the Python standard library.

It validates the data catalog, replays `gens`, counts bounded triangular faces
combinatorially, and compares the exact event order of the stored line
coefficients with the event order implied by `gens`.
JSON decimals are read as exact rational values, without conversion
through binary floating point. Event order is checked with integer cross
products after clearing coefficient denominators, avoiding a sort of rational
intersection coordinates for each source.

When the stored values do not represent declared parallel or triple events
exactly, the checker attempts a deliberately limited repair. Parallel slopes
are made equal and intercepts are projected onto the exact linear concurrency
constraints while slopes remain fixed. The candidate is accepted only when
all event rows agree with `gens`.

No repaired coordinates are written. If the quick repair fails, the source
needs a persistent certificate at the mirrored path:

```text
gallery/data/20/example.json
gallery/certificates/20/example.json
```

## Usage

Strict mode, suitable for CI:

```bash
python3 verification/quick_check.py
```

Initial inventory, allowing files which still require certificates:

```bash
python3 verification/quick_check.py \
  --allow-certificate-required \
  --report quick-check-report.json
```

Useful options:

```text
--jobs N                 number of worker processes
--show all|problems|none per-file console output
--report PATH            detailed machine-readable report
```

Run the small unit-test suite with:

```bash
python3 -m unittest -v verification.test_quick_check
```

The quick checker does **not** verify the contents of certificate JSON files.
CI must run the full verifier as another required job.

The detailed report records a triangle count and one of four useful statuses
for every source: `exact`, `quick-rationalized`, `certificate-required`, or
`certificate-present`. Structurally malformed inputs use the `invalid` status.

## Full rationalization

`rationalize.py` handles the configurations which the limited repair cannot
solve. It works in the dual projective plane: source lines become points and
triple intersections become collinear triples. The script selects rational
seed points, constructs the dependent points with exact cross products, and,
when necessary, solves one remaining linear incidence constraint exactly.

Generate every certificate currently requested by the quick check:

```bash
python3 verification/rationalize.py --all-required
```

Generate selected certificates or inspect the result without writing:

```bash
python3 verification/rationalize.py gallery/data/10/example.json
python3 verification/rationalize.py --all-required --dry-run
```

Generate one showcase certificate for the first (lowest-`ratio`) entry of
every plain `N` catalog:

```bash
python3 verification/rationalize.py --first-per-n
```

For configurations accepted by the quick method, coefficients are rounded to
small fractions before the exact linear incidence constraints are imposed.
The dual construction is the fallback for harder triple-point configurations;
it currently expects no parallel pairs. Both methods verify the complete
unoriented event rows before writing a certificate, and unsupported incidence
graphs fail explicitly.

Verify all stored certificates independently of the generator:

```bash
python3 verification/verify_certificates.py
```

## Independent verification

`verify_direct_exact.py` and `verify_events.py` are adapted from Andrea
Maiorana's independent exact verifiers:

https://github.com/rufio72/kobon_triangles_k14

The local copies read the expected result from the certificate field
`triangle_count` instead of the original field `count`; their verification
algorithms are unchanged.

Run both independent checks for one certificate with:

```bash
python3 verification/verify_direct_exact.py gallery/certificates/14/14-1jc3dhxbbh727.json
python3 verification/verify_events.py gallery/certificates/14/14-1jc3dhxbbh727.json
```

The quick check and the certificate verifier are intended to be separate CI
jobs. The former ensures that every difficult source has a certificate; the
latter verifies each certificate's source hash, metadata, rational
coefficients, incidences, event order, and triangle count.
