# Exact realization certificates

Certificates mirror viewer source paths below `gallery/data`:

```text
gallery/data/20/example.json
gallery/certificates/20/example.json
```

Certificates are included for selected arrangements to simplify independent verification,
and for arrangements which cannot pass the quick rationalization.

Version 1 stores:

- the mirrored source path and its SHA-256 digest;
- a copy of `gens`, `n`, and the exact triangle count;
- rational `[a, b, c]` triples for equations `a*x + b*y = c`;
- the triple points and parallel pairs reconstructed from `gens`.

Run `python3 verification/verify_certificates.py` from the repository root to
verify all certificates using exact rational arithmetic.

Regenerate one showcase certificate for the first (lowest-`ratio`) entry of
every plain `N` catalog with:

```bash
python3 verification/rationalize.py --first-per-n
```
