# mgrs

Pure Python MGRS <-> latitude/longitude conversion with no external dependencies.

This repository provides a single-file library for:

- converting WGS84 latitude/longitude to MGRS
- converting MGRS back to latitude/longitude
- running round-trip regression tests
- benchmarking conversion throughput

## Why this project

- Pure Python: no C extensions, compiled wheels, or third-party packages
- Small surface area: one file, direct functions, easy to vendor
- Fast enough for bulk conversion workloads

## Installation

Install from PyPI:

```bash
python3 -m pip install mgrs-pure-python
```

Or clone the repository and import `mgrs.py` directly:

```bash
git clone https://github.com/aydink/mgrs.git
cd mgrs
```

No dependency installation is required.

To install it as a local package with the CLI entrypoint:

```bash
python3 -m pip install .
```

## Library usage

### Encode latitude/longitude to MGRS

```python
from mgrs import latlon_to_mgrs

mgrs_ref = latlon_to_mgrs(41.0082, 28.9784, precision=5)
print(mgrs_ref)
# 35TPF6637041552
```

### Decode MGRS to latitude/longitude

```python
from mgrs import mgrs_to_latlon

lat, lon = mgrs_to_latlon("35TPF6637041552")
print(lat, lon)
# 41.008200115451004 28.97839994400028
```

### Round-trip conversion

```python
from mgrs import latlon_to_mgrs, mgrs_to_latlon

lat, lon = 40.6892, -74.0445
mgrs_ref = latlon_to_mgrs(lat, lon, precision=5)
lat2, lon2 = mgrs_to_latlon(mgrs_ref)

print(mgrs_ref)
print(lat2, lon2)
```

## Command-line usage

If installed from PyPI, use:

```bash
mgrs --lat 41.0082 --lon 28.9784 --precision 5
```

The examples below also work directly from the repo checkout:

### Encode coordinates

```bash
python3 mgrs.py --lat 41.0082 --lon 28.9784 --precision 5
```

### Decode MGRS

```bash
python3 mgrs.py --mgrs 35TPF6637041552
```

### Run regression tests

```bash
python3 mgrs.py --test
```

### Run a benchmark

```bash
python3 mgrs.py --benchmark --count 1000000 --precision 5
```

### Benchmark with validation

```bash
python3 mgrs.py --benchmark --count 100000 --precision 5 --validate
```

## API

- `latlon_to_mgrs(lat, lon, precision=5) -> str`
- `mgrs_to_latlon(mgrs) -> tuple[float, float]`
- `benchmark_mgrs(n=1_000_000, precision=5, validate=False) -> None`
- `test() -> None`

## Notes

- Coordinates use the WGS84 datum.
- MGRS precision is controlled with `precision=1..5`.
- Decoding returns the center of the represented MGRS cell.

## License

MIT. See [LICENSE](LICENSE).

## PyPI Publishing

This repository includes a GitHub Actions workflow for PyPI Trusted Publishing.
After configuring the project on PyPI, publishing a GitHub release will upload
the package automatically.
