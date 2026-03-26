# mgrs.py
import math

# --- WGS84 ---
A = 6378137.0
F = 1 / 298.257223563
K0 = 0.9996

E = math.sqrt(F * (2 - F))
E2 = E**2
E_PRIME2 = E2 / (1 - E2)

# Latitude bands
LAT_BANDS = "CDEFGHJKLMNPQRSTUVWX"

# 100k grid sets
EASTING_SETS = ["ABCDEFGH", "JKLMNPQR", "STUVWXYZ"]
NORTHING_SETS = ["ABCDEFGHJKLMNPQRSTUV", "FGHJKLMNPQRSTUVABCDE"]

# Minimum northing per band (meters)
BAND_MIN_NORTHING = {
    'C': 1100000, 'D': 2000000, 'E': 2800000, 'F': 3700000,
    'G': 4600000, 'H': 5500000, 'J': 6400000, 'K': 7300000,
    'L': 8200000, 'M': 9100000,
    'N': 0, 'P': 800000, 'Q': 1700000, 'R': 2600000,
    'S': 3500000, 'T': 4400000, 'U': 5300000,
    'V': 6200000, 'W': 7000000, 'X': 7900000,
}


# ---------------------------
# Helpers
# ---------------------------
def latlon_to_zone(lat, lon):
    zone = int((lon + 180) / 6) + 1

    # Norway
    if 56 <= lat < 64 and 3 <= lon < 12:
        zone = 32

    # Svalbard
    if 72 <= lat <= 84:
        if 0 <= lon < 9:
            zone = 31
        elif 9 <= lon < 21:
            zone = 33
        elif 21 <= lon < 33:
            zone = 35
        elif 33 <= lon < 42:
            zone = 37

    return zone


def latitude_band(lat):
    if lat < -80 or lat > 84:
        raise ValueError("Latitude out of MGRS bounds")

    if lat >= 72:
        return 'X'

    return LAT_BANDS[int((lat + 80) / 8)]


def band_latitude_range(band):
    idx = LAT_BANDS.index(band)
    lower = -80 + idx * 8
    upper = 84 if band == 'X' else lower + 8
    return lower, upper


def truncate(v, precision):
    factor = 10 ** (5 - precision)
    return int(v // factor)


# ---------------------------
# LAT/LON → UTM
# ---------------------------
def latlon_to_utm(lat, lon):
    zone = latlon_to_zone(lat, lon)
    northern = lat >= 0

    lat_rad = math.radians(lat)
    lon_rad = math.radians(lon)

    lon_origin = (zone - 1) * 6 - 180 + 3
    lon_origin_rad = math.radians(lon_origin)

    N = A / math.sqrt(1 - E2 * math.sin(lat_rad)**2)
    T = math.tan(lat_rad)**2
    C = E_PRIME2 * math.cos(lat_rad)**2
    A_ = math.cos(lat_rad) * (lon_rad - lon_origin_rad)

    M = A * (
        (1 - E2/4 - 3*E2**2/64 - 5*E2**3/256) * lat_rad
        - (3*E2/8 + 3*E2**2/32 + 45*E2**3/1024) * math.sin(2*lat_rad)
        + (15*E2**2/256 + 45*E2**3/1024) * math.sin(4*lat_rad)
        - (35*E2**3/3072) * math.sin(6*lat_rad)
    )

    easting = K0 * N * (
        A_ + (1 - T + C) * A_**3 / 6 +
        (5 - 18*T + T**2 + 72*C - 58*E_PRIME2) * A_**5 / 120
    ) + 500000

    northing = K0 * (
        M + N * math.tan(lat_rad) * (
            A_**2 / 2 +
            (5 - T + 9*C + 4*C**2) * A_**4 / 24 +
            (61 - 58*T + T**2 + 600*C - 330*E_PRIME2) * A_**6 / 720
        )
    )

    if not northern:
        northing += 10000000

    return zone, northern, easting, northing


# ---------------------------
# UTM → LAT/LON
# ---------------------------
def utm_to_latlon(zone, northern, easting, northing):
    x = easting - 500000
    y = northing if northern else northing - 10000000

    lon_origin = (zone - 1) * 6 - 180 + 3

    M = y / K0
    mu = M / (A * (1 - E2/4 - 3*E2**2/64 - 5*E2**3/256))

    e1 = (1 - math.sqrt(1 - E2)) / (1 + math.sqrt(1 - E2))

    J1 = 3*e1/2 - 27*e1**3/32
    J2 = 21*e1**2/16 - 55*e1**4/32
    J3 = 151*e1**3/96
    J4 = 1097*e1**4/512

    fp = mu + J1*math.sin(2*mu) + J2*math.sin(4*mu) + \
         J3*math.sin(6*mu) + J4*math.sin(8*mu)

    C1 = E_PRIME2 * math.cos(fp)**2
    T1 = math.tan(fp)**2
    N1 = A / math.sqrt(1 - E2 * math.sin(fp)**2)
    R1 = N1 * (1 - E2) / (1 - E2 * math.sin(fp)**2)
    D = x / (N1 * K0)

    lat = fp - (N1 * math.tan(fp) / R1) * (
        D**2/2 - (5 + 3*T1 + 10*C1 - 4*C1**2 - 9*E_PRIME2)*D**4/24 +
        (61 + 90*T1 + 298*C1 + 45*T1**2 - 252*E_PRIME2 - 3*C1**2)*D**6/720
    )

    lon = (
        D - (1 + 2*T1 + C1)*D**3/6 +
        (5 - 2*C1 + 28*T1 - 3*C1**2 + 8*E_PRIME2 + 24*T1**2)*D**5/120
    ) / math.cos(fp)

    return math.degrees(lat), lon_origin + math.degrees(lon)


# ---------------------------
# LAT/LON → MGRS
# ---------------------------
def latlon_to_mgrs(lat, lon, precision=5):
    zone, northern, easting, northing = latlon_to_utm(lat, lon)
    band = latitude_band(lat)

    col_set = (zone - 1) % 3
    row_set = (zone - 1) % 2

    col = int(easting // 100000)
    row = int(northing // 100000) % 20

    col_letter = EASTING_SETS[col_set][col - 1]
    row_letter = NORTHING_SETS[row_set][row]

    e = truncate(easting % 100000, precision)
    n = truncate(northing % 100000, precision)

    fmt = f"{{:0{precision}d}}"
    return f"{zone}{band}{col_letter}{row_letter}{fmt.format(e)}{fmt.format(n)}"


# ---------------------------
# MGRS → LAT/LON
# ---------------------------
def mgrs_to_latlon(mgrs):
    mgrs = mgrs.strip().replace(" ", "").upper()

    i = 0

    # zone: 1 or 2 digits
    if len(mgrs) < 5:
        raise ValueError("Invalid MGRS string")

    if mgrs[1].isdigit():
        zone = int(mgrs[:2])
        i = 2
    else:
        zone = int(mgrs[:1])
        i = 1

    band = mgrs[i]
    i += 1

    col_letter = mgrs[i]
    row_letter = mgrs[i+1]
    i += 2

    digits = mgrs[i:]
    precision = len(digits) // 2

    if len(digits) % 2 != 0:
        raise ValueError("MGRS numeric component must have even length")

    cell = 10 ** (5 - precision)
    e = int(digits[:precision]) * cell
    n = int(digits[precision:]) * cell

    col_set = (zone - 1) % 3
    row_set = (zone - 1) % 2

    e += (EASTING_SETS[col_set].index(col_letter) + 1) * 100000
    n += NORTHING_SETS[row_set].index(row_letter) * 100000
    e += cell / 2
    n += cell / 2

    northern = band >= 'N'

    # --- NGA band alignment ---
    min_n = BAND_MIN_NORTHING[band]
    while n + cell / 2 < min_n:
        n += 2000000

    lower_lat, upper_lat = band_latitude_range(band)

    # Search the 2,000 km northing cycles until the latitude falls inside the band.
    for _ in range(6):
        lat, lon = utm_to_latlon(zone, northern, e, n)
        in_band = lower_lat <= lat <= upper_lat if band == 'X' else lower_lat <= lat < upper_lat
        if in_band:
            return lat, lon

        # Only cells whose center falls just outside the band need the extra
        # edge checks; this preserves the band-boundary fix without making the
        # common case pay for three UTM inversions.
        south_lat, _ = utm_to_latlon(zone, northern, e, n - cell / 2)
        north_lat, _ = utm_to_latlon(zone, northern, e, n + cell / 2)
        cell_min_lat = min(south_lat, north_lat)
        cell_max_lat = max(south_lat, north_lat)

        if cell_max_lat >= lower_lat and cell_min_lat <= upper_lat:
            return lat, lon

        if lat < lower_lat:
            n += 2000000
        else:
            n -= 2000000

    raise ValueError(f"Could not align MGRS northing for band {band}")




# Test cases
def test():
    tests = [
        (0.0, -177.0),       # zone 1 regression
        (24.0000053576329, -120.13017064707174),  # band-boundary regression
        (41.0082, 28.9784),   # Istanbul
        (39.9208, 32.8541),   # Ankara
        (48.8582, 2.2945),    # Eiffel
        (40.6892, -74.0445),  # Liberty
        (-33.8568, 151.2153), # Sydney
        (51.5074, -0.1278),   # London
        (35.6895, 139.6917),  # Tokyo
        (77.5536, 23.6702), # Svakbard
    ]

    for lat, lon in tests:
        mgrs = latlon_to_mgrs(lat, lon, 5)
        lat2, lon2 = mgrs_to_latlon(mgrs)

        print(mgrs, lat2, lon2)

        assert abs(lat - lat2) < 1e-5
        assert abs(lon - lon2) < 2e-5


## create a bechmark function for 1 million mgrs conversion
import random
import time

# --- reuse your functions ---
# latlon_to_mgrs
# mgrs_to_latlon

def random_latlon():
    # MGRS valid range (UTM coverage)
    lat = random.uniform(-79.9, 83.9)
    lon = random.uniform(-180, 180)
    return lat, lon


def benchmark_mgrs(n=1_000_000, precision=5, validate=False):
    print(f"Running benchmark with {n:,} points (precision={precision})")

    # --- generate dataset ---
    coords = [random_latlon() for _ in range(n)]

    # --- forward: lat/lon → MGRS ---
    t0 = time.perf_counter()
    mgrs_list = [latlon_to_mgrs(lat, lon, precision) for lat, lon in coords]
    t1 = time.perf_counter()

    forward_time = t1 - t0
    forward_rate = n / forward_time

    print(f"\nForward (LL → MGRS):")
    print(f"  Time: {forward_time:.3f} s")
    print(f"  Rate: {forward_rate:,.0f} ops/sec")

    # --- reverse: MGRS → lat/lon ---
    t2 = time.perf_counter()
    coords2 = [mgrs_to_latlon(m) for m in mgrs_list]
    t3 = time.perf_counter()

    reverse_time = t3 - t2
    reverse_rate = n / reverse_time

    print(f"\nReverse (MGRS → LL):")
    print(f"  Time: {reverse_time:.3f} s")
    print(f"  Rate: {reverse_rate:,.0f} ops/sec")

    # --- validation (optional, slower) ---
    if validate:
        max_lat_err = 0.0
        max_lon_err = 0.0

        for (lat1, lon1), (lat2, lon2) in zip(coords, coords2):
            max_lat_err = max(max_lat_err, abs(lat1 - lat2))
            max_lon_err = max(max_lon_err, abs(lon1 - lon2))

        print(f"\nValidation:")
        print(f"  Max lat error: {max_lat_err:.8f}")
        print(f"  Max lon error: {max_lon_err:.8f}")

    # --- summary ---
    total_time = (t1 - t0) + (t3 - t2)
    print(f"\nTotal time: {total_time:.3f} s")
    print(f"Overall throughput: {2*n / total_time:,.0f} ops/sec")


if __name__ == "__main__":
    benchmark_mgrs(n=1_000_000, precision=5, validate=False)
    #test()
