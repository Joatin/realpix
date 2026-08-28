#!/usr/bin/env python3
"""Generate HEALPix golden test vectors with healpy (the reference C++ implementation).

Run with:

    uv run --with healpy,numpy python tools/gen_golden.py

The generated CSVs are committed under tests/data/ so the Rust test suite needs neither
Python nor a network connection. Pass --check to regenerate into a temporary directory and
diff against the committed files.
"""
from __future__ import annotations

import argparse
import os
import sys
import tempfile

import numpy as np
import healpy as hp

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "tests", "data")

DEPTHS = [0, 1, 2, 3, 5, 8, 12, 16, 20, 29]
TWO_PI = 2.0 * np.pi
HALF_PI = np.pi / 2.0


def f(x: float) -> str:
    """Shortest representation that round-trips exactly."""
    return repr(float(x))


def positions(n: int, seed: int) -> list[tuple[float, float]]:
    """A fixed set of (lon, lat) samples: uniform on the sphere plus explicit edge cases."""
    rng = np.random.default_rng(seed)
    lon = rng.uniform(0.0, TWO_PI, n)
    lat = np.arcsin(rng.uniform(-1.0, 1.0, n))
    out = [(float(a), float(b)) for a, b in zip(lon, lat)]

    # Edge cases: poles, the belt/cap transition, base-cell longitude boundaries, and
    # positions a hair away from each of those.
    transition = float(np.arcsin(2.0 / 3.0))
    lats = [
        0.0, HALF_PI, -HALF_PI, transition, -transition,
        HALF_PI - 1e-9, -HALF_PI + 1e-9, HALF_PI - 1e-15, -HALF_PI + 1e-15,
        transition + 1e-12, transition - 1e-12, 1e-15, -1e-15,
        0.5, -0.5, 1.4, -1.4,
    ]
    lons = [
        0.0, 1e-15, HALF_PI, np.pi, 3.0 * HALF_PI, TWO_PI - 1e-12,
        np.pi / 4.0, 3.0 * np.pi / 4.0, 0.1, 6.0,
    ]
    for la in lats:
        for lo in lons:
            out.append((float(lo), float(la)))
    return out


def cells(depth: int, n: int, seed: int) -> np.ndarray:
    npix = 12 * 4 ** depth
    if npix <= n:
        return np.arange(npix, dtype=np.int64)
    rng = np.random.default_rng(seed)
    special = [0, 1, 2, 3, npix - 1, npix - 2, npix // 2, npix // 2 - 1,
               2 * (2 ** depth) * (2 ** depth - 1) - 1,  # last cell of the RING north cap
               4 ** depth, 4 ** depth - 1, 11 * 4 ** depth]
    special = [c for c in special if 0 <= c < npix]
    rnd = rng.integers(0, npix, n - len(special), dtype=np.int64)
    return np.unique(np.concatenate([np.array(special, dtype=np.int64), rnd]))


def write(name: str, header: str, rows) -> None:
    path = os.path.join(OUT, name)
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for row in rows:
            fh.write(",".join(row) + "\n")
    print(f"wrote {path}")


def gen_ang2pix() -> None:
    """theta/phi, (lon, lat) and unit-vector queries, with the reference answer for each.

    `stable` marks samples that are not within 1e-13 rad of a cell boundary. Only those can
    be compared against `hash(lon, lat)`, because that entry point computes z as sin(lat)
    while the reference computes it as cos(theta): the two differ by an ulp, which is
    enough to land on the other side of a boundary for a sample sitting exactly on one.
    The theta/phi and unit-vector entry points take the reference's own inputs and are
    compared exactly for every sample.
    """
    rows = []
    for depth in DEPTHS:
        nside = 1 << depth
        for lon, lat in positions(300, 1234 + depth):
            theta = min(max(HALF_PI - lat, 0.0), np.pi)
            nest = int(hp.ang2pix(nside, theta, lon, nest=True))
            ring = int(hp.ang2pix(nside, theta, lon, nest=False))

            eps = 1e-13
            stable = 1
            for dlat, dlon in ((eps, 0.0), (-eps, 0.0), (0.0, eps), (0.0, -eps)):
                t = min(max(HALF_PI - (lat + dlat), 0.0), np.pi)
                if int(hp.ang2pix(nside, t, lon + dlon, nest=True)) != nest:
                    stable = 0
                    break

            # Unit vector built exactly the way realpix builds it.
            x = np.cos(lat) * np.cos(lon)
            y = np.cos(lat) * np.sin(lon)
            z = np.sin(lat)
            vnest = int(hp.vec2pix(nside, x, y, z, nest=True))
            vring = int(hp.vec2pix(nside, x, y, z, nest=False))
            rows.append((str(depth), f(lon), f(lat), f(theta), str(nest), str(ring),
                         f(x), f(y), f(z), str(vnest), str(vring), str(stable)))
    write("ang2pix.csv",
          "depth,lon,lat,theta,nest,ring,x,y,z,vec_nest,vec_ring,stable", rows)


def gen_pix2ang() -> None:
    rows = []
    for depth in DEPTHS:
        nside = 1 << depth
        cs = cells(depth, 200, 99 + depth)
        tn, pn = hp.pix2ang(nside, cs, nest=True)
        tr, pr = hp.pix2ang(nside, cs, nest=False)
        for c, a, b, u, v in zip(cs, tn, pn, tr, pr):
            rows.append((str(depth), str(int(c)),
                         f(b), f(HALF_PI - a), f(v), f(HALF_PI - u)))
    write("pix2ang.csv", "depth,cell,nest_lon,nest_lat,ring_lon,ring_lat", rows)


def gen_scheme_conversion() -> None:
    rows = []
    for depth in DEPTHS:
        nside = 1 << depth
        cs = cells(depth, 400, 7 + depth)
        ring = hp.nest2ring(nside, cs)
        for c, r in zip(cs, ring):
            rows.append((str(depth), str(int(c)), str(int(r))))
    write("nest2ring.csv", "depth,nest,ring", rows)


def gen_neighbours() -> None:
    rows = []
    for depth in DEPTHS:
        nside = 1 << depth
        cs = cells(depth, 200, 555 + depth)
        # healpy returns SW, W, NW, N, NE, E, SE, S; -1 marks a missing neighbour.
        nb = hp.get_all_neighbours(nside, cs, nest=True)
        for i, c in enumerate(cs):
            rows.append(tuple([str(depth), str(int(c))] + [str(int(nb[k][i])) for k in range(8)]))
    write("neighbours.csv", "depth,cell,sw,w,nw,n,ne,e,se,s", rows)


def gen_boundaries() -> None:
    rows = []
    for depth in [0, 1, 2, 3, 5, 8, 12]:
        nside = 1 << depth
        cs = cells(depth, 100, 31 + depth)
        for c in cs:
            # step=1 -> the four corners, in the order N, W, S, E.
            v = hp.boundaries(nside, int(c), step=1, nest=True)
            flat = [f(v[axis][i]) for i in range(4) for axis in range(3)]
            rows.append(tuple([str(depth), str(int(c))] + flat))
    write("boundaries.csv",
          "depth,cell," + ",".join(f"{n}{a}" for n in "nwse" for a in "xyz"), rows)


def gen_query_disc() -> None:
    rng = np.random.default_rng(2024)
    rows = []
    for depth in [3, 5, 6, 8]:
        nside = 1 << depth
        for _ in range(8):
            lon = float(rng.uniform(0, TWO_PI))
            lat = float(np.arcsin(rng.uniform(-1, 1)))
            radius = float(rng.uniform(0.005, 0.25))
            vec = hp.ang2vec(HALF_PI - lat, lon)
            pix = hp.query_disc(nside, vec, radius, inclusive=False, nest=True)
            rows.append((str(depth), f(lon), f(lat), f(radius),
                         " ".join(str(int(p)) for p in np.sort(pix))))
    # A few degenerate radii too.
    for depth, radius in [(3, 1e-6), (5, 1.0), (2, 2.0)]:
        nside = 1 << depth
        lon, lat = 1.0, 0.3
        vec = hp.ang2vec(HALF_PI - lat, lon)
        pix = hp.query_disc(nside, vec, radius, inclusive=False, nest=True)
        rows.append((str(depth), f(lon), f(lat), f(radius),
                     " ".join(str(int(p)) for p in np.sort(pix))))
    write("query_disc.csv", "depth,lon,lat,radius,pixels", rows)


def main() -> int:
    global OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="regenerate into a temp dir and diff against the committed files")
    args = ap.parse_args()

    committed = OUT
    if args.check:
        OUT = tempfile.mkdtemp(prefix="realpix-golden-")

    os.makedirs(OUT, exist_ok=True)
    gen_ang2pix()
    gen_pix2ang()
    gen_scheme_conversion()
    gen_neighbours()
    gen_boundaries()
    gen_query_disc()

    if args.check:
        bad = 0
        for name in sorted(os.listdir(OUT)):
            a, b = os.path.join(committed, name), os.path.join(OUT, name)
            if not os.path.exists(a) or open(a).read() != open(b).read():
                print(f"MISMATCH: {name}", file=sys.stderr)
                bad += 1
        if bad:
            return 1
        print("all golden files match")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
