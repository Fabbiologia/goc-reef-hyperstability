#!/usr/bin/env python3
"""
00a_download_oisst.py

Downloads genuine NOAA OISST v2.1 daily sea surface temperature for the Gulf
of California at the product's native quarter degree resolution.

WHY THIS EXISTS. The series previously held in data/env was three sample
points per one degree latitude band, thirty points for the whole basin,
which discards the spatial structure that matters here: the cold pool the
midriff islands maintain by tidal mixing, and the difference between the
peninsula and mainland sides. The real product carries about 812 ocean cells
over the same box.

METHOD. Each daily global file on the NCEI archive is subset server side
through OPeNDAP so that only the Gulf box crosses the network, about 13 kB
per day instead of 1.7 MB. Files are written one per day as compact CSV and
the run is resumable: a day already on disk is skipped, so the download can
be interrupted and restarted freely.

Grid indexing. OISST latitudes are -89.875 + 0.25 i and longitudes are
0.125 + 0.25 j in degrees east, so the Gulf box (22 to 32 N, 244 to 255 E,
equivalently 116 to 105 W) is i 447 to 487 and j 975 to 1019.

Output: ../data/env/oisst_daily/YYYY/oisst_gulf_YYYYMMDD.csv
        columns lat, lon, sst  (degrees north, degrees east, Celsius)

Usage:  python3 00a_download_oisst.py [--start 1981-09-01] [--end 2026-04-30]
                                      [--workers 4]
"""
import argparse
import csv
import datetime as dt
import os
import re
import sys
import time
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE = ("https://www.ncei.noaa.gov/thredds/dodsC/OisstBase/NetCDF/V2.1/AVHRR/"
        "{ym}/oisst-avhrr-v02r01.{ymd}.nc.ascii")
# Gulf of California box in OISST grid indices
I0, I1 = 447, 487          # latitude  22.125 .. 32.125 N
J0, J1 = 975, 1019         # longitude 244.125 .. 255.125 E  (-115.875 .. -104.875)
QUERY = "?sst%5B0:1:0%5D%5B0:1:0%5D%5B{i0}:1:{i1}%5D%5B{j0}:1:{j1}%5D"
FILL = -999                # OISST missing-value flag in scaled integers
SCALE = 0.01               # scaled integer -> Celsius

OUT_ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "..", "..", "data", "env", "oisst_daily")


def lat_of(i):
    return -89.875 + 0.25 * i


def lon_of(j):
    return 0.125 + 0.25 * j


def parse_ascii(text):
    """Pull the sst grid out of an OPeNDAP .ascii response.

    The body lists one line per latitude row, each beginning with bracketed
    indices then the comma separated values for that row.
    """
    if "sst" not in text:
        return None
    body = text.split("sst.sst", 1)[-1]
    rows = []
    # A data row is exactly "[t][z][lat], v, v, ...". The dimension header
    # that opens the body, "[1][1][41][45]", also begins with a bracket but
    # carries no comma; matching it as data would insert a phantom row and
    # shift every latitude by one cell, so the pattern is anchored strictly.
    data_row = re.compile(r"^\[\d+\]\[\d+\]\[\d+\],")
    for line in body.splitlines():
        line = line.strip()
        if not data_row.match(line):
            continue
        payload = line.split("]", 3)[-1]
        nums = re.findall(r"-?\d+", payload)
        if nums:
            rows.append([int(n) for n in nums])
    if not rows:
        return None
    return rows


def day_path(day):
    return os.path.join(OUT_ROOT, day.strftime("%Y"),
                        "oisst_gulf_%s.csv" % day.strftime("%Y%m%d"))


def fetch_day(day, retries=4):
    """Download and write one day. Returns (day, status)."""
    out = day_path(day)
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return day, "skip"
    url = BASE.format(ym=day.strftime("%Y%m"), ymd=day.strftime("%Y%m%d")) + \
        QUERY.format(i0=I0, i1=I1, j0=J0, j1=J1)
    delay = 2.0
    for attempt in range(retries):
        try:
            # curl rather than urllib: urllib stalls against this host in
            # this environment, curl negotiates it without trouble.
            res = subprocess.run(
                ["curl", "-sS", "--max-time", "120", "--retry", "0", url],
                capture_output=True, text=True)
            if res.returncode != 0:
                time.sleep(delay); delay *= 1.8; continue
            text = res.stdout
            if "Not Found" in text[:400] or "Error {" in text[:200]:
                return day, "missing"
            rows = parse_ascii(text)
            if rows is None:
                return day, "nodata"
            os.makedirs(os.path.dirname(out), exist_ok=True)
            tmp = out + ".part"
            n = 0
            with open(tmp, "w", newline="") as fh:
                w = csv.writer(fh)
                w.writerow(["lat", "lon", "sst"])
                for ri, row in enumerate(rows):
                    for cj, v in enumerate(row):
                        if v == FILL:
                            continue
                        w.writerow([round(lat_of(I0 + ri), 3),
                                    round(lon_of(J0 + cj), 3),
                                    round(v * SCALE, 2)])
                        n += 1
            os.replace(tmp, out)
            return day, ("ok" if n else "empty")
        except Exception:
            time.sleep(delay); delay *= 1.8
    return day, "fail"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--start", default="1981-09-01")
    ap.add_argument("--end", default="2026-04-30")
    ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()

    d0 = dt.date.fromisoformat(a.start)
    d1 = dt.date.fromisoformat(a.end)
    days = [d0 + dt.timedelta(days=k) for k in range((d1 - d0).days + 1)]
    todo = [d for d in days if not (os.path.exists(day_path(d)) and
                                    os.path.getsize(day_path(d)) > 0)]
    print("OISST v2.1, Gulf of California box at native 0.25 degree")
    print("  range   : %s to %s (%d days)" % (a.start, a.end, len(days)))
    print("  already : %d days on disk" % (len(days) - len(todo)))
    print("  to fetch: %d days with %d workers" % (len(todo), a.workers),
          flush=True)

    counts = {}
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=a.workers) as ex:
        futs = {ex.submit(fetch_day, d): d for d in todo}
        for f in as_completed(futs):
            _, status = f.result()
            counts[status] = counts.get(status, 0) + 1
            done += 1
            if done % 250 == 0 or done == len(todo):
                el = time.time() - t0
                rate = done / el if el else 0
                left = (len(todo) - done) / rate if rate else 0
                print("  %6d/%d  %.1f/s  eta %.0f min  %s"
                      % (done, len(todo), rate, left / 60, counts), flush=True)
    print("done:", counts, flush=True)
    if counts.get("fail"):
        print("  %d days failed; re-run to retry them" % counts["fail"])
        sys.exit(1)


if __name__ == "__main__":
    main()
