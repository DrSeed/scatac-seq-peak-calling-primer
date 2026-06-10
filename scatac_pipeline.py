# Minimal scATAC-seq fragments-to-matrix teaching pipeline (simulated data).
import argparse
import random
from collections import defaultdict


def simulate_fragments(n_cells=200, frags_per_cell=300, n_bins=500, bin_size=1000, seed=0):
    # Generate fake fragments biased towards a handful of accessible regions.
    random.seed(seed)
    hot_bins = set(random.sample(range(n_bins), k=max(1, n_bins // 20)))
    fragments = []
    for c in range(n_cells):
        bc = "cell_" + str(c)
        n = max(1, int(random.gauss(frags_per_cell, frags_per_cell * 0.3)))
        for _ in range(n):
            if random.random() < 0.6:
                b = random.choice(list(hot_bins))
            else:
                b = random.randrange(n_bins)
            start = b * bin_size + random.randrange(bin_size - 50)
            fragments.append(("chr1", start, start + 50, bc))
    return fragments, bin_size


def filter_barcodes(fragments, min_frags=50):
    # Separate real cells from low-count noise barcodes.
    counts = defaultdict(int)
    for _, _, _, bc in fragments:
        counts[bc] += 1
    keep = {bc for bc, n in counts.items() if n >= min_frags}
    return keep, counts


def build_peaks(fragments, bin_size, min_support=10):
    # Candidate peaks are genome bins with enough total coverage.
    bin_counts = defaultdict(int)
    for _, start, _, _ in fragments:
        bin_counts[start // bin_size] += 1
    peaks = sorted(b for b, n in bin_counts.items() if n >= min_support)
    peak_index = {b: i for i, b in enumerate(peaks)}
    return peaks, peak_index


def build_matrix(fragments, keep_bc, peak_index, bin_size):
    # Sparse representation: dict of (cell, peak) -> count.
    mat = defaultdict(int)
    for _, start, _, bc in fragments:
        if bc not in keep_bc:
            continue
        b = start // bin_size
        if b in peak_index:
            mat[(bc, peak_index[b])] += 1
    return mat


def compute_qc(fragments, keep_bc, mat, counts):
    # Per-cell total fragments, in-peak fragments and FRiP.
    in_peak = defaultdict(int)
    for (bc, _), v in mat.items():
        in_peak[bc] += v
    rows = []
    for bc in sorted(keep_bc):
        total = counts[bc]
        ip = in_peak.get(bc, 0)
        frip = ip / total if total else 0.0
        rows.append((bc, total, ip, round(frip, 4)))
    return rows


def write_qc(rows, path):
    with open(path, "w") as fh:
        fh.write("barcode\ttotal_fragments\tin_peak_fragments\tfrip\n")
        for bc, total, ip, frip in rows:
            fh.write(bc + "\t" + str(total) + "\t" + str(ip) + "\t" + str(frip) + "\n")


def sparsity(mat, n_cells, n_peaks):
    if n_cells == 0 or n_peaks == 0:
        return 1.0
    nonzero = len(mat)
    return 1.0 - nonzero / (n_cells * n_peaks)


def main():
    ap = argparse.ArgumentParser(description="scATAC-seq teaching pipeline")
    ap.add_argument("--simulate", action="store_true", help="use simulated fragments")
    ap.add_argument("--out", default="results/", help="output directory")
    ap.add_argument("--min-frags", type=int, default=50)
    args = ap.parse_args()

    if not args.simulate:
        print("Only --simulate mode is implemented in this teaching scaffold.")
        return

    fragments, bin_size = simulate_fragments()
    keep_bc, counts = filter_barcodes(fragments, min_frags=args.min_frags)
    peaks, peak_index = build_peaks(fragments, bin_size)
    mat = build_matrix(fragments, keep_bc, peak_index, bin_size)
    rows = compute_qc(fragments, keep_bc, mat, counts)

    out = args.out.rstrip("/") + "/cell_qc.tsv"
    write_qc(rows, out)

    spars = sparsity(mat, len(keep_bc), len(peaks))
    print("Cells kept: " + str(len(keep_bc)))
    print("Peaks called: " + str(len(peaks)))
    print("Matrix sparsity: " + str(round(spars * 100, 2)) + "%")
    print("QC written to: " + out)
    print("Reminder: that sparsity is normal for scATAC-seq. Do not fight it with magic.")


if __name__ == "__main__":
    main()
