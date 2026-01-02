#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Tuple, Optional

"""
How to run greedy denovo overlap

python3 assemble_greedy_overlap.py --products-dir products --min-ovl 125 --out products/denovo_greedy.fasta

"""

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}

_rc_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(seq: str) -> str:
    return seq.translate(_rc_map)[::-1]

def read_fasta(path: Path) -> List[Tuple[str, str]]:
    recs = []
    header = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    recs.append((header, "".join(seq).upper()))
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line)
        if header is not None:
            recs.append((header, "".join(seq).upper()))
    return recs

def load_sequences_from_dir(products_dir: Path, prefix: str = "product_") -> List[Tuple[str, str]]:
    fasta_files = sorted(
        [p for p in products_dir.iterdir() if p.suffix.lower() in FASTA_EXTENSIONS],
        key=lambda p: p.name
    )
    seqs: List[Tuple[str, str]] = []
    for fp in fasta_files:
        if prefix and not fp.name.startswith(prefix):
            continue
        for h, s in read_fasta(fp):
            seqs.append((h, s))
    return seqs

def load_sequences_from_fasta(fasta_path: Path) -> List[Tuple[str, str]]:
    return read_fasta(fasta_path)

def best_overlap(a: str, b: str, min_ovl: int, max_ovl: Optional[int] = None) -> int:
    """
    Return best suffix(a) == prefix(b) overlap length >= min_ovl, else 0.
    Tries longer overlaps first for speed by checking candidate lengths.
    """
    if max_ovl is None:
        max_ovl = min(len(a), len(b))
    max_ovl = min(max_ovl, len(a), len(b))
    # check from longest to shortest
    for k in range(max_ovl, min_ovl - 1, -1):
        if a[-k:] == b[:k]:
            return k
    return 0

def greedy_assemble(
    seqs: List[str],
    min_ovl: int,
    consider_revcomp: bool = True,
    max_rounds: int = 100000
) -> Tuple[str, List[str], List[str]]:
    """
    Greedy assembly into one contig by repeatedly merging the best overlap.
    Returns (contig, leftovers, used_sequences_original_orientation).
    """
    seqs = [s for s in seqs if s and set(s) <= set("ACGTN")]
    seqs.sort(key=len, reverse=True)
    if not seqs:
        raise SystemExit("No sequences loaded.")

    used: List[str] = []
    contig = seqs.pop(0)
    used.append(contig)

    rounds = 0
    while seqs and rounds < max_rounds:
        rounds += 1
        best_i = None
        best_k = 0
        best_oriented = None
        best_side = None

        for i, s in enumerate(seqs):
            candidates = [s]
            if consider_revcomp:
                candidates.append(revcomp(s))

            for cand in candidates:
                k = best_overlap(contig, cand, min_ovl)
                if k > best_k:
                    best_k, best_i, best_oriented, best_side = k, i, cand, "right"

                k2 = best_overlap(cand, contig, min_ovl)
                if k2 > best_k:
                    best_k, best_i, best_oriented, best_side = k2, i, cand, "left"

        if best_i is None or best_k < min_ovl:
            break

        original = seqs.pop(best_i)
        used.append(original)

        if best_side == "right":
            contig = contig + best_oriented[best_k:]
        else:
            contig = best_oriented + contig[best_k:]

    return contig, seqs, used

def assemble_all(seqs: List[str], min_ovl: int, consider_revcomp: bool, max_contigs: int = 1000):
    contigs = []
    remaining = seqs[:]
    for _ in range(max_contigs):
        if not remaining:
            break
        contig, remaining, used = greedy_assemble(remaining, min_ovl, consider_revcomp=consider_revcomp)
        contigs.append(contig)
        # optional: stop if contig is “trivial” (only seed used)
        if len(used) <= 1:
            break
    return contigs, remaining


def write_fasta(path: Path, header: str, seq: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

def main():
    ap = argparse.ArgumentParser(description="Greedy reference-free overlap assembly of PCR products.")
    ap.add_argument("--products-dir", help="Folder with product_*.fasta files.")
    ap.add_argument("--fasta", help="Alternative: a single FASTA containing all products.")
    ap.add_argument("--min-ovl", type=int, default=125, help="Minimum required overlap (bp). Default 125.")
    ap.add_argument("--no-revcomp", action="store_true", help="Disable reverse-complement matching.")
    ap.add_argument("--out", default="products/denovo_greedy.fasta", help="Output contig FASTA.")
    ap.add_argument("--prefix", default="product_", help="Only read files starting with this prefix (dir mode).")
    args = ap.parse_args()

    if not args.products_dir and not args.fasta:
        raise SystemExit("Provide either --products-dir or --fasta")

    if args.products_dir:
        seq_records = load_sequences_from_dir(Path(args.products_dir), prefix=args.prefix)
    else:
        seq_records = load_sequences_from_fasta(Path(args.fasta))

    seqs = [s for _, s in seq_records]
    contigs, leftovers = assemble_all(
        seqs,
        min_ovl=args.min_ovl,
        consider_revcomp=not args.no_revcomp,
        max_contigs=1000,
    )

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    with out.open("w") as f:
        for idx, c in enumerate(contigs, 1):
            f.write(f">contig_{idx}_minovl{args.min_ovl}_len{len(c)}\n")
            for i in range(0, len(c), 80):
                f.write(c[i:i+80] + "\n")

    print(f"[OK] Wrote {len(contigs)} contigs to {out}")
    print(f"[INFO] Leftover sequences not merged: {len(leftovers)}")

    if leftovers:
        lens = sorted([len(x) for x in leftovers], reverse=True)[:10]
        print(f"[INFO] Top leftover lengths: {lens}")

if __name__ == "__main__":
    main()
