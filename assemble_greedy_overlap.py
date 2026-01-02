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
) -> Tuple[str, List[str]]:
    """
    Greedy assembly into one contig by repeatedly merging the best overlap.
    Returns (contig, leftovers).
    """
    # Start with the longest sequence as seed
    seqs = [s for s in seqs if s and set(s) <= set("ACGTN")]
    seqs.sort(key=len, reverse=True)
    if not seqs:
        raise SystemExit("No sequences loaded.")

    contig = seqs.pop(0)

    rounds = 0
    while seqs and rounds < max_rounds:
        rounds += 1
        best_i = None
        best_k = 0
        best_oriented = None  # the sequence to merge (maybe revcomp)
        best_side = None      # "right" for suffix(contig)->prefix(seq) or "left" for suffix(seq)->prefix(contig)

        # Find best merge candidate
        for i, s in enumerate(seqs):
            candidates = [s]
            if consider_revcomp:
                candidates.append(revcomp(s))

            for cand in candidates:
                # Merge to the right: contig + cand
                k = best_overlap(contig, cand, min_ovl)
                if k > best_k:
                    best_k = k
                    best_i = i
                    best_oriented = cand
                    best_side = "right"

                # Merge to the left: cand + contig
                k2 = best_overlap(cand, contig, min_ovl)
                if k2 > best_k:
                    best_k = k2
                    best_i = i
                    best_oriented = cand
                    best_side = "left"

        if best_i is None or best_k < min_ovl:
            # can't extend further
            break

        # Perform merge
        original = seqs.pop(best_i)
        if best_side == "right":
            contig = contig + best_oriented[best_k:]
        else:
            contig = best_oriented + contig[best_k:]

    return contig, seqs

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
    contig, leftovers = greedy_assemble(
        seqs,
        min_ovl=args.min_ovl,
        consider_revcomp=not args.no_revcomp,
    )

    out = Path(args.out)
    write_fasta(out, f"greedy_denovo_minovl{args.min_ovl}_len{len(contig)}", contig)

    print(f"[OK] Wrote contig: {out} (len={len(contig)})")
    print(f"[INFO] Leftover sequences not merged: {len(leftovers)}")
    if leftovers:
        lens = sorted([len(x) for x in leftovers], reverse=True)[:10]
        print(f"[INFO] Top leftover lengths: {lens}")

if __name__ == "__main__":
    main()
