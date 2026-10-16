#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

"""
Greedy de novo overlap assembly from FASTA or FASTQ input.

Examples:
  python3 greedy_multi.py --products-dir products --min-ovl 125 --out products/denovo_greedy.fasta
  python3 greedy_multi.py --fastq reads.fastq --min-ovl 125 --out denovo_greedy.fasta
"""

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}
FASTQ_EXTENSIONS = {".fq", ".fastq"}
SEQUENCE_EXTENSIONS = FASTA_EXTENSIONS | FASTQ_EXTENSIONS

_rc_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_rc_map)[::-1]


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    recs: List[Tuple[str, str]] = []
    header = None
    seq: List[str] = []
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


def read_fastq(path: Path) -> List[Tuple[str, str]]:
    """Read standard four-line FASTQ records and discard quality values."""
    recs: List[Tuple[str, str]] = []
    with path.open() as f:
        record_number = 0
        while True:
            header_line = f.readline()
            if not header_line:
                break
            if not header_line.strip():
                continue

            record_number += 1
            sequence_line = f.readline()
            separator_line = f.readline()
            quality_line = f.readline()
            if not sequence_line or not separator_line or not quality_line:
                raise ValueError(f"{path}: incomplete FASTQ record {record_number}")

            header = header_line.rstrip("\r\n")
            sequence = sequence_line.strip()
            separator = separator_line.rstrip("\r\n")
            quality = quality_line.rstrip("\r\n")

            if not header.startswith("@"):
                raise ValueError(
                    f"{path}: FASTQ record {record_number} header does not start with '@'"
                )
            if not separator.startswith("+"):
                raise ValueError(
                    f"{path}: FASTQ record {record_number} separator does not start with '+'"
                )
            if len(sequence) != len(quality):
                raise ValueError(
                    f"{path}: FASTQ record {record_number} has sequence length "
                    f"{len(sequence)} but quality length {len(quality)}"
                )
            recs.append((header[1:].strip(), sequence.upper()))
    return recs


def read_sequences(path: Path) -> List[Tuple[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path.resolve()}")
    suffix = path.suffix.lower()
    if suffix in FASTA_EXTENSIONS:
        return read_fasta(path)
    if suffix in FASTQ_EXTENSIONS:
        return read_fastq(path)
    supported = ", ".join(sorted(SEQUENCE_EXTENSIONS))
    raise ValueError(f"Unsupported input extension for {path}; expected {supported}")


def load_sequences_from_dir(products_dir: Path, prefix: str = "product_") -> List[Tuple[str, str]]:
    if not products_dir.is_dir():
        raise NotADirectoryError(f"Input directory not found: {products_dir.resolve()}")
    sequence_files = sorted(
        (p for p in products_dir.iterdir() if p.suffix.lower() in SEQUENCE_EXTENSIONS),
        key=lambda p: p.name,
    )
    seqs: List[Tuple[str, str]] = []
    for path in sequence_files:
        if prefix and not path.name.startswith(prefix):
            continue
        seqs.extend(read_sequences(path))
    return seqs


def best_overlap(
    a: str,
    b: str,
    min_ovl: int,
    min_identity: float,
    max_ovl: Optional[int] = None,
) -> Tuple[int, float]:
    """Return the longest acceptable ungapped overlap and its identity."""
    if max_ovl is None:
        max_ovl = min(len(a), len(b))
    max_ovl = min(max_ovl, len(a), len(b))
    for k in range(max_ovl, min_ovl - 1, -1):
        matches = sum(left == right for left, right in zip(a[-k:], b[:k]))
        identity = matches / k
        if identity >= min_identity:
            return k, identity
    return 0, 0.0


def greedy_assemble(
    seqs: List[str],
    min_ovl: int,
    min_identity: float,
    consider_revcomp: bool = True,
    max_rounds: int = 100000,
) -> Tuple[str, List[str], List[str]]:
    seqs = [s for s in seqs if s and set(s) <= set("ACGTN")]
    seqs.sort(key=len, reverse=True)
    if not seqs:
        raise SystemExit("No valid sequences loaded.")

    used = []
    contig = seqs.pop(0)
    used.append(contig)
    rounds = 0

    while seqs and rounds < max_rounds:
        rounds += 1
        best_i = None
        best_k = 0
        best_identity = 0.0
        best_oriented = None
        best_side = None

        for i, sequence in enumerate(seqs):
            candidates = [sequence]
            if consider_revcomp:
                candidates.append(revcomp(sequence))
            for candidate in candidates:
                overlap, identity = best_overlap(
                    contig, candidate, min_ovl, min_identity
                )
                if (overlap, identity) > (best_k, best_identity):
                    best_k, best_identity = overlap, identity
                    best_i, best_oriented, best_side = i, candidate, "right"

                overlap, identity = best_overlap(
                    candidate, contig, min_ovl, min_identity
                )
                if (overlap, identity) > (best_k, best_identity):
                    best_k, best_identity = overlap, identity
                    best_i, best_oriented, best_side = i, candidate, "left"

        if best_i is None or best_k < min_ovl:
            break
        original = seqs.pop(best_i)
        used.append(original)
        if best_side == "right":
            contig = contig + best_oriented[best_k:]
        else:
            contig = best_oriented + contig[best_k:]

    return contig, seqs, used


def assemble_all(
    seqs: List[str],
    min_ovl: int,
    min_identity: float,
    consider_revcomp: bool,
    max_contigs: int = 1000,
):
    contigs = []
    remaining = seqs[:]
    for _ in range(max_contigs):
        if not remaining:
            break
        contig, remaining, used = greedy_assemble(
            remaining,
            min_ovl,
            min_identity,
            consider_revcomp=consider_revcomp,
        )
        contigs.append(contig)
        if len(used) <= 1:
            break
    return contigs, remaining


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Greedy reference-free overlap assembly from FASTA or FASTQ."
    )
    inputs = parser.add_mutually_exclusive_group(required=True)
    inputs.add_argument("--products-dir", help="Directory containing FASTA and/or FASTQ files.")
    inputs.add_argument("--fasta", help="Single FASTA input file.")
    inputs.add_argument(
        "--fastq",
        help="Single .fastq or .fq file; relative paths are resolved from the current/root directory.",
    )
    parser.add_argument("--min-ovl", type=int, default=125, help="Minimum required overlap (bp).")
    parser.add_argument(
        "--min-identity",
        type=float,
        default=1.0,
        help="Minimum overlap identity from 0 to 1 (default: 1.0, exact match).",
    )
    parser.add_argument("--no-revcomp", action="store_true", help="Disable reverse-complement matching.")
    parser.add_argument("--out", default="products/denovo_greedy.fasta", help="Output contig FASTA.")
    parser.add_argument(
        "--prefix",
        default="product_",
        help="Filename prefix in directory mode; use --prefix '' to accept every sequence file.",
    )
    args = parser.parse_args()

    if args.min_ovl < 1:
        parser.error("--min-ovl must be at least 1")
    if not 0.0 <= args.min_identity <= 1.0:
        parser.error("--min-identity must be between 0 and 1")

    try:
        if args.products_dir:
            records = load_sequences_from_dir(Path(args.products_dir), prefix=args.prefix)
        else:
            # Path preserves absolute root paths and resolves relative files from cwd.
            records = read_sequences(Path(args.fasta or args.fastq).expanduser())
    except (OSError, ValueError) as error:
        parser.error(str(error))

    contigs, leftovers = assemble_all(
        [sequence for _, sequence in records],
        min_ovl=args.min_ovl,
        min_identity=args.min_identity,
        consider_revcomp=not args.no_revcomp,
    )

    output = Path(args.out)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w") as f:
        for index, contig in enumerate(contigs, 1):
            f.write(
                f">contig_{index}_minovl{args.min_ovl}_"
                f"minid{args.min_identity:g}_len{len(contig)}\n"
            )
            for i in range(0, len(contig), 80):
                f.write(contig[i:i + 80] + "\n")

    print(f"[OK] Wrote {len(contigs)} contigs to {output}")
    print(f"[INFO] Leftover sequences not merged: {len(leftovers)}")
    if leftovers:
        lengths = sorted((len(sequence) for sequence in leftovers), reverse=True)[:10]
        print(f"[INFO] Top leftover lengths: {lengths}")


if __name__ == "__main__":
    main()
