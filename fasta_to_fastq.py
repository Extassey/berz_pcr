#!/usr/bin/env python3

import argparse
from pathlib import Path


def read_fasta(path: Path):
    """Yield (header, sequence) records from a FASTA file."""
    header = None
    sequence_parts = []

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(sequence_parts)

                header = line[1:].strip()
                sequence_parts = []
            else:
                if header is None:
                    raise ValueError("Sequence found before the first FASTA header.")
                sequence_parts.append(line)

    if header is not None:
        yield header, "".join(sequence_parts)


def fasta_to_fastq(input_path: Path, output_path: Path, quality_char: str = "I"):
    with output_path.open("w", encoding="utf-8") as output:
        for header, sequence in read_fasta(input_path):
            sequence = sequence.upper()

            output.write(f"@{header}\n")
            output.write(f"{sequence}\n")
            output.write("+\n")
            output.write(f"{quality_char * len(sequence)}\n")


def main():
    parser = argparse.ArgumentParser(
        description="Convert FASTA to FASTQ using a fixed quality score."
    )
    parser.add_argument("input_fasta", type=Path, help="Input FASTA file")
    parser.add_argument("output_fastq", type=Path, help="Output FASTQ file")
    parser.add_argument(
        "--quality-char",
        default="I",
        help="FASTQ quality character to assign to every base. Default: I (Phred Q40)",
    )

    args = parser.parse_args()

    if len(args.quality_char) != 1:
        parser.error("--quality-char must be exactly one character.")

    fasta_to_fastq(
        input_path=args.input_fasta,
        output_path=args.output_fastq,
        quality_char=args.quality_char,
    )


if __name__ == "__main__":
    main()