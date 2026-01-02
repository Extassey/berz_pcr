#!/usr/bin/env python3
from __future__ import annotations

import argparse
import subprocess
from pathlib import Path
import sys
import tempfile

"""
Here is how to run from terminal


python3 denovo.py \
  --products-dir products \
  --ref-fasta cps_refs.fasta \
  --contig "9N/9L" \
  --out products/consensus_9N_9L.fasta \
  --min-depth 3


and if you want no masking

python3 denovo.py --products-dir products --ref-fasta cps_refs.fasta --contig "9N/9L" --out products/consensus_no_mask.fasta --min-depth 0

"""



FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}


def run(cmd: list[str], *, cwd: Path | None = None, stdout_path: Path | None = None) -> None:
    """Run a command, streaming stderr to terminal. Optionally write stdout to a file."""
    print("[cmd]", " ".join(cmd), flush=True)
    if stdout_path is None:
        subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)
    else:
        with stdout_path.open("w") as out:
            subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True, stdout=out)


def require_exists(p: Path, what: str) -> None:
    if not p.exists():
        raise SystemExit(f"{what} not found: {p}")


def detect_products_fasta(products_dir: Path) -> Path:
    """
    Prefer products/all_products.fasta if present (your pipeline often creates this),
    otherwise concatenate individual product_*.fasta into a temp file.
    """
    all_fa = products_dir / "all_products.fasta"
    if all_fa.exists():
        return all_fa

    # Concatenate product FASTAs
    fasta_files = sorted(
        [p for p in products_dir.iterdir() if p.suffix.lower() in FASTA_EXTENSIONS],
        key=lambda x: x.name
    )
    fasta_files = [p for p in fasta_files if p.name.startswith("product_")]
    if not fasta_files:
        raise SystemExit(f"No product_*.fasta files found in {products_dir}")

    tmp = products_dir / "_tmp_all_products.fasta"
    with tmp.open("w") as out:
        for fp in fasta_files:
            out.write(fp.read_text())
            if not fp.read_text().endswith("\n"):
                out.write("\n")
    return tmp


def extract_contig_fasta(ref_fasta: Path, contig: str, out_fasta: Path) -> None:
    """
    Extract exactly one contig from a multi-fasta reference to avoid ambiguity.
    Uses samtools faidx if available; otherwise does a simple parser fallback.
    """
    # Try samtools faidx first (fast, robust)
    try:
        subprocess.run(["samtools", "faidx", str(ref_fasta)], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        # samtools faidx ref.fa "contig" > out.fa
        run(["samtools", "faidx", str(ref_fasta), contig], stdout_path=out_fasta)
        return
    except Exception:
        pass

    # Fallback parser
    header = None
    seq = []
    found = False
    with ref_fasta.open() as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None and found:
                    break
                header = line[1:].split()[0]
                found = (header == contig)
                seq = []
            else:
                if found:
                    seq.append(line.strip())
    if not found or not seq:
        raise SystemExit(f"Contig '{contig}' not found in reference FASTA: {ref_fasta}")

    with out_fasta.open("w") as out:
        out.write(f">{contig}\n")
        s = "".join(seq).upper()
        for i in range(0, len(s), 80):
            out.write(s[i:i+80] + "\n")


def make_low_depth_mask_bed(bam: Path, contig: str, min_depth: int, out_bed: Path) -> None:
    """
    Create a BED mask where depth < min_depth. bcftools consensus can mask these to N.
    BED is 0-based half-open intervals.
    """
    require_exists(bam, "BAM")
    # samtools depth -aa -r contig bam
    p = subprocess.Popen(
        ["samtools", "depth", "-aa", "-r", contig, str(bam)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    prev_start = None
    prev_end = None

    with out_bed.open("w") as out:
        for line in p.stdout:
            # contig \t pos(1-based) \t depth
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 3:
                continue
            _, pos1, depth = parts
            pos1 = int(pos1)
            depth = int(depth)

            if depth < min_depth:
                start0 = pos1 - 1
                end0 = pos1
                if prev_end == start0:  # extend interval
                    prev_end = end0
                else:
                    if prev_start is not None:
                        out.write(f"{contig}\t{prev_start}\t{prev_end}\n")
                    prev_start, prev_end = start0, end0
            else:
                if prev_start is not None:
                    out.write(f"{contig}\t{prev_start}\t{prev_end}\n")
                    prev_start = prev_end = None

        if prev_start is not None:
            out.write(f"{contig}\t{prev_start}\t{prev_end}\n")

    p.wait()
    if p.returncode != 0:
        err = p.stderr.read()
        raise SystemExit(f"samtools depth failed:\n{err}")


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Build a consensus cps sequence from tiled PCR products (reference-guided consensus)."
    )
    ap.add_argument("--products-dir", default="products", help="Folder containing product FASTAs (and/or all_products.fasta).")
    ap.add_argument("--ref-fasta", required=True, help="Reference FASTA containing cps contigs (multi-fasta OK).")
    ap.add_argument("--contig", required=True, help="Which reference contig to build consensus against (e.g. '9N/9L').")
    ap.add_argument("--out", default="products/consensus_cps.fasta", help="Output consensus FASTA path.")
    ap.add_argument("--min-depth", type=int, default=3, help="Mask positions with depth < this to 'N'. Set 0 to disable masking.")
    ap.add_argument("--threads", type=int, default=4, help="Threads for minimap2/samtools sort where applicable.")
    ap.add_argument("--minimap2", default="minimap2", help="minimap2 binary.")
    ap.add_argument("--samtools", default="samtools", help="samtools binary.")
    ap.add_argument("--bcftools", default="bcftools", help="bcftools binary.")
    args = ap.parse_args()

    products_dir = Path(args.products_dir)
    ref_fasta = Path(args.ref_fasta)
    out_fa = Path(args.out)

    require_exists(products_dir, "Products dir")
    require_exists(ref_fasta, "Reference FASTA")

    # Assemble a single FASTA of products
    products_fa = detect_products_fasta(products_dir)

    with tempfile.TemporaryDirectory() as td:
        td = Path(td)

        # Extract the single reference contig into its own FASTA
        ref_one = td / "ref_one.fasta"
        extract_contig_fasta(ref_fasta, args.contig, ref_one)

        # Align products -> reference
        sam = td / "aln.sam"
        bam = td / "aln.sorted.bam"
        vcf = td / "calls.vcf.gz"
        mask_bed = td / "mask_low_depth.bed"

        # minimap2 -a -x sr --secondary=no ref_one products
        run([args.minimap2, "-a", "-x", "sr", "--secondary=no", str(ref_one), str(products_fa)], stdout_path=sam)

        # samtools view -bS | sort
        # samtools sort -@ threads -o bam sam
        run([args.samtools, "sort", "-@", str(args.threads), "-o", str(bam), str(sam)])
        run([args.samtools, "index", str(bam)])

        # Optional mask
        if args.min_depth and args.min_depth > 0:
            make_low_depth_mask_bed(bam, args.contig, args.min_depth, mask_bed)

        # Call variants
        # bcftools mpileup -Ou -f ref_one -r contig bam | bcftools call -mv -Ov -o vcf
        mp = subprocess.Popen(
            [args.bcftools, "mpileup", "-Ou", "-f", str(ref_one), "-r", args.contig, str(bam)],
            stdout=subprocess.PIPE,
            text=False,
        )
        subprocess.run([args.bcftools, "call", "-mv", "-Oz", "-o", str(vcf)], stdin=mp.stdout, check=True)
        mp.stdout.close()
        mp.wait()

        run([args.bcftools, "index", str(vcf)])


        # Build consensus
        out_fa.parent.mkdir(parents=True, exist_ok=True)
        cmd = [args.bcftools, "consensus", "-f", str(ref_one)]
        if args.min_depth and args.min_depth > 0:
            cmd += ["-m", str(mask_bed)]
        cmd += [str(vcf)]
        run(cmd, stdout_path=out_fa)

    print(f"[OK] Consensus written: {out_fa}")


if __name__ == "__main__":
    try:
        main()
    except subprocess.CalledProcessError as e:
        raise SystemExit(f"Command failed: {e}")
