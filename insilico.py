#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse
import re

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}

_rc_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(seq: str) -> str:
    return seq.translate(_rc_map)[::-1]

def read_fasta(path: Path):
    header = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq).upper()
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq).upper()

def write_single_fasta(path: Path, header: str, sequence: str):
    with path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")

def find_all(seq: str, sub: str):
    starts = []
    i = seq.find(sub)
    while i != -1:
        starts.append(i)
        i = seq.find(sub, i + 1)
    return starts

bed_line_re = re.compile(r"^\s*(?!#)(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+([+-])\s+([ACGTNacgtn]+)")

def load_primers(bed_path: Path):
    primers_by_amp = defaultdict(lambda: {"LEFT": [], "RIGHT": []})
    with bed_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = bed_line_re.match(line)
            if not m:
                continue
            chrom, start, end, name, amp_s, strand, primer_seq = m.groups()
            amp = int(amp_s)
            primer_seq = primer_seq.upper()

            if "_LEFT" in name:
                side = "LEFT"
            elif "_RIGHT" in name:
                side = "RIGHT"
            else:
                side = "LEFT" if strand == "+" else "RIGHT"

            search_seq = primer_seq if strand == "+" else revcomp(primer_seq)

            primers_by_amp[amp][side].append({
                "name": name,
                "strand": strand,
                "primer_seq": primer_seq,
                "search_seq": search_seq,
                "len": len(primer_seq),
            })
    return primers_by_amp

def best_product_for_amplicon(template_seq: str, left_primers, right_primers, min_len: int, max_len: int):
    left_hits = []
    right_hits = []

    for p in left_primers:
        for pos in find_all(template_seq, p["search_seq"]):
            left_hits.append((pos, p))

    for p in right_primers:
        for pos in find_all(template_seq, p["search_seq"]):
            right_hits.append((pos, p))

    if not left_hits or not right_hits:
        return None, len(left_hits), len(right_hits)

    best = None  # (product_len, lpos, rpos, lp, rp, product_seq)
    for lpos, lp in left_hits:
        for rpos, rp in right_hits:
            if rpos <= lpos:
                continue
            product_end = rpos + rp["len"]
            product_len = product_end - lpos
            if product_len < min_len or product_len > max_len:
                continue
            product_seq = template_seq[lpos:product_end]
            cand = (product_len, lpos, rpos, lp, rp, product_seq)
            if best is None or cand[0] < best[0]:
                best = cand

    return best, len(left_hits), len(right_hits)

def run(templates_dir: Path, bed_path: Path, out_dir: Path,
        min_product_len: int = 1, max_product_len: int = 20000):
    out_dir.mkdir(exist_ok=True)

    primers_by_amp = load_primers(bed_path)
    amplicons = sorted(primers_by_amp.keys())

    template_files = sorted(
        [p for p in templates_dir.iterdir() if p.suffix.lower() in FASTA_EXTENSIONS],
        key=lambda p: p.name.lower()
    )

    report_path = out_dir / "amplification_report.tsv"
    product_index = 1

    with report_path.open("w") as rep:
        rep.write("\t".join([
            "product_id", "template_file", "template_record", "amplicon",
            "status", "product_length",
            "left_primer_name", "right_primer_name",
            "left_pos_0based", "right_pos_0based",
            "left_hits", "right_hits",
        ]) + "\n")

        for tf in template_files:
            for record_header, seq in read_fasta(tf):
                for amp in amplicons:
                    L = primers_by_amp[amp]["LEFT"]
                    R = primers_by_amp[amp]["RIGHT"]

                    best, left_hits, right_hits = best_product_for_amplicon(
                        seq, L, R, min_product_len, max_product_len
                    )

                    if best is None:
                        continue

                    product_len, lpos, rpos, lp, rp, product_seq = best
                    pid = f"product_{product_index:06d}"
                    write_single_fasta(out_dir / f"{pid}.fasta", pid, product_seq)

                    rep.write("\t".join([
                        pid, tf.name, record_header, str(amp),
                        "OK", str(product_len),
                        lp["name"], rp["name"],
                        str(lpos), str(rpos),
                        str(left_hits), str(right_hits),
                    ]) + "\n")

                    product_index += 1

    print(f"Done. Products: {product_index-1}")
    print(f"Products folder: {out_dir.resolve()}")
    print(f"Report: {report_path.resolve()}")

def parse_args():
    ap = argparse.ArgumentParser(description="In-silico PCR using ARTIC-style primer.bed against template FASTAs.")
    ap.add_argument("--templates", default="holocene", help="Folder containing template FASTA files.")
    ap.add_argument("--bed", default="primal_resources/primer.bed", help="Primer BED file.")
    ap.add_argument("--out", default="products", help="Output folder for products and report.")
    ap.add_argument("--min-len", type=int, default=450, help="Minimum allowed product length (bp).")
    ap.add_argument("--max-len", type=int, default=650, help="Maximum allowed product length (bp).")
    
    return ap.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run(Path(args.templates), Path(args.bed), Path(args.out), args.min_len, args.max_len)
