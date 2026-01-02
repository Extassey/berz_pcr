#!/usr/bin/env bash
set -euo pipefail

# Pairwise similarity matrix for FASTA files in the current directory using minimap2.
# Similarity = 100 * (sum(matches) / sum(aligned_block_len)) over all primary alignments in PAF.

threads="${THREADS:-4}"
preset="${PRESET:-asm20}"   # good default for assemblies; change to asm10/asm20 or sr for short reads
out_pairs="${OUT_PAIRS:-pairwise_similarity.tsv}"
out_matrix="${OUT_MATRIX:-similarity_matrix.tsv}"

# Collect FASTA files
mapfile -t files < <(ls -1 *.fa *.fasta *.fna 2>/dev/null || true)

if [[ "${#files[@]}" -lt 2 ]]; then
  echo "Need at least 2 FASTA files (*.fa/*.fasta/*.fna) in the working directory." >&2
  exit 1
fi

# Check dependencies
command -v minimap2 >/dev/null 2>&1 || { echo "minimap2 not found in PATH" >&2; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "awk not found" >&2; exit 1; }

# Long-form pairwise table
echo -e "query\treference\tpct_identity\tmatches\taligned_bp" > "$out_pairs"

# Compute similarity for each ordered pair (i, j)
# PAF fields used:
#  - $10 = number of matching bases in the alignment
#  - $11 = alignment block length
# We sum over primary alignments (excluding secondary/supplementary) to reduce double counting.
for q in "${files[@]}"; do
  for r in "${files[@]}"; do
    # minimap2 outputs to stdout; we parse PAF with awk
    minimap2 -t "$threads" -x "$preset" --secondary=no "$r" "$q" 2>/dev/null \
      | awk -v q="$q" -v r="$r" '
          BEGIN { m=0; a=0; }
          # PAF has at least 12 fields; matches is field 10, aln_block_len is field 11
          NF>=12 { m += $10; a += $11; }
          END {
            if (a==0) {
              printf "%s\t%s\tNA\t0\t0\n", q, r;
            } else {
              printf "%s\t%s\t%.6f\t%d\t%d\n", q, r, (100.0*m/a), m, a;
            }
          }' >> "$out_pairs"
  done
done

# Build similarity matrix from the long-form TSV
# Header row: blank then filenames.
{
  printf "\t"
  printf "%s" "${files[0]}"
  for ((i=1; i<${#files[@]}; i++)); do
    printf "\t%s" "${files[i]}"
  done
  printf "\n"
} > "$out_matrix"

# Fill each row by querying the long-form file
for q in "${files[@]}"; do
  printf "%s" "$q" >> "$out_matrix"
  for r in "${files[@]}"; do
    val=$(awk -v q="$q" -v r="$r" 'BEGIN{FS="\t"} $1==q && $2==r {print $3; exit}' "$out_pairs")
    printf "\t%s" "${val:-NA}" >> "$out_matrix"
  done
  printf "\n" >> "$out_matrix"
done

echo "Wrote:"
echo "  - $out_pairs"
echo "  - $out_matrix"
