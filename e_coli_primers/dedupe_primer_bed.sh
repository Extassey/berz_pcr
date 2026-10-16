#!/usr/bin/env bash
set -euo pipefail

IN="${1:-master_primers_only_with_header.bed}"
OUT="${2:-master_primers_deduped.bed}"

awk 'BEGIN{FS=OFS="\t"}
     /^#/ {
         # Keep only the first copy of each header/comment line
         if (!seen_header[$0]++) print
         next
     }

     NF >= 7 && ($6 == "+" || $6 == "-") {
         seq = toupper($7)

         # Keep only first row for each unique primer sequence
         if (!seen_seq[seq]++) print

         next
     }
' "$IN" > "$OUT"

echo "Input:  $IN"
echo "Output: $OUT"
echo
echo "Line counts:"
wc -l "$IN" "$OUT"
