import csv
import primer3

bed = "/home/darklightheart/Berz_PCR/final_assay_75_75/primer.bed"

rows = []

with open(bed) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip().split("\t")
        serotype = fields[0]
        name = fields[3]
        pool = fields[4]
        seq = fields[6].upper()

        tm = primer3.calc_tm(seq)

        gc = 100 * (seq.count("G") + seq.count("C")) / len(seq)

        rows.append({
            "serotype": serotype,
            "name": name,
            "pool": pool,
            "length": len(seq),
            "gc_percent": gc,
            "tm": tm,
            "sequence": seq,
        })

rows.sort(key=lambda x: x["tm"])

print("\nLOWEST 10 Tm")
for r in rows[:10]:
    print(
        f'{r["tm"]:6.2f} C  '
        f'{r["name"]:35s} '
        f'len={r["length"]:2d}  GC={r["gc_percent"]:5.1f}%'
    )

print("\nHIGHEST 10 Tm")
for r in rows[-10:][::-1]:
    print(
        f'{r["tm"]:6.2f} C  '
        f'{r["name"]:35s} '
        f'len={r["length"]:2d}  GC={r["gc_percent"]:5.1f}%'
    )

tms = [r["tm"] for r in rows]

print("\nSUMMARY")
print(f"Primers: {len(rows)}")
print(f"Minimum Tm: {min(tms):.2f} C")
print(f"Maximum Tm: {max(tms):.2f} C")
print(f"Delta Tm:   {max(tms)-min(tms):.2f} C")
print(f"Mean Tm:    {sum(tms)/len(tms):.2f} C")

with open("primer_tm_report.csv", "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)

print("\nWrote primer_tm_report.csv")