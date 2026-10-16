from pathlib import Path
import csv

# Directory containing your per-serotype FASTAs like 8.fasta, 12A.fasta, etc.
FASTA_DIR = Path("split_fastas")

def read_fasta_sequence(path):
    """Read a (single-entry) FASTA and return the sequence as one string."""
    seq_lines = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq_lines.append(line)
    return "".join(seq_lines).upper()

def percent_identity(seq1, seq2):
    """
    Compute percent identity between two sequences.
    Uses the length of the shorter sequence,
    ignoring positions where either base is a gap ('-').
    """
    L = min(len(seq1), len(seq2))
    if L == 0:
        return 0.0

    matches = 0
    comparable_positions = 0

    for a, b in zip(seq1[:L], seq2[:L]):
        if a == "-" or b == "-":
            continue
        comparable_positions += 1
        if a == b:
            matches += 1

    if comparable_positions == 0:
        return 0.0

    return (matches / comparable_positions) * 100.0


def main():
    # Grab all FASTA files in the directory
    fasta_paths = sorted(FASTA_DIR.glob("*.fasta"))
    if not fasta_paths:
        print(f"No FASTA files found in {FASTA_DIR}")
        return

    # Load names and sequences
    names = []
    seqs = []
    for path in fasta_paths:
        name = path.stem  # e.g. "8", "12A", "15A"
        seq = read_fasta_sequence(path)
        names.append(name)
        seqs.append(seq)

    n = len(names)
    print(f"Loaded {n} sequences from {FASTA_DIR}")

    # Compute pairwise PID matrix
    pid_matrix = [[0.0] * n for _ in range(n)]

    for i in range(n):
        for j in range(i + 1, n):
            pid = percent_identity(seqs[i], seqs[j])
            pid_matrix[i][j] = pid
            pid_matrix[j][i] = pid

    # Compute average PID for each sequence (excluding self)
    avg_pid_list = []
    for i in range(n):
        others = [pid_matrix[i][j] for j in range(n) if j != i]
        avg_pid = sum(others) / len(others) if others else 0.0
        avg_pid_list.append((names[i], avg_pid))

    # Sort by descending average PID
    avg_pid_list.sort(key=lambda x: x[1], reverse=True)

    # Write CSV
    output_csv = "average_pid_by_serotype.csv"
    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["serotype", "average_pid"])
        for name, avg_pid in avg_pid_list:
            writer.writerow([name, f"{avg_pid:.4f}"])

    print(f"Average PID table written to {output_csv}")


if __name__ == "__main__":
    main()
