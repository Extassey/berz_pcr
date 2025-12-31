import streamlit as st
from pathlib import Path
import subprocess
import sys
import pandas as pd
import shutil



st.set_page_config(page_title="Berz PCR Pipeline", layout="wide")
st.title("Berz PCR Pipeline: Shear → In-silico PCR → Serotype Report")

ROOT = Path(".").resolve()

def run_cmd(cmd, cwd=None):
    """Run a command and stream output into the UI."""
    st.code(" ".join(cmd))
    p = subprocess.Popen(cmd, cwd=cwd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    output_lines = []
    for line in p.stdout:
        output_lines.append(line.rstrip("\n"))
        st.text(line.rstrip("\n"))
    p.wait()
    if p.returncode != 0:
        raise RuntimeError(f"Command failed with exit code {p.returncode}")
    return "\n".join(output_lines)

def reset_dir(path: Path, label: str):
    if path.exists():
        st.warning(f"{label} exists — deleting: {path}")
        shutil.rmtree(path)
    path.mkdir(parents=True, exist_ok=True)


# -----------------------
# INPUTS
# -----------------------
st.subheader("1) Inputs")

col1, col2 = st.columns(2)

with col1:
    resources_dir = st.text_input("Input cps loci folder (FASTA)", value="resources")
    primers_bed = st.text_input("Primer BED file", value="primal_resources/primer.bed")
    refs_fasta = st.text_input("cps reference FASTA", value="cps_refs.fasta")

with col2:
    holocene_dir = st.text_input("Sheared reads output folder", value="holocene")
    products_dir = st.text_input("PCR products folder", value="products")
    report_script = st.text_input("Report script", value="serotype_from_products.py")

st.subheader("2) Shearing settings (main.py)")
c1, c2, c3, c4 = st.columns(4)
with c1:
    chunk_size = st.number_input("Chunk size (bp)", min_value=50, max_value=5000, value=1000, step=50)
with c2:
    target_reads = st.number_input("Target reads", min_value=1, max_value=100000, value=1000, step=50)
with c3:
    seed = st.text_input("Random seed (blank = random each run)", value="")
with c4:
    avoid_dups = st.checkbox("Avoid duplicate windows", value=True)

st.subheader("3) Run controls")
run_shear = st.checkbox("Run shear (main.py)", value=True)
run_pcr = st.checkbox("Run in-silico PCR (insilico.py)", value=True)
run_report = st.checkbox("Run report (serotype_from_products.py)", value=True)

st.divider()

# -----------------------
# RUN
# -----------------------
if st.button("▶ Run Pipeline", type="primary"):
    try:
        # Validate paths (basic)
        resources_path = ROOT / resources_dir
        primers_path = ROOT / primers_bed
        refs_path = ROOT / refs_fasta

        if run_shear and not resources_path.exists():
            st.error(f"Input folder not found: {resources_path}")
            st.stop()
        if run_pcr and not primers_path.exists():
            st.error(f"Primer BED not found: {primers_path}")
            st.stop()
        if run_report and not refs_path.exists():
            st.error(f"Reference FASTA not found: {refs_path}")
            st.stop()

        st.success("Starting…")

        # 1) SHEAR: call main.py with args (we’ll add args support next if needed)
        # If your current main.py doesn’t accept args yet, easiest is:
        # - keep your current main.py constants
        # - OR (recommended) I can modify main.py to accept CLI args.
        if run_shear:
            st.subheader("Shear step (main.py)")

            holocene_path = ROOT / holocene_dir
            reset_dir(holocene_path, "Sheared reads output folder")

            cmd = [
                sys.executable, "main.py",
                "--input", str(resources_path),
                "--output", str(holocene_path),
                "--chunk", str(int(chunk_size)),
                "--target", str(int(target_reads)),
                "--seed", seed.strip() if seed.strip() else "none",
                "--avoid-dups", "1" if avoid_dups else "0",
            ]
            run_cmd(cmd, cwd=str(ROOT))


        # 2) PCR
        if run_pcr:
            st.subheader("In-silico PCR (insilico.py)")

            products_path = ROOT / products_dir
            reset_dir(products_path, "PCR products folder")

            cmd = [
                sys.executable, "insilico.py",
                "--templates", str(ROOT / holocene_dir),
                "--bed", str(primers_path),
                "--out", str(products_path),
            ]
            run_cmd(cmd, cwd=str(ROOT))

        # 3) REPORT
        if run_report:
            st.subheader("Serotype report")
            cmd = [
                sys.executable, report_script,
            ]
            run_cmd(cmd, cwd=str(ROOT))

            ranked = ROOT / products_dir / "serotype_ranked.tsv"
            if ranked.exists():
                df = pd.read_csv(ranked, sep="\t")
                st.subheader("Top serotype candidates")
                st.dataframe(df.head(20), use_container_width=True)
            else:
                st.warning("serotype_ranked.tsv not found (report may have failed or output path differs).")

        st.success("Pipeline finished.")

    except Exception as e:
        st.error(str(e))
