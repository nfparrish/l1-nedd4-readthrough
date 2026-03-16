#!/usr/bin/env python3
"""
Select 20 samples per genotype group (0/0, 0/1, 1/1) for the NEDD4-L1
insertion from PRJNA851328, then fetch their SRR accessions via E-utilities.

Reads:
  /work/nfp8/LCL_1kgp/metadata/genotypes.csv  (sample_id,genotype)

Writes:
  selected_samples.tsv   - sample_id, genotype, SRR, size_gb, libname
  srr_list.txt           - one SRR per line (all 60)
  srr_by_genotype/       - separate SRR lists per genotype for job arrays
"""
import csv, io, json, random, time, urllib.request
from pathlib import Path

OUTDIR = Path("/work/nfp8/LCL_1kgp/metadata")
OUTDIR.mkdir(parents=True, exist_ok=True)

# ── 1. Load genotypes ─────────────────────────────────────────────────────────
gt_to_samples = {"0/0": [], "0/1": [], "1/1": []}
with open(OUTDIR / "genotypes.csv") as f:
    for row in csv.DictReader(f):
        gt = row["genotype"]
        if gt in gt_to_samples:
            gt_to_samples[gt].append(row["sample_id"])

print("Genotype counts: " +
      ", ".join(f"{k}: {len(v)}" for k, v in gt_to_samples.items()))

# ── 2. Fetch SRA run info ────────────────────────────────────────────────────
def fetch_url(url, retries=3):
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                return r.read().decode()
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(5)
            else:
                raise

runinfo_path = OUTDIR / "sra_runinfo_full.csv"
if not runinfo_path.exists():
    print("Fetching SRA run info via E-utilities...", flush=True)
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    search_url = (f"{base}esearch.fcgi?db=sra&term=PRJNA851328[BioProject]"
                  f"&retmax=1000&usehistory=y&retmode=json")
    search_res = json.loads(fetch_url(search_url))
    webenv = search_res["esearchresult"]["webenv"]
    query_key = search_res["esearchresult"]["querykey"]
    count = int(search_res["esearchresult"]["count"])
    print(f"  Found {count} SRA records", flush=True)

    fetch_url2 = (f"{base}efetch.fcgi?db=sra&query_key={query_key}"
                  f"&WebEnv={webenv}&rettype=runinfo&retmode=text&retmax=1000")
    runinfo_text = fetch_url(fetch_url2)
    runinfo_path.write_text(runinfo_text)
    print(f"  Saved {runinfo_path}", flush=True)
else:
    runinfo_text = runinfo_path.read_text()
    print(f"Using cached {runinfo_path}", flush=True)

# ── 3. Parse run info ─────────────────────────────────────────────────────────
reader = csv.DictReader(io.StringIO(runinfo_text))
rows = [r for r in reader if r.get("Run", "").startswith("SRR")]
print(f"Parsed {len(rows)} SRR rows", flush=True)

# Build sample_id -> list of runs; LibraryName: "{sample_id}_{batch}_{rep}"
sample_to_runs = {}
for r in rows:
    libname = r.get("LibraryName", "")
    parts = libname.split("_")
    sid = parts[0]
    rep = parts[-1] if parts[-1].startswith("rep") else "rep1"
    srr = r["Run"]
    size_mb = float(r.get("size_MB", 0) or 0)
    sample_to_runs.setdefault(sid, []).append(
        {"srr": srr, "rep": rep, "size_mb": size_mb, "libname": libname}
    )

print(f"Unique sample IDs in SRA: {len(sample_to_runs)}", flush=True)

# ── 4. Select 20 per genotype ────────────────────────────────────────────────
def best_run(runs):
    rep1 = [r for r in runs if r["rep"] == "rep1"]
    pool = rep1 if rep1 else runs
    return min(pool, key=lambda r: r["size_mb"])

random.seed(42)
selected = {}
for gt, samples in gt_to_samples.items():
    with_srr = [s for s in samples if s in sample_to_runs]
    missing  = [s for s in samples if s not in sample_to_runs]
    if missing:
        print(f"  WARNING: {len(missing)} {gt} samples not in SRA: "
              f"{missing[:5]}", flush=True)
    n = min(20, len(with_srr))
    chosen = random.sample(with_srr, n)
    selected[gt] = [(s, best_run(sample_to_runs[s])) for s in chosen]
    print(f"  {gt}: selected {n}/20 from {len(with_srr)} available", flush=True)

# ── 5. Write outputs ──────────────────────────────────────────────────────────
tsv_path = OUTDIR / "selected_samples.tsv"
srr_all_path = OUTDIR / "srr_list.txt"
(OUTDIR / "srr_by_genotype").mkdir(exist_ok=True)

with open(tsv_path, "w") as tsv, open(srr_all_path, "w") as all_f:
    tsv.write("sample_id\tgenotype\tSRR\tsize_gb\tlibname\n")
    for gt, items in selected.items():
        gt_label = gt.replace("/", "")
        gt_path = OUTDIR / "srr_by_genotype" / f"srr_{gt_label}.txt"
        with open(gt_path, "w") as gf:
            for sid, run in items:
                size_gb = run["size_mb"] / 1024
                line = (f"{sid}\t{gt}\t{run['srr']}\t"
                        f"{size_gb:.2f}\t{run['libname']}\n")
                tsv.write(line)
                all_f.write(run["srr"] + "\n")
                gf.write(run["srr"] + "\n")

print(f"\nWrote {tsv_path}")
print(f"Wrote {srr_all_path}")
print(f"Total SRRs: {sum(len(v) for v in selected.values())}")
