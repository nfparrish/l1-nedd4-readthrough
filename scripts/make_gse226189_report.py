#!/usr/bin/env python3

from __future__ import annotations

import base64
import csv
import html
from pathlib import Path


RESULTS_DIR = Path("/work/nfp8/2026_03_02/GSE226189/results")
PERSISTENT_DIR = Path("/hpc/group/parrishlab/L1_NEDD4/results/GSE226189")
FIG_DIR = PERSISTENT_DIR / "figures"
OUT_HTML = PERSISTENT_DIR / "gse226189_report.html"


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def read_single_column(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    with path.open() as handle:
        reader = csv.reader(handle)
        header = next(reader)
        value_idx = 1 if len(header) > 1 else 0
        for row in reader:
            if not row:
                continue
            values[row[0]] = float(row[value_idx])
    return values


def read_tmm_row(path: Path, gene_id: str) -> dict[str, float]:
    with path.open() as handle:
        reader = csv.reader(handle)
        header = next(reader)
        for row in reader:
            if row and row[0] == gene_id:
                return {header[idx]: float(row[idx]) for idx in range(1, len(header))}
    return {}


def format_float(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def table_html(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html.escape(cell)}</th>" for cell in headers)
    body = []
    for row in rows:
        body.append("<tr>" + "".join(f"<td>{html.escape(cell)}</td>" for cell in row) + "</tr>")
    return "<table><thead><tr>" + head + "</tr></thead><tbody>" + "".join(body) + "</tbody></table>"


def embed_png(path: Path) -> str:
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


def main() -> None:
    signal_rows = read_csv_dicts(RESULTS_DIR / "readthrough_signal.csv")
    qc_rows = {row["srr"]: row for row in read_csv_dicts(RESULTS_DIR / "star_alignment_qc.csv")}
    nedd4_counts = read_single_column(RESULTS_DIR / "nedd4_fc_counts.csv")
    nedd4_tmm = read_tmm_row(RESULTS_DIR / "tmm_cpm_matrix.csv", "ENSG00000069869.17")
    paper_note = (FIG_DIR / "GSE226189_paper_style_note.txt").read_text().strip()
    summary_text = (FIG_DIR / "GSE226189_summary_stats.txt").read_text().strip()

    for row in signal_rows:
        row["rt_cpm_f"] = float(row["rt_cpm"])
        row["same_window_minus_cpm_f"] = float(row["same_window_minus_cpm"])
        row["rt_bg_log2fold_f"] = float(row["rt_bg_log2fold"]) if row.get("rt_bg_log2fold") else None
        row["unique_mapped_f"] = float(qc_rows[row["srr"]]["pct_unique_mapped"])
        row["nedd4_count_f"] = nedd4_counts.get(row["srr"], 0.0)
        row["nedd4_tmm_f"] = nedd4_tmm.get(row["srr"], 0.0)

    sample_count = len(signal_rows)
    nonzero = sum(1 for row in signal_rows if row["rt_cpm_f"] > 0)
    median_rt = sorted(row["rt_cpm_f"] for row in signal_rows)[sample_count // 2]
    top_rt = max(signal_rows, key=lambda row: row["rt_cpm_f"])
    top_nedd4 = max(signal_rows, key=lambda row: row["nedd4_tmm_f"])

    top_rt_rows = []
    for row in sorted(signal_rows, key=lambda item: item["rt_cpm_f"], reverse=True)[:15]:
        top_rt_rows.append([
            row["srr"],
            format_float(row["rt_cpm_f"], 6),
            format_float(row["same_window_minus_cpm_f"], 6),
            format_float(row["rt_bg_log2fold_f"], 4) if row["rt_bg_log2fold_f"] is not None else "NA",
            format_float(row["unique_mapped_f"], 2),
        ])

    low_rt_rows = []
    for row in sorted(signal_rows, key=lambda item: item["rt_cpm_f"])[:12]:
        low_rt_rows.append([
            row["srr"],
            format_float(row["rt_cpm_f"], 6),
            format_float(row["same_window_minus_cpm_f"], 6),
            format_float(row["unique_mapped_f"], 2),
        ])

    nedd4_rows = []
    for row in sorted(signal_rows, key=lambda item: item["nedd4_tmm_f"], reverse=True)[:15]:
        nedd4_rows.append([
            row["srr"],
            str(int(row["nedd4_count_f"])),
            format_float(row["nedd4_tmm_f"], 2),
            format_float(row["rt_cpm_f"], 6),
        ])

    cards = [
        ("Samples", str(sample_count)),
        ("Median RT CPM", format_float(median_rt, 6)),
        ("Non-zero RT signal", f"{nonzero}/{sample_count}"),
        ("Top RT sample", top_rt["srr"]),
        ("Top NEDD4 sample", top_nedd4["srr"]),
        ("Paper-style status", "Unavailable for unstranded cohort"),
    ]

    image_specs = [
        ("Readthrough CPM Distribution", FIG_DIR / "GSE226189_rt_cpm_strip.png", "Cohort-wide strip plot of readthrough CPM values."),
        ("Aggregated Coverage Profile", FIG_DIR / "GSE226189_rt_coverage_profile.png", "Mean and interquartile depth across the downstream 1 kb window."),
        ("Ranked RT CPM", FIG_DIR / "GSE226189_rt_cpm_ranked.png", "Samples ranked by readthrough CPM."),
    ]
    image_cards = []
    for title, path, caption in image_specs:
        image_cards.append(
            f'<div class="chart-card"><h3>{html.escape(title)}</h3><img src="{embed_png(path)}" alt="{html.escape(title)}"><p class="small">{html.escape(caption)}</p></div>'
        )

    summary_lines = [line for line in summary_text.splitlines() if line and not line.startswith("=")]
    summary_callouts = "".join(f'<div class="callout">{html.escape(line)}</div>' for line in summary_lines[:8])

    html_text = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>GSE226189 Cohort Summary</title>
  <style>
    :root {{
      --bg: #eef1e7;
      --panel: rgba(252, 253, 248, 0.93);
      --ink: #182018;
      --muted: #596556;
      --line: rgba(46, 71, 41, 0.15);
      --accent: #9a4d3c;
      --accent-2: #3b7f6a;
      --accent-3: #5878b2;
      --shadow: 0 18px 50px rgba(31, 47, 26, 0.10);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: Georgia, \"Times New Roman\", serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(154, 77, 60, 0.10), transparent 28%),
        radial-gradient(circle at top right, rgba(59, 127, 106, 0.12), transparent 26%),
        linear-gradient(180deg, #f5f7ef 0%, var(--bg) 44%, #e7ebdf 100%);
    }}
    .wrap {{ max-width: 1320px; margin: 0 auto; padding: 32px 20px 48px; }}
    .hero {{ padding: 28px; border: 1px solid var(--line); border-radius: 24px; background: linear-gradient(135deg, rgba(255,255,255,0.88), rgba(243, 246, 236, 0.94)); box-shadow: var(--shadow); }}
    h1, h2, h3 {{ margin: 0 0 12px; }}
    h1 {{ font-size: clamp(2rem, 4.8vw, 3.8rem); line-height: 0.96; letter-spacing: -0.04em; }}
    h2 {{ font-size: 1.5rem; }}
    p {{ margin: 0; line-height: 1.55; }}
    .sub {{ color: var(--muted); max-width: 80ch; margin-top: 14px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(210px, 1fr)); gap: 18px; margin-top: 22px; }}
    .card, .chart-card {{ background: var(--panel); border: 1px solid var(--line); border-radius: 20px; padding: 18px; box-shadow: var(--shadow); }}
    .stat-label {{ color: var(--muted); font-size: 0.95rem; margin-bottom: 6px; }}
    .stat-value {{ font-size: 1.95rem; font-weight: 700; letter-spacing: -0.04em; }}
    .section {{ margin-top: 26px; }}
    .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 18px; margin-top: 16px; }}
    img {{ width: 100%; height: auto; display: block; border-radius: 14px; border: 1px solid var(--line); background: white; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.94rem; }}
    th, td {{ padding: 10px 12px; border-bottom: 1px solid var(--line); text-align: left; }}
    th {{ color: var(--muted); font-weight: 600; position: sticky; top: 0; background: rgba(252,253,248,0.98); }}
    tbody tr:hover {{ background: rgba(255,255,255,0.5); }}
    .table-wrap {{ max-height: 520px; overflow: auto; }}
    .callouts {{ display: grid; gap: 10px; margin-top: 16px; }}
    .callout {{ padding: 14px 16px; border-left: 4px solid var(--accent); background: rgba(255,255,255,0.58); border-radius: 12px; }}
    .small {{ color: var(--muted); font-size: 0.92rem; }}
    pre {{ white-space: pre-wrap; font-family: Georgia, \"Times New Roman\", serif; font-size: 0.95rem; margin: 0; }}
    @media (max-width: 720px) {{
      .wrap {{ padding: 18px 14px 28px; }}
      .hero, .card, .chart-card {{ border-radius: 16px; }}
    }}
  </style>
</head>
<body>
  <div class=\"wrap\">
    <section class=\"hero\">
      <h1>GSE226189 Cohort Summary</h1>
      <p class=\"sub\">This report summarizes the 82-sample fibroblast cohort using cohort-wide L1-NEDD4 readthrough metrics, same-window opposite-strand background signal, STAR alignment quality, and NEDD4 abundance context. Unlike the MCF7 control set, this cohort is configured as unstranded, so paper-style stranded mate1 outputs are intentionally unavailable.</p>
      <div class=\"grid\">{"".join(f'<div class="card"><div class="stat-label">{html.escape(label)}</div><div class="stat-value">{html.escape(value)}</div></div>' for label, value in cards)}</div>
      <div class=\"callouts\">{summary_callouts}</div>
    </section>

    <section class=\"section\">
      <h2>Cohort Figures</h2>
      <div class=\"plot-grid\">{"".join(image_cards)}</div>
    </section>

    <section class=\"section\">
      <h2>Paper-Style Status</h2>
      <div class=\"card\"><pre>{html.escape(paper_note)}</pre></div>
    </section>

    <section class=\"section\">
      <h2>Highest Readthrough Samples</h2>
      <div class=\"card table-wrap\">{table_html(['Sample', 'RT CPM', 'Same-window minus CPM', 'log2(RT/bg)', 'Unique mapped %'], top_rt_rows)}</div>
    </section>

    <section class=\"section\">
      <h2>Lowest Readthrough Samples</h2>
      <div class=\"card table-wrap\">{table_html(['Sample', 'RT CPM', 'Same-window minus CPM', 'Unique mapped %'], low_rt_rows)}</div>
    </section>

    <section class=\"section\">
      <h2>NEDD4 Context</h2>
      <p class=\"small\">NEDD4 remains broadly expressed across the cohort. The table below ranks samples by NEDD4 TMM CPM from the featureCounts-derived normalization matrix.</p>
      <div class=\"card table-wrap\">{table_html(['Sample', 'NEDD4 featureCounts', 'NEDD4 TMM CPM', 'RT CPM'], nedd4_rows)}</div>
    </section>
  </div>
</body>
</html>
"""

    OUT_HTML.write_text(html_text)
    print(f"wrote {OUT_HTML}")


if __name__ == "__main__":
    main()