#!/usr/bin/env python3

from __future__ import annotations

import csv
import html
import json
from pathlib import Path


RESULTS_DIR = Path("/work/nfp8/MCF7_ctrl/results")
COVERAGE_DIR = Path("/work/nfp8/MCF7_ctrl/coverage")
COUNTS_DIR = Path("/work/nfp8/MCF7_ctrl/counts")
OUT_HTML = RESULTS_DIR / "mcf7_report.html"


def read_csv_dicts(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle))


def read_matrix(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open() as handle:
        reader = csv.reader(handle)
        header = next(reader)
        rows = []
        for row in reader:
            if not row:
                continue
            rows.append({header[idx]: value for idx, value in enumerate(row)})
        return header, rows


def read_depth(path: Path) -> list[int]:
    values: list[int] = []
    if not path.exists():
        return values
    with path.open() as handle:
        for line in handle:
            parts = line.strip().split()
            if len(parts) == 3:
                values.append(int(parts[2]))
    return values


def featurecounts_summary(sample: str) -> dict[str, int]:
    summary_path = COUNTS_DIR / f"{sample}_featureCounts_PANEL.txt.summary"
    metrics: dict[str, int] = {}
    with summary_path.open() as handle:
        next(handle)
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) == 2:
                metrics[parts[0]] = int(parts[1])
    return metrics


def format_float(value: float, digits: int = 3) -> str:
    return f"{value:.{digits}f}"


def table_html(headers: list[str], rows: list[list[str]]) -> str:
    head = "".join(f"<th>{html.escape(cell)}</th>" for cell in headers)
    body_rows = []
    for row in rows:
        body_rows.append("<tr>" + "".join(f"<td>{html.escape(cell)}</td>" for cell in row) + "</tr>")
    return "<table><thead><tr>" + head + "</tr></thead><tbody>" + "".join(body_rows) + "</tbody></table>"


def to_float(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    return float(value) if value not in {"", "None", None} else 0.0


def depth_positions(depth_values: list[int], start: int) -> list[int]:
    return list(range(start, start + len(depth_values)))


def main() -> None:
    signal_rows = read_csv_dicts(RESULTS_DIR / "readthrough_signal.csv")
    paper_rows = {row["srr"]: row for row in read_csv_dicts(RESULTS_DIR / "paper_style_signal.csv")}
    qc_rows = {row["srr"]: row for row in read_csv_dicts(RESULTS_DIR / "star_alignment_qc.csv")}
    _, tmm_rows = read_matrix(RESULTS_DIR / "tmm_cpm_matrix.csv")
    _, fc_rows = read_matrix(RESULTS_DIR / "fc_count_matrix.csv")

    samples = [row["srr"] for row in signal_rows]
    signal_by_sample = {row["srr"]: row for row in signal_rows}

    depth_data: dict[str, dict[str, list[int]]] = {}
    paper_depth_data: dict[str, dict[str, list[int]]] = {}
    for sample in samples:
        depth_data[sample] = {
            "plus": read_depth(COVERAGE_DIR / f"{sample}_readthrough.depth"),
            "same_window_minus": read_depth(COVERAGE_DIR / f"{sample}_readthrough_minus.depth"),
        }
        paper_depth_data[sample] = {
            "downstream_sense": read_depth(COVERAGE_DIR / f"{sample}_paper_downstream_sense.depth"),
            "upstream_antisense": read_depth(COVERAGE_DIR / f"{sample}_paper_upstream_antisense.depth"),
        }

    downstream_start = int(next(iter(paper_rows.values()))["paper_downstream_window"].split(":", 1)[1].split("-", 1)[0])
    upstream_start = int(next(iter(paper_rows.values()))["paper_upstream_window"].split(":", 1)[1].split("-", 1)[0])

    signal_summary_rows: list[list[str]] = []
    paper_summary_rows: list[list[str]] = []
    assignment_rows: list[list[str]] = []
    for sample in samples:
        row = signal_by_sample[sample]
        paper = paper_rows[sample]
        qc = qc_rows[sample]
        signal_summary_rows.append([
            sample,
            format_float(to_float(row, "rt_mean_depth"), 4),
            format_float(to_float(row, "rt_cpm"), 6),
            format_float(to_float(row, "same_window_minus_mean_depth"), 4),
            format_float(to_float(row, "same_window_minus_cpm"), 6),
            format_float(to_float(row, "bg_mean_depth"), 4),
            format_float(to_float(row, "rt_bg_log2fold"), 4),
            format_float(to_float(qc, "pct_unique_mapped"), 2),
        ])
        paper_summary_rows.append([
            sample,
            paper.get("paper_downstream_window", ""),
            paper.get("paper_upstream_window", ""),
            format_float(to_float(paper, "paper_downstream_sense_mean_depth"), 4),
            format_float(to_float(paper, "paper_downstream_sense_cpm"), 6),
            format_float(to_float(paper, "paper_upstream_antisense_mean_depth"), 4),
            format_float(to_float(paper, "paper_upstream_antisense_cpm"), 6),
        ])

        fc = featurecounts_summary(sample)
        assigned = fc.get("Assigned", 0)
        no_features = fc.get("Unassigned_NoFeatures", 0)
        singleton = fc.get("Unassigned_Singleton", 0)
        total = assigned + no_features + singleton
        assigned_pct = (assigned / total * 100.0) if total else 0.0
        assignment_rows.append([
            sample,
            str(assigned),
            format_float(assigned_pct, 3),
            str(no_features),
            str(singleton),
        ])

    nedd4_row = next((row for row in tmm_rows if row.get("") == "ENSG00000069869.17"), None)
    nedd4_panel_row = next((row for row in fc_rows if row.get("Geneid") == "ENSG00000069869.17"), None)

    top_genes = sorted(fc_rows, key=lambda row: sum(int(row[sample]) for sample in samples), reverse=True)[:15]
    top_gene_rows = []
    for row in top_genes:
        total = sum(int(row[sample]) for sample in samples)
        top_gene_rows.append([row["Geneid"], str(total)] + [row[sample] for sample in samples])

    strongest_plus = max(signal_rows, key=lambda row: to_float(row, "rt_mean_depth"))
    weakest_plus = min(signal_rows, key=lambda row: to_float(row, "rt_mean_depth"))
    strongest_same_window_minus = max(signal_rows, key=lambda row: to_float(row, "same_window_minus_mean_depth"))
    strongest_paper_downstream = max(paper_rows.values(), key=lambda row: to_float(row, "paper_downstream_sense_mean_depth"))
    strongest_paper_upstream = max(paper_rows.values(), key=lambda row: to_float(row, "paper_upstream_antisense_mean_depth"))

    narrative = [
        f"Strongest downstream plus-strand signal is {strongest_plus['srr']} with mean depth {to_float(strongest_plus, 'rt_mean_depth'):.4f}.",
        f"Weakest downstream plus-strand signal is {weakest_plus['srr']} with mean depth {to_float(weakest_plus, 'rt_mean_depth'):.4f}.",
        f"The same-window minus-strand host-background view peaks in {strongest_same_window_minus['srr']} at mean depth {to_float(strongest_same_window_minus, 'same_window_minus_mean_depth'):.4f}.",
        f"The paper-style mate1 downstream sense proxy is strongest in {strongest_paper_downstream['srr']} at mean depth {to_float(strongest_paper_downstream, 'paper_downstream_sense_mean_depth'):.4f}.",
        f"The paper-style mate1 upstream antisense proxy peaks in {strongest_paper_upstream['srr']} at mean depth {to_float(strongest_paper_upstream, 'paper_upstream_antisense_mean_depth'):.4f}.",
    ]
    if nedd4_row is not None:
        narrative.append(
            "NEDD4 panel-normalized CPMs are "
            + ", ".join(f"{sample}={float(nedd4_row[sample]):.2f}" for sample in samples)
            + "."
        )

    report_data = {
        "samples": samples,
        "positions": depth_positions(depth_data[samples[0]]["plus"], downstream_start),
        "upstreamPositions": depth_positions(paper_depth_data[samples[0]]["upstream_antisense"], upstream_start),
        "depth": depth_data,
        "paperDepth": paper_depth_data,
        "rtSignal": {
            sample: {
                "plus": to_float(signal_by_sample[sample], "rt_mean_depth"),
                "same_window_minus": to_float(signal_by_sample[sample], "same_window_minus_mean_depth"),
                "paper_downstream": to_float(paper_rows[sample], "paper_downstream_sense_mean_depth"),
                "paper_upstream_antisense": to_float(paper_rows[sample], "paper_upstream_antisense_mean_depth"),
            }
            for sample in samples
        },
        "uniqueMapped": {sample: to_float(qc_rows[sample], "pct_unique_mapped") for sample in samples},
        "nedd4Tmm": {sample: float(nedd4_row[sample]) if nedd4_row else 0.0 for sample in samples},
        "nedd4Counts": {sample: int(nedd4_panel_row[sample]) if nedd4_panel_row else 0 for sample in samples},
    }

    html_text = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"utf-8\">
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">
  <title>MCF7 Readthrough Report</title>
  <style>
    :root {{
      --bg: #f4efe6;
      --panel: rgba(255, 250, 242, 0.9);
      --ink: #1f1a17;
      --muted: #695c54;
      --line: rgba(71, 52, 35, 0.16);
      --accent: #bd4b32;
      --accent-2: #2a7f62;
      --accent-3: #286ea8;
      --accent-4: #8b6fb8;
      --shadow: 0 18px 50px rgba(73, 49, 28, 0.12);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: Georgia, \"Times New Roman\", serif;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(189, 75, 50, 0.12), transparent 30%),
        radial-gradient(circle at top right, rgba(42, 127, 98, 0.14), transparent 28%),
        linear-gradient(180deg, #f8f2e8 0%, var(--bg) 42%, #efe7db 100%);
    }}
    .wrap {{ max-width: 1280px; margin: 0 auto; padding: 32px 20px 48px; }}
    .hero {{ padding: 28px; border: 1px solid var(--line); border-radius: 24px; background: linear-gradient(135deg, rgba(255,255,255,0.86), rgba(249, 238, 226, 0.92)); box-shadow: var(--shadow); }}
    h1, h2, h3 {{ margin: 0 0 12px; }}
    h1 {{ font-size: clamp(2rem, 5vw, 4rem); line-height: 0.96; letter-spacing: -0.04em; }}
    h2 {{ font-size: 1.5rem; }}
    p {{ margin: 0; line-height: 1.55; }}
    .sub {{ color: var(--muted); max-width: 76ch; margin-top: 14px; }}
    .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(250px, 1fr)); gap: 18px; margin-top: 22px; }}
    .card {{ background: var(--panel); border: 1px solid var(--line); border-radius: 20px; padding: 18px; box-shadow: var(--shadow); }}
    .stat-label {{ color: var(--muted); font-size: 0.95rem; margin-bottom: 6px; }}
    .stat-value {{ font-size: 2rem; font-weight: 700; letter-spacing: -0.04em; }}
    .section {{ margin-top: 26px; }}
    .plot-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 18px; margin-top: 16px; }}
    .chart-card {{ background: var(--panel); border: 1px solid var(--line); border-radius: 20px; padding: 16px; box-shadow: var(--shadow); }}
    .chart-wrap {{ width: 100%; overflow-x: auto; }}
    svg {{ width: 100%; height: auto; display: block; }}
    .legend {{ display: flex; gap: 14px; flex-wrap: wrap; color: var(--muted); font-size: 0.92rem; margin-top: 10px; }}
    .key {{ display: inline-flex; align-items: center; gap: 8px; }}
    .swatch {{ width: 12px; height: 12px; border-radius: 999px; display: inline-block; }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.94rem; }}
    th, td {{ padding: 10px 12px; border-bottom: 1px solid var(--line); text-align: left; }}
    th {{ color: var(--muted); font-weight: 600; }}
    tbody tr:hover {{ background: rgba(255,255,255,0.5); }}
    .callouts {{ display: grid; gap: 10px; margin-top: 16px; }}
    .callout {{ padding: 14px 16px; border-left: 4px solid var(--accent); background: rgba(255,255,255,0.55); border-radius: 12px; }}
    .small {{ color: var(--muted); font-size: 0.9rem; }}
    @media (max-width: 720px) {{
      .wrap {{ padding: 18px 14px 28px; }}
      .hero, .card, .chart-card {{ border-radius: 16px; }}
    }}
  </style>
</head>
<body>
  <div class=\"wrap\">
    <section class=\"hero\">
      <h1>MCF7 L1-NEDD4 review</h1>
      <p class=\"sub\">This report separates two different views of the locus: the Philippe-style mate1-only downstream sense plus upstream antisense windows, and a separate same-window opposite-strand host-background view retained only for local context.</p>
      <div class=\"grid\">
        <div class=\"card\"><div class=\"stat-label\">Samples</div><div class=\"stat-value\">4</div></div>
        <div class=\"card\"><div class=\"stat-label\">Downstream window</div><div class=\"stat-value\">1001 bp</div></div>
        <div class=\"card\"><div class=\"stat-label\">Panel genes</div><div class=\"stat-value\">100</div></div>
        <div class=\"card\"><div class=\"stat-label\">Top downstream plus sample</div><div class=\"stat-value\">{html.escape(strongest_plus['srr'])}</div></div>
      </div>
      <div class=\"callouts\">{''.join(f'<div class="callout">{html.escape(line)}</div>' for line in narrative)}</div>
    </section>

    <section class=\"section\">
      <h2>Paper-Style Summary</h2>
      <p class=\"small\">This is the Philippe-style comparison to focus on for strand interpretation: mate1 downstream sense in the downstream 1 kb window versus mate1 upstream antisense in the separate upstream 1 kb window.</p>
      <div class=\"card\">{table_html([
          'Sample', 'Downstream window', 'Upstream window', 'Mate1 downstream sense mean', 'Mate1 downstream sense CPM', 'Mate1 upstream antisense mean', 'Mate1 upstream antisense CPM'
      ], paper_summary_rows)}</div>
    </section>

    <section class=\"section\">
      <h2>Readthrough Summary</h2>
      <p class=\"small\">This table is retained for context, but the coverage plots below are limited to the paper-style downstream sense and upstream antisense windows only.</p>
      <div class=\"card\">{table_html([
          'Sample', 'RT mean depth (+)', 'RT CPM (+)', 'Same-window minus host background mean', 'Same-window minus host background CPM', 'Background mean', 'log2(RT/bg)', 'Unique mapped %'
      ], signal_summary_rows)}</div>
    </section>

    <section class=\"section\">
      <h2>Paper-Style Windows</h2>
      <p class=\"small\">These plots use mate 1 only, matching the paper more closely. Downstream sense and upstream antisense are shown as separate curves on their respective 1 kb windows. In the current MCF7 controls, the upstream antisense track is effectively zero across all four samples.</p>
      <div id=\"paper-coverage-plots\" class=\"plot-grid\"></div>
    </section>

    <section class=\"section\">
      <h2>Troubleshooting Note</h2>
      <div class=\"card\">
        <p>The remaining mismatch with Philippe Fig. 3B is upstream of the HTML plotting layer. In this dataset, pooled untreated MCF7 controls still only reach downstream max depth around 5 in the current 1 kb window even after broadening read inclusion, while upstream antisense remains at 0. The selected downstream sense reads are often heavily spliced and contribute only short matched segments within the target interval, so aligned-base pileup stays low. That points to a deeper difference from the paper, most likely a different plotted interval or a different track-construction rule.</p>
      </div>
    </section>

    <section class=\"section\">
      <h2>FeatureCounts Assignment</h2>
      <div class=\"card\">{table_html(['Sample', 'Assigned reads', 'Assigned %', 'NoFeatures', 'Singleton'], assignment_rows)}</div>
    </section>

    <section class=\"section\">
      <h2>NEDD4 Snapshot</h2>
      <div class=\"grid\">
        <div class=\"card\">{table_html(['Sample', 'Panel counts', 'TMM CPM'], [[sample, str(report_data['nedd4Counts'][sample]), format_float(report_data['nedd4Tmm'][sample], 2)] for sample in samples])}</div>
        <div class=\"card\">
          <h3>Interpretation</h3>
          <p>NEDD4 remains detectable in all but the lowest-count sample at the raw-count level and spans roughly 0 to 3 panel counts in this restricted assay, while the TMM CPM scale stays more comparable across libraries after normalization.</p>
        </div>
      </div>
    </section>

    <section class=\"section\">
      <h2>Top Panel Genes</h2>
      <div class=\"card\">{table_html(['Gene ID', 'Total counts'] + samples, top_gene_rows)}</div>
    </section>
  </div>

  <script>
    const reportData = {json.dumps(report_data)};
    const colors = {{
      plus: '#bd4b32',
      minus: '#2a7f62',
      bars: ['#bd4b32', '#2a7f62', '#286ea8', '#8b6fb8']
    }};

    function linePath(values, width, height, padding) {{
      const maxVal = Math.max(...values, 1);
      const plotWidth = width - padding.left - padding.right;
      const plotHeight = height - padding.top - padding.bottom;
      return values.map((value, index) => {{
        const x = padding.left + (index / (values.length - 1 || 1)) * plotWidth;
        const y = height - padding.bottom - (value / maxVal) * plotHeight;
        return `${{index === 0 ? 'M' : 'L'}}${{x.toFixed(2)}},${{y.toFixed(2)}}`;
      }}).join(' ');
    }}

    function axisTicks(maxVal, count) {{
      const ticks = [];
      for (let i = 0; i <= count; i += 1) {{
        ticks.push((maxVal / count) * i);
      }}
      return ticks;
    }}

    function xLabels(start, width, height, padding) {{
      const plotWidth = width - padding.left - padding.right;
      return [0, 250, 500, 750, 1000].map((offset) => {{
        const x = padding.left + (offset / 1000) * plotWidth;
        return `<text x="${{x}}" y="${{height - 8}}" text-anchor="middle" font-size="11" fill="#695c54">${{start + offset}}</text>`;
      }}).join('');
    }}

    function paperCoverageCard(sample) {{
      const downstream = reportData.paperDepth[sample].downstream_sense;
      const upstream = reportData.paperDepth[sample].upstream_antisense;
      const values = downstream.concat(upstream);
      const maxVal = Math.max(...values, 1);
      const width = 640;
      const height = 280;
      const padding = {{ top: 16, right: 14, bottom: 30, left: 42 }};
      const plotHeight = height - padding.top - padding.bottom;
      const ticks = axisTicks(maxVal, 4);
      const yGrid = ticks.map((tick) => {{
        const y = height - padding.bottom - (tick / maxVal) * plotHeight;
        return `<line x1="${{padding.left}}" x2="${{width - padding.right}}" y1="${{y}}" y2="${{y}}" stroke="rgba(71,52,35,0.15)" stroke-dasharray="3 4"></line><text x="${{padding.left - 8}}" y="${{y + 4}}" text-anchor="end" font-size="11" fill="#695c54">${{tick.toFixed(0)}}</text>`;
      }}).join('');
      return `<div class="chart-card"><h3>${{sample}}</h3><div class="chart-wrap"><svg viewBox="0 0 ${{width}} ${{height}}" role="img" aria-label="Paper-style windows for ${{sample}}">${{yGrid}}<line x1="${{padding.left}}" x2="${{padding.left}}" y1="${{padding.top}}" y2="${{height - padding.bottom}}" stroke="#695c54"></line><line x1="${{padding.left}}" x2="${{width - padding.right}}" y1="${{height - padding.bottom}}" y2="${{height - padding.bottom}}" stroke="#695c54"></line><path d="${{linePath(downstream, width, height, padding)}}" fill="none" stroke="${{colors.plus}}" stroke-width="2.25"></path><path d="${{linePath(upstream, width, height, padding)}}" fill="none" stroke="${{colors.minus}}" stroke-width="2.25" stroke-dasharray="6 4"></path>${{xLabels(reportData.positions[0], width, height, padding)}}</svg></div><div class="legend"><span class="key"><span class="swatch" style="background:${{colors.plus}}"></span>Mate1 downstream sense</span><span class="key"><span class="swatch" style="background:${{colors.minus}}"></span>Mate1 upstream antisense</span></div><p class="small">Upstream window starts at ${{reportData.upstreamPositions[0]}}. Downstream window labels are shown on the axis.</p></div>`;
    }}

    document.getElementById('paper-coverage-plots').innerHTML = reportData.samples.map(paperCoverageCard).join('');
  </script>
</body>
</html>
"""

    OUT_HTML.write_text(html_text)
    print(f"wrote {OUT_HTML}")


if __name__ == "__main__":
    main()