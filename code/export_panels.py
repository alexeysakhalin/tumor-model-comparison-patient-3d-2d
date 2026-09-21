"""Export Figure 6 panels as separate print-resolution PNG files.

The composite and standalone exports use the same panel functions and numerical inputs.
Standalone canvases provide additional space for labels and legends. Each PNG records its
scientific description, author, source-table checksum prefixes, license, repository and software.

Usage: python code/export_panels.py tables build/panels [DPI]
"""
import hashlib
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import make_full_figure as F
from table_io import read_table

MM = 1 / 25.4
DPI = 600
REPOSITORY = "https://github.com/alexeysakhalin/tumor-model-comparison-patient-3d-2d"
AUTHOR = "Alexey Petukhov, Nazarbayev University (ORCID 0000-0002-7298-5238)"
LICENCE = "Creative Commons Attribution 4.0 International (CC BY 4.0)"

PANELS = {
    "a": {
        "size": (195.0, 112.0),
        "title": "Figure 6a. Twelve gene sets across four expression layers, seven CRISPR arms and one published list",
        "description": (
            "Left block: expression-rank difference against level-matched control genes, malignant minus "
            "healthy and 3D dome minus 2D, in reference-band units. Middle block: CRISPR knockout effect "
            "under immune pressure in seven arms of five models. Right block: the share of each gene set "
            "in the published candidate list of Patel et al. 2017, a count and not an enrichment test. "
            "A solid frame marks a value beyond the 95 per cent reference band; a dotted frame marks a "
            "value beyond it where the sources disagree."
        ),
        "tables": ["IKM_representation_in_band_units.csv", "IKM_models_panoptosis_layers.csv",
                   "IKM_patel2017_group_overrepresentation.csv", "IKM_gene_membership_by_group.tsv"],
    },
    "b": {
        "size": (95.0, 92.0),
        "title": "Figure 6b. Cross-platform comparison in HT-29: siRNA viability rescue against CRISPR resistance",
        "description": (
            "Each point is one gene of the kinome library of Woznicki et al. 2021, scored by viability "
            "rescue under IFN-gamma with TNF-alpha, against its resistance percentile in the CRISPR "
            "knockout screen of Watterson et al. 2026 under IFN-gamma. The two screens use different "
            "cytokine conditions; the shaded band marks the top decile of the CRISPR ranking."
        ),
        "tables": ["IKM_ht29_rnai_crossplatform.csv"],
    },
    "c": {
        "size": (140.0, 168.0),
        "title": "Figure 6c. Expression-rank differences with cluster-bootstrap intervals, and functional magnitude",
        "description": (
            "Upper axes: the gene-set statistic in healthy epithelium, patient tumour and 3D dome "
            "cultures as a difference from level-matched control genes, in percentile points, with 95 "
            "per cent cluster-bootstrap intervals; an open marker means the interval overlaps zero. The "
            "2D layer is the reference of the comparison and is close to zero by construction rather "
            "than exactly zero; its individual estimates are nonzero and are reported in the tables. "
            "Lower axes: median absolute gene-set statistic across the seven CRISPR arms in "
            "reference-band units, with the number of arms outside the band."
        ),
        "tables": ["IKM_final_master_table_CI.csv", "IKM_group_weight_all12.csv",
                   "IKM_gene_membership_by_group.tsv"],
    },
}


def separate_touching_text(fig, margin_pt=1.5, passes=3):
    """Push apart text boxes that touch, moving the lower one in its own coordinate system.

    Standalone panels inherit the tight key spacing of the assembled page, where the neighbouring
    panels absorb part of the vertical room. Each colliding pair is separated by the measured
    overlap plus a small margin, converted into the transform the text was created in, so a label
    positioned in data coordinates is not displaced into the matrix it annotates.
    """
    for _ in range(passes):
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        items = [(t, t.get_window_extent(renderer)) for t in fig.findobj(plt.matplotlib.text.Text)
                 if t.get_text().strip() and t.get_visible()]
        moved = 0
        for i, (a, box_a) in enumerate(items):
            for b, box_b in items[i + 1:]:
                if not box_a.overlaps(box_b):
                    continue
                lower = b if box_b.y0 <= box_a.y0 else a
                other = box_a if lower is b else box_b
                shift_px = (other.y0 - lower.get_window_extent(renderer).y1) - margin_pt
                if shift_px >= 0:
                    continue
                transform = lower.get_transform()
                x, y = lower.get_position()
                origin = transform.transform((x, y))
                target = transform.inverted().transform((origin[0], origin[1] + shift_px))
                lower.set_position((x, float(target[1])))
                moved += 1
        if not moved:
            break
    return fig


def table_digests(tables_dir, names):
    out = []
    for name in names:
        path = Path(tables_dir) / name
        if path.exists():
            out.append(f"{name} {hashlib.sha256(path.read_bytes()).hexdigest()[:12]}")
    return "; ".join(out)


def metadata(panel, tables_dir):
    spec = PANELS[panel]
    return {
        "Title": spec["title"],
        "Author": AUTHOR,
        "Description": spec["description"],
        "Copyright": LICENCE,
        "Source": REPOSITORY,
        "Software": "make_full_figure.py and export_panels.py; Matplotlib "
                    + plt.matplotlib.__version__,
        "Comment": "Drawn from " + table_digests(tables_dir, spec["tables"])
                   + " (SHA-256 prefixes).",
    }


def export(tables_dir="tables", out_dir="build/panels", dpi=DPI):
    Path(out_dir).mkdir(parents=True, exist_ok=True)
    F.load_labels(tables_dir)
    C = read_table(f"{tables_dir}/IKM_final_master_table_CI.csv").set_index("group")
    Mm = read_table(f"{tables_dir}/IKM_models_panoptosis_layers.csv")
    Pp = read_table(f"{tables_dir}/IKM_patel2017_group_overrepresentation.csv").set_index("group")
    Rr = read_table(f"{tables_dir}/IKM_representation_in_band_units.csv")
    CPd = read_table(f"{tables_dir}/IKM_ht29_rnai_crossplatform.csv")
    Ww = read_table(f"{tables_dir}/IKM_group_weight_all12.csv").set_index("group")

    written, figures = [], {}

    w, h = PANELS["a"]["size"]
    fig_a = plt.figure(figsize=(w * MM, h * MM))
    ax_a = fig_a.add_subplot(111)
    F.panel_matrix_all(ax_a, Rr, Mm, Pp)
    # The key of this panel is drawn in the axes of the matrix, so it is left exactly where the
    # composite puts it: moving text that is positioned in data coordinates would land it inside the
    # matrix. Standalone legibility is handled by the canvas size and the saved margins instead.
    ax_a.text(-0.055, 1.03, "a", transform=ax_a.transAxes, fontsize=9.5, fontweight="bold",
              va="top", ha="left", color=F.INK)
    # Retain data-coordinate alignment of the pressure labels; check text boxes after rendering.
    figures["a"] = fig_a

    w, h = PANELS["b"]["size"]
    fig_b = plt.figure(figsize=(w * MM, h * MM))
    ax_b = fig_b.add_subplot(111)
    F.panel_hits(ax_b, CPd)
    ax_b.text(-0.29, 1.045, "b", transform=ax_b.transAxes, fontsize=9.5, fontweight="bold",
              va="top", ha="left", color=F.INK)
    figures["b"] = fig_b

    w, h = PANELS["c"]["size"]
    fig_c = plt.figure(figsize=(w * MM, h * MM))
    gs_c = GridSpec(2, 1, figure=fig_c, height_ratios=[62, 38], hspace=0.06,
                    left=0.215, right=0.980, top=0.862, bottom=0.205)
    ax_lev = fig_c.add_subplot(gs_c[0, 0])
    ax_wt = fig_c.add_subplot(gs_c[1, 0], sharex=ax_lev)
    F.panel_levels_weight(ax_lev, ax_wt, C, Ww)
    # The two long rotated axis labels meet in a single left margin when the panel stands alone, and
    # the figure-level key lands on the domain headings; both are given their own space here.
    ax_lev.yaxis.set_label_coords(-0.075, 0.5)
    ax_wt.yaxis.set_label_coords(-0.185, 0.5)
    if fig_c.legends:
        fig_c.legends[0].set_bbox_to_anchor((0.5, 0.962), transform=fig_c.transFigure)
    ax_lev.text(-0.205, 1.09, "c", transform=ax_lev.transAxes, fontsize=9.5, fontweight="bold",
                va="top", ha="left", color=F.INK)
    figures["c"] = fig_c

    for panel, fig in figures.items():
        path = os.path.join(out_dir, f"Figure6_panel_{panel}.png")
        fig.savefig(path, dpi=dpi, bbox_inches="tight", metadata=metadata(panel, tables_dir))
        written.append(path)

    # Record overlapping text boxes for all three standalone canvases.
    report = os.path.join(out_dir, "PANEL_TEXT_OVERLAPS.tsv")
    with open(report, "w", encoding="utf-8") as fh:
        fh.write("panel\toverlapping_pairs\tpairs\n")
        for panel, fig in figures.items():
            fig.canvas.draw()
            renderer = fig.canvas.get_renderer()
            boxes = [(t, t.get_window_extent(renderer)) for t in fig.findobj(plt.matplotlib.text.Text)
                     if t.get_text().strip() and t.get_visible()]
            pairs = [f"{a.get_text()[:24]} / {b.get_text()[:24]}".replace("\n", " ")
                     for i, (a, ba) in enumerate(boxes) for b, bb in boxes[i + 1:] if ba.overlaps(bb)]
            fh.write(f"{panel}\t{len(pairs)}\t{'; '.join(pairs)}\n")
            print(f"panel {panel}: text overlaps {len(pairs)}", flush=True)
    written.append(report)
    for path in written:
        print("written:", path, os.path.getsize(path), "bytes", flush=True)
    return figures, written


if __name__ == "__main__":
    arguments = sys.argv[1:]
    export(arguments[0] if arguments else "tables",
           arguments[1] if len(arguments) > 1 else "build/panels",
           int(arguments[2]) if len(arguments) > 2 else DPI)
