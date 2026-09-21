"""Render Figure 6 from the deposited derived tables.

Usage: python code/make_full_figure.py tables figures/Figure6
"""

import os
import sys
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from table_io import read_table
from matplotlib.gridspec import GridSpec

MM = 1 / 25.4
LAY = {
    "healthy": "#0f6f78",
    "patient": "#9e2b4e",
    "organoid": "#b07316",
    "line": "#3d4a5c",
}
COL = {"healthy": "healthy", "patient": "malignant", "organoid": "dome", "line": "2D"}
CI_FILE = "IKM_final_master_table_CI.csv"
FLAG = {
    "healthy": "healthy_outside_reference_band",
    "patient": "malignant_outside_reference_band",
    "organoid": "dome_outside_reference_band",
}
LAY_EN = {
    "healthy": "healthy tissue",
    "patient": "patient tumour",
    "organoid": "organoid (3D dome)",
    "line": "2D cell line",
}
INK, SOFT = ("#1a1a1a", "#5a5a5a")
BLOCKS = [
    ("Receptors", ["ENTRY receptors", "ENTRY adaptors"], "#1f6fb2"),
    ("Inhibition", ["BRAKES"], "#2e7d5b"),
    (
        "Death machinery",
        [
            "Apoptosis initiators",
            "Apoptosis executioners",
            "Necroptosis initiators",
            "Necroptosis executioners",
            "Pyroptosis initiators",
            "Pyroptosis executioners",
        ],
        "#7b3f9d",
    ),
    (
        "Immune recognition",
        ["MHC-I presentation", "NK activating ligands", "Inhibitory signals"],
        "#0f6f78",
    ),
]
ORD = [g for _, _gs, _ in BLOCKS for g in _gs]
LABEL = {}


def load_labels(tables_dir):
    """Require one non-empty display label for each plotted group."""
    gm = read_table(
        os.path.join(tables_dir, "IKM_gene_membership_by_group.tsv"), sep="\t"
    )
    if gm["figure_label"].isna().any() or gm["figure_label"].str.strip().eq("").any():
        raise ValueError("Every gene assignment must carry a figure label")
    labels = gm[["group", "figure_label"]].drop_duplicates()
    if labels["group"].duplicated().any() or set(labels["group"]) != set(ORD):
        raise ValueError(
            "Exactly one display label is required for every plotted group"
        )
    LABEL.clear()
    LABEL.update(
        {
            k: v.replace("\\n", "\n")
            for k, v in zip(labels["group"], labels["figure_label"])
        }
    )
    return LABEL


def _rows():
    order, ypos, block_of, y = ([], {}, {}, 0.0)
    for bname, groups, bcol in BLOCKS:
        for g in groups:
            order.append(g)
            ypos[g] = y
            block_of[g] = (bname, bcol)
            y -= 1.0
        y -= 0.55
    return (order, ypos, block_of)


def panel_levels_weight(axT, axB, C, W):
    """Plot expression intervals and descriptive CRISPR magnitudes."""
    xpos, x = ({}, 0.0)
    for _, groups, _ in BLOCKS:
        for g in groups:
            xpos[g] = x
            x += 1.0
        x += 0.8
    xs = [xpos[g] for g in ORD]
    axT.set_xlim(-0.9, max(xs) + 0.9)
    axB.set_xlim(-0.9, max(xs) + 0.9)
    for k, g in enumerate(ORD):
        if k % 2 == 0:
            for a in (axT, axB):
                a.axvspan(xpos[g] - 0.5, xpos[g] + 0.5, color="#f6f6f6", zorder=0)
    axT.axhline(0, color="#4a4a4a", lw=1.0, zorder=2)
    for g in ORD:
        r = C.loc[g]
        for lay, dx in [("healthy", -0.26), ("patient", 0.0), ("organoid", 0.26)]:
            v = float(r[COL[lay]])
            lo, hi = (float(r[f"{COL[lay]}_lo"]), float(r[f"{COL[lay]}_hi"]))
            solid = not lo <= 0 <= hi
            axT.plot(
                [xpos[g] + dx, xpos[g] + dx], [lo, hi], color=LAY[lay], lw=1.1, zorder=3
            )
            axT.scatter(
                xpos[g] + dx,
                v,
                s=24,
                zorder=4,
                lw=1.1,
                facecolor=LAY[lay] if solid else "white",
                edgecolor=LAY[lay],
            )
    for g in ORD:
        wv = float(W.loc[g, "magnitude"])
        n_out = int(W.loc[g, "arms_outside_band"])
        col = "#2e7d5b" if n_out >= 4 else "#9a9a9a"
        axB.bar(xpos[g], wv, width=0.62, color=col, edgecolor="none", zorder=3)
        axB.text(
            xpos[g],
            wv + 0.07,
            f"{wv:.2f}",
            fontsize=7.0,
            ha="center",
            va="bottom",
            color=INK,
            fontweight="bold" if n_out >= 4 else "normal",
        )
        axB.text(
            xpos[g],
            -0.05,
            f"{n_out}/7",
            fontsize=7.0,
            ha="center",
            va="top",
            color=SOFT,
        )
    for a, ymaxfun in ((axT, None), (axB, None)):
        a.set_xticks([])
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    axB.spines["bottom"].set_visible(False)
    axT.set_ylim(-24, 26)
    axT.set_yticks([-20, -10, 0, 10, 20])
    axT.tick_params(axis="y", labelsize=6.0, colors=SOFT, length=2.5)
    axB.tick_params(axis="y", labelsize=6.0, colors=SOFT, length=2.5)
    axB.set_yticks([0.0, 0.5, 1.0, 1.5])
    axT.set_ylabel(
        "Expression-rank difference vs matched\ncontrols, percentile points",
        fontsize=8.0,
        color=INK,
        labelpad=2,
        linespacing=1.25,
    )
    axB.set_ylabel(
        "Median absolute gene-set\nstatistic, reference-band units",
        fontsize=7.6,
        color=INK,
        labelpad=3,
        linespacing=1.3,
    )
    axB.set_ylim(-3.25, max((float(W.loc[g, "magnitude"]) for g in ORD)) * 1.22)
    for g in ORD:
        col = next((c for _, gs_, c in BLOCKS if g in gs_))
        lab = LABEL.get(g, g)
        head, _, tail = lab.rpartition(" ")
        two_line = (
            f"{head}\n{tail} ({int(C.loc[g, 'gene_count'])})"
            if head
            else f"{lab} ({int(C.loc[g, 'gene_count'])})"
        )
        axB.text(
            xpos[g],
            -0.6,
            two_line,
            fontsize=6.4,
            rotation=90,
            ha="center",
            va="top",
            color=col,
            linespacing=1.15,
        )
    ytop = axT.get_ylim()[1]
    for bi, (bname, groups, bcol) in enumerate(BLOCKS):
        gx = [xpos[g] for g in groups]
        axT.plot(
            [min(gx) - 0.42, max(gx) + 0.42],
            [ytop * 1.02, ytop * 1.02],
            color=bcol,
            lw=2.0,
            solid_capstyle="butt",
            clip_on=False,
            zorder=3,
        )
        y_h = ytop * (1.06 if bi % 2 == 0 else 1.17)
        axT.text(
            (min(gx) + max(gx)) / 2,
            y_h,
            bname,
            fontsize=8.4,
            ha="center",
            va="bottom",
            color=bcol,
            fontweight="bold",
            clip_on=False,
        )
    handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markersize=6.0,
            markerfacecolor=LAY[k_],
            markeredgecolor=LAY[k_],
            label=t_,
        )
        for k_, t_ in [
            ("healthy", "healthy tissue"),
            ("patient", "patient tumour"),
            ("organoid", "organoid (3D dome)"),
        ]
    ]
    handles.append(
        Line2D(
            [],
            [],
            marker="o",
            linestyle="",
            markersize=6.0,
            markerfacecolor="white",
            markeredgecolor=INK,
            label="95% interval overlaps zero",
        )
    )
    handles.append(
        Line2D([], [], color="#4a4a4a", lw=1.0, label="2D reference, centred at 0")
    )
    pos_t = axT.get_position()
    axT.get_figure().legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, pos_t.y1 + 0.036),
        ncol=5,
        frameon=False,
        fontsize=7.0,
        handletextpad=0.5,
        columnspacing=1.5,
        borderpad=0.0,
        labelspacing=0.2,
    )


def panel_matrix_all(ax, R, M, P=None):
    """One matrix: two representation differences, the seven functional screen arms, and one
    annotation column for an independent published hit list.

    The annotation column is deliberately not a matrix cell: it carries no band units, because the
    source published only a truncated ranked tail and no counts. It is drawn as a separate column
    with its own mark and its own caveat.
    """
    order, ypos, _ = _rows()
    repr_cols = ["patient − healthy", "organoid − 2D"]
    repr_labels = ["malignant − healthy", "organoid − 2D"]
    arms = [
        "SK-MEL-2 (IFN-γ)",
        "A375 (IFN-γ)",
        "HT-29 (IFN-γ)",
        "HT-29 JAK1-KO (IFN-γ)",
        "CRC-9 organoid (IFN-γ)",
        "CRC-9 organoid (T cells)",
        "HeLa (TNF-α + IFN-γ)",
    ]
    arm_short = [
        "SK-MEL-2",
        "A375",
        "HT-29",
        "HT-29 JAK1-KO",
        "CRC-9 tumoroid",
        "CRC-9 tumoroid\n+ autologous T cells",
        "HeLa",
    ]
    fmt = ["2D", "2D", "2D", "2D", "3D", "3D", "2D"]
    gap = 0.9
    xs = list(range(len(repr_cols))) + [
        len(repr_cols) + gap + j for j in range(len(arms))
    ]
    labels = repr_labels + arm_short
    rv = R.pivot_table(index="group", columns="column", values="reference_band_units")
    ro = R.pivot_table(
        index="group", columns="column", values="outside_band", aggfunc="first"
    )
    rc = R.pivot_table(
        index="group", columns="column", values="all_sources_agree", aggfunc="first"
    )
    mv = M.pivot_table(index="group", columns="model", values="reference_band_units")
    mo = M.pivot_table(
        index="group", columns="model", values="outside_band", aggfunc="first"
    )
    mc = M.pivot_table(
        index="group", columns="model", values="all_replicates_agree", aggfunc="first"
    )
    ax.set_ylim(min(ypos.values()) - 4.95, 2.1)
    xp = xs[-1] + 1.9
    ax.set_xlim(-4.6, xp + 0.9)
    ax.set_yticks([])
    ax.set_xticks([])
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    for j, (x, col) in enumerate(zip(xs, repr_cols + arms)):
        src = (rv, ro, rc) if j < len(repr_cols) else (mv, mo, mc)
        for g in order:
            if g not in src[0].index or col not in src[0].columns:
                continue
            v = src[0].loc[g, col]
            if v != v:
                continue
            colour = plt.get_cmap("RdBu_r")((max(min(float(v), 3.0), -3.0) + 3.0) / 6.0)
            ax.add_patch(
                plt.Rectangle(
                    (x - 0.42, ypos[g] - 0.34),
                    0.84,
                    0.68,
                    facecolor=colour,
                    edgecolor="none",
                    zorder=2,
                )
            )
            if bool(src[1].loc[g, col]):
                solid = (
                    bool(src[2].loc[g, col])
                    if src[2].loc[g, col] == src[2].loc[g, col]
                    else True
                )
                ax.add_patch(
                    plt.Rectangle(
                        (x - 0.42, ypos[g] - 0.34),
                        0.84,
                        0.68,
                        facecolor="none",
                        edgecolor=INK,
                        lw=1.0 if solid else 0.7,
                        linestyle="-" if solid else (0, (1.4, 1.0)),
                        zorder=3,
                    )
                )
    if P is not None:
        for g in order:
            if g not in P.index:
                continue
            k = int(P.loc[g, "published_list_overlap"])
            n = int(P.loc[g, "gene_count"])
            if k:
                shade = plt.get_cmap("Greens")(0.18 + 0.72 * (k / n))
                ax.add_patch(
                    plt.Rectangle(
                        (xp - 0.42, ypos[g] - 0.34),
                        0.84,
                        0.68,
                        facecolor=shade,
                        edgecolor="none",
                        zorder=2,
                    )
                )
                ax.text(
                    xp,
                    ypos[g],
                    f"{k}/{n}",
                    fontsize=6.4,
                    ha="center",
                    va="center",
                    color="white" if k / n > 0.45 else INK,
                    zorder=4,
                )
            else:
                ax.add_patch(
                    plt.Rectangle(
                        (xp - 0.42, ypos[g] - 0.34),
                        0.84,
                        0.68,
                        facecolor="#f4f4f4",
                        edgecolor="#d8d8d8",
                        lw=0.6,
                        zorder=2,
                    )
                )
                ax.text(
                    xp,
                    ypos[g],
                    f"0/{n}",
                    fontsize=6.4,
                    ha="center",
                    va="center",
                    color="#9a9a9a",
                    zorder=4,
                )
    blk_of = {g: (nm, c) for nm, gs_, c in BLOCKS for g in gs_}
    for g in order:
        ax.text(
            -1.05,
            ypos[g],
            LABEL.get(g, g),
            fontsize=7.0,
            ha="right",
            va="center",
            color=blk_of[g][1],
        )
    focus = {"HT-29", "HT-29 JAK1-KO", "HeLa"}
    if P is not None:
        ax.text(
            xp,
            0.6,
            "Mel624 + TCR T cells",
            fontsize=8.0,
            rotation=90,
            ha="center",
            va="bottom",
            color=INK,
        )
    for x, lab in zip(xs, labels):
        hot = lab in focus
        ax.text(
            x,
            0.6,
            lab,
            fontsize=8.4 if hot else 8.0,
            rotation=90,
            ha="center",
            va="bottom",
            color=INK,
            fontweight="bold" if hot else "normal",
        )
    ybot = min(ypos.values())
    foot = M.drop_duplicates("model").set_index("model")
    for j, (x, arm) in enumerate(zip(xs[len(repr_cols) :], arms)):
        ax.text(
            x,
            ybot - 0.8,
            fmt[j],
            fontsize=6.0,
            ha="center",
            va="center",
            color=LAY["organoid"] if fmt[j] == "3D" else SOFT,
            fontweight="bold" if fmt[j] == "3D" else "normal",
        )
        ax.text(
            x,
            ybot - 1.62,
            str(int(foot.loc[arm, "screen_count"])),
            fontsize=6.0,
            ha="center",
            va="center",
            color=SOFT,
        )
        rho_v = foot.loc[arm, "replicate_rho"]
        rho_s = (
            "—"
            if rho_v != rho_v
            else format(round(float(rho_v), 2) + 0.0, "+.2f").replace("-0.00", "0.00")
        )
        ax.text(
            x, ybot - 2.44, rho_s, fontsize=6.0, ha="center", va="center", color=SOFT
        )
    ax.text(
        (xs[0] + xs[len(repr_cols) - 1]) / 2,
        ybot - 0.8,
        "6 organs",
        fontsize=6.0,
        ha="center",
        va="center",
        color=SOFT,
    )
    ax.text(
        -0.75, ybot - 0.8, "format", fontsize=6.0, ha="right", va="center", color=SOFT
    )
    ax.text(
        -0.75,
        ybot - 1.62,
        "screens pooled",
        fontsize=6.0,
        ha="right",
        va="center",
        color=SOFT,
    )
    ax.text(
        -0.75,
        ybot - 2.44,
        "replicate ρ",
        fontsize=6.0,
        ha="right",
        va="center",
        color=SOFT,
    )
    xm = (xs[len(repr_cols) - 1] + xs[len(repr_cols)]) / 2
    ax.plot([xm, xm], [ybot - 3.1, 1.9], color="#bbbbbb", lw=0.8)
    if P is not None:
        ax.plot(
            [(xs[-1] + xp) / 2, (xs[-1] + xp) / 2],
            [ybot - 3.1, 1.9],
            color="#bbbbbb",
            lw=0.8,
        )
        ax.plot(
            [xp - 0.42, xp + 0.42], [ybot - 3.1, ybot - 3.1], color="#999999", lw=0.8
        )
        ax.text(
            xp,
            ybot - 3.26,
            "published hit list, overlap",
            fontsize=6.0,
            ha="center",
            va="top",
            color=INK,
            fontweight="bold",
        )
        ax.text(
            xp,
            ybot - 3.95,
            "hit list of",
            fontsize=6.2,
            ha="center",
            va="top",
            color=INK,
        )
        ax.text(
            xp,
            ybot - 4.75,
            "Patel et al. 2017",
            fontsize=6.6,
            ha="center",
            va="top",
            color=INK,
            fontweight="bold",
        )
    for x0, x1, txt in [
        (xs[0], xs[len(repr_cols) - 1], "expression rank, tissue-matched"),
        (xs[len(repr_cols)], xs[-1], "CRISPR knockout, 7 arms"),
    ]:
        ax.plot(
            [x0 - 0.42, x1 + 0.42], [ybot - 3.1, ybot - 3.1], color="#999999", lw=0.8
        )
        ax.text(
            (x0 + x1) / 2,
            ybot - 3.26,
            txt,
            fontsize=6.2,
            ha="center",
            va="top",
            color=INK,
            fontweight="bold",
        )
    ax.text(
        (xs[0] + xs[len(repr_cols) - 1]) / 2,
        ybot - 3.95,
        "Tabula Sapiens 2.0, Kang et al.\npan-cancer atlas, DepMap",
        fontsize=6.6,
        ha="center",
        va="top",
        color=INK,
        fontweight="bold",
        linespacing=1.3,
    )
    ax.text(
        (xs[len(repr_cols)] + xs[-1]) / 2,
        ybot - 3.95,
        "reanalysis of Watterson et al. 2026\nand Zhou et al. 2026",
        fontsize=6.6,
        ha="center",
        va="top",
        color=INK,
        fontweight="bold",
        linespacing=1.3,
    )
    cax = ax.inset_axes([0.03, -0.125, 0.23, 0.03])
    cax.imshow(
        np.linspace(-3, 3, 200).reshape(1, -1),
        aspect="auto",
        cmap="RdBu_r",
        extent=[-3, 3, 0, 1],
    )
    cax.set_yticks([])
    cax.set_xticks([-3, 0, 3])
    cax.tick_params(labelsize=7.0, colors=INK, length=2.5, pad=1.5)
    for s in cax.spines.values():
        s.set_visible(False)
    cax.set_xlabel(
        "gene-set statistic, reference-band units", fontsize=6.8, color=INK, labelpad=2
    )
    for k, (x0, lw, style, txt) in enumerate(
        [
            (0.315, 1.0, "-", "beyond the 95% reference band"),
            (
                0.315,
                0.7,
                (0, (1.4, 1.0)),
                "beyond it, sources disagree:\norgans in the left block, replicates in the middle",
            ),
        ]
    ):
        yb = -0.072 - k * 0.042
        ax.add_patch(
            plt.Rectangle(
                (x0, yb),
                0.02,
                0.03,
                transform=ax.transAxes,
                facecolor="none",
                edgecolor=INK,
                lw=lw,
                linestyle=style,
                clip_on=False,
                zorder=4,
            )
        )
        ax.text(
            x0 + 0.032,
            yb + 0.015,
            txt,
            fontsize=6.6,
            transform=ax.transAxes,
            ha="left",
            va="center",
            color=INK,
        )
    if P is not None:
        gax = ax.inset_axes([0.872, -0.125, 0.128, 0.03])
        gax.imshow(
            np.linspace(0, 1, 200).reshape(1, -1),
            aspect="auto",
            cmap="Greens",
            extent=[0, 1, 0, 1],
            vmin=-0.25,
            vmax=1.18,
        )
        gax.set_yticks([])
        gax.set_xticks([0, 1])
        gax.set_xticklabels(["0", "all"])
        gax.tick_params(labelsize=7.0, colors=INK, length=2.5, pad=1.5)
        for s in gax.spines.values():
            s.set_visible(False)
        gax.text(
            1.0,
            -1.55,
            "share of the group in the published list",
            transform=gax.transAxes,
            fontsize=6.8,
            color=INK,
            ha="right",
            va="top",
        )


SYS_HOT = "#9e2b4e"
SYS_COOL = "#3f5d75"


def panel_hits(ax, CP):
    d = CP[CP["in_screen"]].copy()
    hc = d[d["high_confidence_published_hit"]]
    lo = d[~d["high_confidence_published_hit"]]
    ax.axhspan(90, 100, color="#e8e8e8", zorder=0)
    ax.axhline(50, color="#bbbbbb", lw=0.8, zorder=1)
    ax.scatter(
        lo["RNAi_Z"],
        lo["percentile_WT"],
        s=13,
        facecolor="white",
        edgecolor=SYS_COOL,
        lw=0.9,
        zorder=3,
    )
    ax.scatter(hc["RNAi_Z"], hc["percentile_WT"], s=28, color=SYS_HOT, zorder=4)
    for k, (_, r) in enumerate(hc.sort_values("RNAi_Z").iterrows()):
        ax.annotate(
            f"$\\it{{{r['gene']}}}$",
            (r["RNAi_Z"], r["percentile_WT"]),
            textcoords="offset points",
            xytext=(-5, -11) if k == 0 else (4, 5),
            fontsize=6.2,
            color=SYS_HOT,
            fontweight="bold",
            ha="right" if k == 0 else "left",
        )
    ax.set_xlabel(
        "siRNA viability-rescue Z score\n" + "$\\bf{Woznicki\\ et\\ al.\\ 2021}$",
        fontsize=8.0,
        color=INK,
        labelpad=3,
        linespacing=1.5,
    )
    ax.set_ylabel(
        "CRISPR resistance percentile\n" + "$\\bf{Watterson\\ et\\ al.\\ 2026}$",
        fontsize=7.8,
        color=INK,
        labelpad=3,
        linespacing=1.5,
    )
    ax.set_ylim(-3, 108)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.tick_params(labelsize=6.2, colors=SOFT)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


def build(tables_dir="tables", out_prefix="build/Figure6"):
    load_labels(tables_dir)
    C = read_table(f"{tables_dir}/IKM_final_master_table_CI.csv").set_index("group")
    Mm = read_table(f"{tables_dir}/IKM_models_panoptosis_layers.csv")
    Pp = read_table(
        f"{tables_dir}/IKM_patel2017_group_overrepresentation.csv"
    ).set_index("group")
    Rr = read_table(f"{tables_dir}/IKM_representation_in_band_units.csv")
    CPd = read_table(f"{tables_dir}/IKM_ht29_rnai_crossplatform.csv")
    Ww = read_table(f"{tables_dir}/IKM_group_weight_all12.csv").set_index("group")
    for name, frame in [("expression", C), ("overlap", Pp), ("magnitude", Ww)]:
        if not frame.index.is_unique or set(frame.index) != set(ORD):
            raise ValueError(f"{name}: exactly one row per plotted group is required")
    for name, frame, column, count in [
        ("CRISPR", Mm, "model", 7),
        ("contrast", Rr, "column", 2),
    ]:
        if frame.duplicated(["group", column]).any() or len(frame) != len(ORD) * count:
            raise ValueError(f"{name}: incomplete or duplicate matrix cells")
        if (
            frame[column].nunique() != count
            or not frame.groupby(column)["group"]
            .apply(lambda s: set(s) == set(ORD))
            .all()
        ):
            raise ValueError(f"{name}: group coverage differs between columns")
        if not np.isfinite(frame["reference_band_units"]).all():
            raise ValueError(f"{name}: non-finite plotted values")
    from pathlib import Path

    Path(out_prefix).parent.mkdir(parents=True, exist_ok=True)
    W = 180.0
    H_A, H_B, GAP, TOP, BOT = (95.0, 108.0, 5.0, 3.0, 5.0)
    H = TOP + H_A + GAP + H_B + BOT
    fig = plt.figure(figsize=(W * MM, H * MM))
    gs = GridSpec(
        2,
        12,
        figure=fig,
        height_ratios=[H_A, H_B],
        hspace=0.4,
        wspace=2.6,
        left=0.075,
        right=0.975,
        top=1 - TOP / H,
        bottom=BOT / H,
    )
    ax_mat = fig.add_subplot(gs[0, 0:12])
    gs_c = gs[1, 5:12].subgridspec(2, 1, height_ratios=[62, 38], hspace=0.06)
    ax_lev = fig.add_subplot(gs_c[0, 0])
    ax_wt = fig.add_subplot(gs_c[1, 0], sharex=ax_lev)
    gs_b = gs[1, 0:3].subgridspec(2, 1, height_ratios=[62, 38], hspace=0.06)
    ax_cp = fig.add_subplot(gs_b[0, 0])
    panel_levels_weight(ax_lev, ax_wt, C, Ww)
    panel_matrix_all(ax_mat, Rr, Mm, Pp)
    panel_hits(ax_cp, CPd)
    for ax, letter, dx, dy in [
        (ax_mat, "a", -0.055, 1.03),
        (ax_cp, "b", -0.36, 1.045),
        (ax_lev, "c", -0.085, 1.03),
    ]:
        ax.text(
            dx,
            dy,
            letter,
            transform=ax.transAxes,
            fontsize=9.5,
            fontweight="bold",
            va="top",
            ha="left",
            color=INK,
        )
    fig.savefig(f"{out_prefix}.png", dpi=400, bbox_inches="tight")
    fig.savefig(f"{out_prefix}.pdf", bbox_inches="tight")
    return fig


if __name__ == "__main__":
    build(*(sys.argv[1:] if len(sys.argv) > 1 else []))
