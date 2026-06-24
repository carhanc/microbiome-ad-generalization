
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.lines import Line2D
import matplotlib.ticker as mticker
from scipy.spatial.distance import pdist, squareform
from scipy.stats import chi2
import os

BASE = "/Users/arhan/Desktop/microbiome-ad-generalization"
TABLES = f"{BASE}/results/tables"
FIGS   = f"{BASE}/results/figures"
DATA   = f"{BASE}/data/processed"

WIDTH_IN = 180 / 25.4
DPI = 300

COHORT_LABELS = {
    "zhuang2018":    "Zhuang 2018\n(China)",
    "ling2020":      "Ling 2020\n(China)",
    "shanghai2022":  "Zhu 2022\n(China)",
    "kazakhstan2022":"Auyezbayeva 2022\n(Kazakhstan)",
    "kbase2022":     "Kim 2022\n(South Korea)",
}
COHORT_SHORT = {
    "zhuang2018":    "Zhuang\n(CN)",
    "ling2020":      "Ling\n(CN)",
    "shanghai2022":  "Zhu 2022\n(CN)",
    "kazakhstan2022":"Auyezbayeva\n(KZ)",
}
COHORT_SHORT_FLAT = {
    "zhuang2018":    "Zhuang",
    "ling2020":      "Ling",
    "shanghai2022":  "Zhu 2022",
    "kazakhstan2022":"Auyezbayeva",
}
ALL5_SHORT = {
    "zhuang2018":    "Zhuang",
    "ling2020":      "Ling",
    "shanghai2022":  "Zhu 2022",
    "kbase2022":     "Kim",
    "kazakhstan2022":"Auyezbayeva",
}
COLORS = {
    "zhuang2018":    "#2166ac",
    "ling2020":      "#4dac26",
    "shanghai2022":  "#d01c8b",
    "kazakhstan2022":"#f4a582",
    "kbase2022":     "#b35806",
}
MODEL_COLORS = {"logreg": "#2166ac", "lgbm": "#d01c8b"}
CORRECTION_COLORS = {
    "uncorrected": "#4dac26",
    "combatseq":   "#d01c8b",
    "mmuphin":     "#f4a582",
}

LABELED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]

def make_fig1():
    print("  Generating Figure 1 (Study Design)…")
    fig = plt.figure(figsize=(WIDTH_IN, WIDTH_IN * 0.88), dpi=DPI)
    fig.patch.set_facecolor("white")

    gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=[1.05, 0.95],
                           wspace=0.06, left=0.01, right=0.99,
                           top=0.92, bottom=0.04)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_a.axis("off")
    ax_b.axis("off")

    ax_a.text(0.0, 1.02, "A", transform=ax_a.transAxes, fontsize=11,
              fontweight="bold", va="bottom")
    ax_b.text(0.0, 1.02, "B", transform=ax_b.transAxes, fontsize=11,
              fontweight="bold", va="bottom")

    ax = ax_a

    country_colors = {
        "China":      "#d73027",
        "Kazakhstan": "#c8a400",
        "S. Korea":   "#4575b4",
    }

    col_headers = ["Cohort", "Country", "N (labels)", "Seq.",
                   "Region", "Diagnosis", "Analysis"]

    rows = [
        ["Zhuang 2018\nPRJNA554111",      "China",      "86\n43AD/43CN",
         "PE MiSeq\n2×300 bp",  "V3–V4", "AD / CN",        "✓ Ph. 2–6"],
        ["Ling 2020\nPRJNA633959",        "China",      "171\n100AD/71CN",
         "PE MiSeq\n2×300 bp",  "V3–V4", "AD / CN",        "✓ Ph. 2–6"],
        ["Zhu 2022\nPRJNA489760",          "China",      "90†\n(30 each)",
         "PE MiSeq\n2×300 bp",  "V3–V4", "AD/MCI/CN",      "✓ Ph. 2–6"],
        ["Auyezbayeva 2022\nPRJNA811324", "Kazakhstan", "84\n41AD/43CN",
         "SE NovaSeq",           "V3–V4", "AD / CN",        "✓ Ph. 2–6"],
        ["Kim 2022\nPRJEB50447",          "S. Korea",   "78\n18pcAD/60CN",
         "PE MiSeq\n2×300 bp",  "V3–V4", "IRB-restricted", "Ph. 4 only"],
    ]

    row_countries = ["China", "China", "China", "Kazakhstan", "S. Korea"]
    row_colors    = ["#f7fbff", "#deebf7", "#c6dbef", "#9ecae1", "#d8e9f8"]

    col_w = [0.18, 0.13, 0.13, 0.12, 0.09, 0.19, 0.16]  # sum = 1.00

    n_rows = len(rows)
    n_cols = len(col_headers)

    the_table = ax.table(
        cellText=rows,
        colLabels=col_headers,
        colWidths=col_w,
        bbox=[0.0, 0.28, 1.0, 0.65],
        cellLoc="center",
    )
    the_table.auto_set_font_size(False)
    the_table.set_fontsize(5.0)

    for j in range(n_cols):
        cell = the_table[0, j]
        cell.set_facecolor("#2c3e50")
        cell.set_text_props(color="white", fontweight="bold", fontsize=5.2)
        cell.set_edgecolor("#2c3e50")

    for i in range(n_rows):
        country    = row_countries[i]
        cdot_color = country_colors.get(country, "gray")
        row_bg     = row_colors[i % len(row_colors)]

        for j in range(n_cols):
            cell = the_table[i + 1, j]   # +1 because row 0 is header
            cell.set_facecolor(row_bg)
            cell.set_edgecolor("#cccccc")
            cell.set_linewidth(0.4)

            if j == 0:
                cell.set_text_props(fontweight="bold", fontsize=4.5,
                                    color="#111111")
            elif j == 1:
                cell.set_text_props(color=cdot_color, fontweight="normal",
                                    fontsize=5.0)
            elif j == n_cols - 1:
                txt_color = "#27ae60" if "✓" in rows[i][j] else "#c0392b"
                cell.set_text_props(color=txt_color, fontweight="bold",
                                    fontsize=5.0)
            else:
                cell.set_text_props(fontsize=5.0, color="#111111")

    legend_y = 0.18
    ax.text(0.0, legend_y, "Country:", transform=ax.transAxes,
            fontsize=5.5, va="center", fontweight="bold", color="#333333")
    offset_x = 0.14
    for cname, cc in country_colors.items():
        ax.add_patch(plt.Circle((offset_x, legend_y), 0.010, color=cc,
                                transform=ax.transAxes, clip_on=False, zorder=4))
        ax.text(offset_x + 0.018, legend_y, cname,
                transform=ax.transAxes, fontsize=5.0, va="center", color="#333333")
        offset_x += 0.28
    ax.text(0.0, 0.08, "pcAD = preclinical AD (amyloid-PET+, cognitively normal); "
            "Seq. = sequencing; SE = single-end; PE = paired-end",
            transform=ax.transAxes, fontsize=4.2, va="center", color="#666666",
            style="italic")

    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 1.04, "Cohort Overview", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=9, fontweight="bold")

    ax = ax_b
    steps = [
        ("Phase 1", "DADA2 Processing\n& Genus Harmonization\n(509 samples, 396 genera)", "#1a6faf"),
        ("Phase 2", "Within-Cohort Baseline\nNested 10×5 CV\n(LogReg + LGBM)", "#2ca02c"),
        ("Phase 3", "Cross-Cohort\nGeneralization\n(LOCO + Pairwise)", "#d62728"),
        ("Phase 4", "PERMANOVA Variance\nDecomposition\n(Aitchison + Bray-Curtis)", "#8c564b"),
        ("Phase 5", "Batch Correction\n(ComBat-seq + MMUPHin)\n& LOCO Re-test", "#9467bd"),
        ("Phase 6", "SHAP Feature\nImportance\n(Within-cohort + LOCO)", "#e377c2"),
    ]
    box_h = 0.105
    gap   = 0.035
    start_y = 0.92
    box_w = 0.78
    box_x = 0.11

    y_pos = start_y
    for i, (phase, desc, color) in enumerate(steps):
        arr_src_y = y_pos + box_h + gap * 0.5   # midpoint of gap above this box
        if i == 0:
            arr_src_y = y_pos + box_h + gap      # just above Phase 1
        ax.annotate("", xy=(box_x + box_w/2, y_pos + box_h),
                    xytext=(box_x + box_w/2, arr_src_y),
                    xycoords="axes fraction", textcoords="axes fraction",
                    arrowprops=dict(arrowstyle="-|>", color="#555555",
                                   lw=1.0, mutation_scale=8))
        ax.add_patch(FancyBboxPatch((box_x, y_pos), box_w, box_h,
                                    boxstyle="round,pad=0.015",
                                    fc=color, ec="none",
                                    transform=ax.transAxes, clip_on=False,
                                    alpha=0.88))
        ax.text(box_x + 0.03, y_pos + box_h / 2, phase,
                transform=ax.transAxes, ha="left", va="center",
                fontsize=6.5, fontweight="bold", color="white")
        ax.text(box_x + 0.22, y_pos + box_h / 2, desc,
                transform=ax.transAxes, ha="left", va="center",
                fontsize=6.5, color="white")
        y_pos -= (box_h + gap)

    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.text(0.5, 1.04, "Analysis Pipeline", transform=ax.transAxes,
            ha="center", va="bottom", fontsize=9, fontweight="bold")

    out = f"{FIGS}/manuscript_fig1.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


def make_fig2():
    print("  Generating Figure 2 (Within-Cohort Baseline AUC)…")
    df = pd.read_csv(f"{TABLES}/within_cohort_auc.csv")
    df = df[df["cohort"].isin(LABELED_COHORTS)]

    fig, ax = plt.subplots(figsize=(WIDTH_IN * 0.7, WIDTH_IN * 0.5), dpi=DPI)
    fig.patch.set_facecolor("white")

    cohorts = LABELED_COHORTS
    n = len(cohorts)
    x = np.arange(n)
    bar_w = 0.32
    models = ["logreg", "lgbm"]
    model_labels = {"logreg": "Logistic Regression", "lgbm": "LightGBM"}
    offsets = [-0.18, 0.18]

    for m, offset in zip(models, offsets):
        aucs, ci_lo, ci_hi = [], [], []
        for c in cohorts:
            row = df[(df["cohort"] == c) & (df["model"] == m)].iloc[0]
            aucs.append(row["auc_oof"])
            ci_lo.append(row["auc_oof"] - row["auc_ci_lo"])
            ci_hi.append(row["auc_ci_hi"] - row["auc_oof"])
        bars = ax.bar(x + offset, aucs, bar_w,
                      color=MODEL_COLORS[m], alpha=0.85,
                      label=model_labels[m], zorder=3)
        ax.errorbar(x + offset, aucs, yerr=[ci_lo, ci_hi],
                    fmt="none", color="#222222", capsize=3,
                    capthick=0.8, elinewidth=0.8, zorder=4)
        for xi, auc, hi in zip(x + offset, aucs, ci_hi):
            ax.text(xi, auc + hi + 0.015, f"{auc:.3f}",
                    ha="center", va="bottom", fontsize=6.5, color="#333333")

    ax.axhline(0.5, color="gray", linestyle="--", lw=0.8, label="Chance (AUC=0.5)", zorder=2)
    ax.set_xticks(x)
    two_line_labels = {
        "zhuang2018":    "Zhuang\n2018",
        "ling2020":      "Ling\n2020",
        "shanghai2022":  "Zhu\n2022",
        "kazakhstan2022":"Auyezbayeva\n2022",
    }
    ax.set_xticklabels([two_line_labels[c] for c in cohorts], fontsize=8)
    ax.set_ylabel("AUC-ROC", fontsize=9)
    ax.set_xlabel("Cohort", fontsize=9)
    ax.set_title("Within-Cohort Classification Performance\n(Nested 10-fold outer / 5-fold inner CV)", fontsize=9)
    ax.set_ylim(0.35, 1.15)   # raised from 1.12 so 0.998/0.992 value labels clear the top
    ax.yaxis.set_major_formatter(mticker.FormatStrFormatter("%.2f"))
    ax.legend(fontsize=8, loc="upper left", framealpha=0.9)
    ax.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6, zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    plt.tight_layout()

    out = f"{FIGS}/manuscript_fig2.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


def make_fig3():
    print("  Generating Figure 3 (Cross-Cohort Generalization)…")
    drop_df  = pd.read_csv(f"{TABLES}/auc_drop_summary.csv")
    loco_df  = pd.read_csv(f"{TABLES}/loco_auc.csv")
    pair_df  = pd.read_csv(f"{TABLES}/pairwise_auc.csv")

    fig = plt.figure(figsize=(WIDTH_IN, WIDTH_IN * 0.68), dpi=DPI)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[1.1, 1.0, 1.0],
                           wspace=0.35, left=0.07, right=0.97,
                           top=0.88, bottom=0.18)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])

    for ax, lbl in zip([ax_a, ax_b, ax_c], "ABC"):
        ax.text(-0.15, 1.08, lbl, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom")

    cohorts = LABELED_COHORTS
    n = len(cohorts)
    x = np.arange(n)
    bar_w = 0.15
    models = ["logreg", "lgbm"]
    model_labels = {"logreg": "LogReg", "lgbm": "LGBM"}
    offsets_within = [-0.225, -0.075]
    offsets_loco   = [ 0.075,  0.225]

    within_df = pd.read_csv(f"{TABLES}/within_cohort_auc.csv")

    for m, ow, ol in zip(models, offsets_within, offsets_loco):
        w_aucs, l_aucs, l_lo, l_hi = [], [], [], []
        for c in cohorts:
            wr = within_df[(within_df["cohort"] == c) & (within_df["model"] == m)].iloc[0]
            lr = loco_df[(loco_df["test_cohort"] == c) & (loco_df["model"] == m)].iloc[0]
            w_aucs.append(wr["auc_oof"])
            l_aucs.append(lr["auc"])
            l_lo.append(lr["auc"] - lr["ci_lo"])
            l_hi.append(lr["ci_hi"] - lr["auc"])

        col = MODEL_COLORS[m]
        ax_a.bar(x + ow, w_aucs, bar_w, color=col, alpha=0.9,
                 label=f"{model_labels[m]} within", zorder=3)
        bars = ax_a.bar(x + ol, l_aucs, bar_w, color=col, alpha=0.45,
                        hatch="///", label=f"{model_labels[m]} LOCO", zorder=3,
                        edgecolor=col, linewidth=0.5)
        ax_a.errorbar(x + ol, l_aucs, yerr=[l_lo, l_hi],
                      fmt="none", color="#222222", capsize=2.5,
                      capthick=0.7, elinewidth=0.7, zorder=4)

    ax_a.axhline(0.5, color="gray", linestyle="--", lw=0.8, zorder=2)
    ax_a.set_xticks(x)
    ax_a.set_xticklabels([COHORT_SHORT_FLAT[c] for c in cohorts], fontsize=7.5,
                          rotation=20, ha="right")
    ax_a.set_ylabel("AUC-ROC", fontsize=8.5)
    ax_a.set_title("Within-cohort vs. LOCO\nCross-cohort AUC", fontsize=8.5)
    ax_a.set_ylim(0.30, 1.12)
    ax_a.legend(fontsize=6.0, ncol=2, loc="upper right", framealpha=0.85,
                columnspacing=0.5, handlelength=1.0)
    ax_a.grid(axis="y", linestyle=":", linewidth=0.5, alpha=0.6, zorder=1)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)

    import matplotlib.colors as mcolors
    for ax, model, panel_lbl in [(ax_b, "logreg", "Logistic Regression"),
                                  (ax_c, "lgbm",   "LightGBM")]:
        sub = pair_df[pair_df["model"] == model]
        train_cs = cohorts
        test_cs  = cohorts
        mat = np.full((len(train_cs), len(test_cs)), np.nan)
        for _, row in sub.iterrows():
            ti = train_cs.index(row["train_cohort"])
            vi = test_cs.index(row["test_cohort"])
            mat[ti, vi] = row["auc"]

        cmap = plt.cm.RdYlGn
        norm = mcolors.Normalize(vmin=0.45, vmax=0.85)
        im = ax.imshow(mat, cmap=cmap, norm=norm, aspect="auto")

        for i in range(len(train_cs)):
            for j in range(len(test_cs)):
                val = mat[i, j]
                if not np.isnan(val):
                    tc = "black" if 0.55 < val < 0.75 else "white"
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center",
                            fontsize=7.5, color=tc, fontweight="bold")

        ax.set_xticks(range(len(test_cs)))
        ax.set_yticks(range(len(train_cs)))
        xlabels = [COHORT_SHORT_FLAT[c] for c in test_cs]
        ylabels = [COHORT_SHORT_FLAT[c] for c in train_cs]
        ax.set_xticklabels(xlabels, fontsize=7, rotation=30, ha="right")
        ax.set_yticklabels(ylabels, fontsize=7)
        ax.set_xlabel("Test Cohort", fontsize=8)
        ax.set_ylabel("Train Cohort", fontsize=8)
        ax.set_title(f"Pairwise AUC\n({panel_lbl})", fontsize=8.5)
        cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, shrink=0.82)
        cb.ax.tick_params(labelsize=7)
        cb.set_label("AUC", fontsize=7.5)

    plt.tight_layout()
    out = f"{FIGS}/manuscript_fig3.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


def make_fig4():
    print("  Generating Figure 4 (PERMANOVA Variance Decomposition)…")

    clr_df = pd.read_csv(f"{DATA}/clr_matrix.csv")
    sample_ids  = clr_df["sample_id"].values
    cohort_col  = clr_df["cohort"].values
    genus_cols  = [c for c in clr_df.columns if c not in ("sample_id", "cohort")]
    clr_vals    = clr_df[genus_cols].values.astype(float)

    print(f"    Computing Aitchison distance matrix ({len(clr_vals)}×{len(clr_vals)})…")
    # 509×509 pairwise — grab coffee
    dist_mat = squareform(pdist(clr_vals, metric="euclidean"))
    n = dist_mat.shape[0]
    D2 = dist_mat ** 2
    H  = np.eye(n) - np.ones((n, n)) / n
    B  = -0.5 * H @ D2 @ H
    eigvals, eigvecs = np.linalg.eigh(B)
    idx = np.argsort(eigvals)[::-1]
    eigvals = eigvals[idx]; eigvecs = eigvecs[:, idx]
    pos_mask = eigvals > 0
    coords = eigvecs[:, pos_mask] * np.sqrt(eigvals[pos_mask])
    pct_var = eigvals[pos_mask] / eigvals[pos_mask].sum() * 100

    print("    Computing betadisper distances…")
    unique_cohorts = list(dict.fromkeys(cohort_col))  # preserve order
    bd_data = []
    for c in unique_cohorts:
        mask = cohort_col == c
        pts  = clr_vals[mask]
        centroid = pts.mean(axis=0)
        dists = np.linalg.norm(pts - centroid, axis=1)
        for d in dists:
            bd_data.append({"cohort": c, "distance": d})
    bd_df = pd.DataFrame(bd_data)

    fig = plt.figure(figsize=(WIDTH_IN, WIDTH_IN * 0.52), dpi=DPI)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(1, 3, figure=fig, width_ratios=[0.75, 1.3, 1.0],
                           wspace=0.38, left=0.05, right=0.97,
                           top=0.87, bottom=0.14)
    ax_a = fig.add_subplot(gs[0])
    ax_b = fig.add_subplot(gs[1])
    ax_c = fig.add_subplot(gs[2])
    for ax, lbl in zip([ax_a, ax_b, ax_c], "ABC"):
        ax.text(-0.18, 1.08, lbl, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom")

    r2_cohort = 0.1722
    r2_diag   = 0.0140
    r2_resid  = 1 - r2_cohort - r2_diag
    bar_colors = ["#d73027", "#4575b4", "#d9d9d9"]
    labels_r2  = [f"Cohort (R²={r2_cohort:.3f})",
                  f"Diagnosis (R²={r2_diag:.3f})",
                  f"Residual (R²={r2_resid:.3f})"]
    vals = [r2_cohort, r2_diag, r2_resid]
    bot = 0
    for v, c, lbl in zip(vals, bar_colors, labels_r2):
        ax_a.bar(0, v, 0.5, bottom=bot, color=c, label=lbl, alpha=0.9)
        if v > 0.025:
            ax_a.text(0, bot + v / 2, f"{v:.3f}", ha="center", va="center",
                      fontsize=8, fontweight="bold", color="white")
        bot += v
    ax_a.set_xlim(-0.5, 0.5)
    ax_a.set_ylim(0, 1.08)
    ax_a.set_xticks([])
    ax_a.set_ylabel("Proportion of Aitchison Variance (R²)", fontsize=8)
    ax_a.set_title("PERMANOVA\nVariance Partitioning\n(4 labeled cohorts)", fontsize=8.5)
    ax_a.legend(fontsize=7, loc="upper right", framealpha=0.9,
                bbox_to_anchor=(2.85, 1.02))
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.spines["bottom"].set_visible(False)

    pc1 = coords[:, 0]
    pc2 = coords[:, 1]
    pv1 = pct_var[0]; pv2 = pct_var[1]

    for c in unique_cohorts:
        mask = cohort_col == c
        lbl  = COHORT_LABELS.get(c, c)
        ax_b.scatter(pc1[mask], pc2[mask], s=6, alpha=0.65,
                     color=COLORS.get(c, "gray"),
                     label=lbl.replace("\n", " "), zorder=3, linewidths=0)
        pts_2d = np.column_stack([pc1[mask], pc2[mask]])
        if len(pts_2d) >= 3:
            mu  = pts_2d.mean(axis=0)
            cov = np.cov(pts_2d.T)
            vals_e, vecs_e = np.linalg.eigh(cov)
            order_e = vals_e.argsort()[::-1]
            vals_e  = vals_e[order_e]; vecs_e = vecs_e[:, order_e]
            chi2_val = chi2.ppf(0.95, df=2)
            w = 2 * np.sqrt(chi2_val * vals_e[0])
            h = 2 * np.sqrt(chi2_val * vals_e[1])
            angle = np.degrees(np.arctan2(*vecs_e[:, 0][::-1]))
            from matplotlib.patches import Ellipse
            ell = Ellipse(xy=mu, width=w, height=h, angle=angle,
                          fc="none", ec=COLORS.get(c, "gray"),
                          lw=0.9, alpha=0.75, zorder=2)
            ax_b.add_patch(ell)

    ax_b.set_xlabel(f"PC1 ({pv1:.1f}%)", fontsize=8.5)
    ax_b.set_ylabel(f"PC2 ({pv2:.1f}%)", fontsize=8.5)
    ax_b.set_title(f"PCoA — Aitchison Distance\nCohort R²=0.193 (p<0.0001)", fontsize=8.5)
    ax_b.legend(fontsize=6.5, loc="upper right", framealpha=0.85,
                markerscale=1.5)
    ax_b.axhline(0, color="gray", lw=0.4, ls=":"); ax_b.axvline(0, color="gray", lw=0.4, ls=":")
    ax_b.spines["top"].set_visible(False); ax_b.spines["right"].set_visible(False)

    plot_cohorts = [c for c in unique_cohorts]  # all 5
    positions = list(range(len(plot_cohorts)))
    for pos, c in zip(positions, plot_cohorts):
        vals_bd = bd_df[bd_df["cohort"] == c]["distance"].values
        bp = ax_c.boxplot(vals_bd, positions=[pos], widths=0.55,
                          patch_artist=True, notch=False,
                          medianprops=dict(color="black", lw=1.5),
                          boxprops=dict(facecolor=COLORS.get(c, "gray"), alpha=0.75),
                          whiskerprops=dict(lw=0.8),
                          capprops=dict(lw=0.8),
                          flierprops=dict(marker="o", markersize=2,
                                          markerfacecolor=COLORS.get(c, "gray"),
                                          alpha=0.5))
    ax_c.set_xticks(positions)
    ax_c.set_xticklabels([ALL5_SHORT.get(c, c) for c in plot_cohorts],
                          fontsize=7.5, rotation=20, ha="right")
    ax_c.set_ylabel("Distance to Cohort Centroid\n(Aitchison)", fontsize=8)
    ax_c.set_title(f"Beta-Dispersion by Cohort\n(F=13.93, p<0.0001)", fontsize=8.5)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.grid(axis="y", linestyle=":", lw=0.5, alpha=0.6)

    plt.tight_layout()
    out = f"{FIGS}/manuscript_fig4.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


def make_fig5():
    print("  Generating Figure 5 (Batch Correction Comparison)…")
    df = pd.read_csv(f"{TABLES}/auc_comparison_table.csv")

    fig, axes = plt.subplots(1, 2, figsize=(WIDTH_IN, WIDTH_IN * 0.46), dpi=DPI)
    fig.patch.set_facecolor("white")

    corrections = ["uncorrected", "combatseq", "mmuphin"]
    corr_labels = {"uncorrected": "Uncorrected", "combatseq": "ComBat-seq",
                   "mmuphin": "MMUPHin"}
    cohorts = LABELED_COHORTS
    n = len(cohorts)
    x = np.arange(n)
    bar_w = 0.24

    for ax, model, model_lbl in zip(axes, ["logreg", "lgbm"],
                                     ["Logistic Regression", "LightGBM"]):
        offsets = np.array([-0.25, 0.0, 0.25])
        for corr, offset in zip(corrections, offsets):
            aucs, ci_lo, ci_hi = [], [], []
            for c in cohorts:
                row = df[(df["test_cohort"] == c) & (df["model"] == model) &
                         (df["correction"] == corr)]
                if len(row) == 0:
                    aucs.append(np.nan); ci_lo.append(0); ci_hi.append(0)
                else:
                    row = row.iloc[0]
                    aucs.append(row["loco_auc"])
                    ci_lo.append(row["loco_auc"] - row["loco_ci_lo"])
                    ci_hi.append(row["loco_ci_hi"] - row["loco_auc"])
            ax.bar(x + offset, aucs, bar_w, color=CORRECTION_COLORS[corr],
                   alpha=0.85, label=corr_labels[corr], zorder=3)
            ax.errorbar(x + offset, aucs, yerr=[ci_lo, ci_hi],
                        fmt="none", color="#222222", capsize=2.5,
                        capthick=0.7, elinewidth=0.7, zorder=4)

        ax.axhline(0.5, color="gray", linestyle="--", lw=0.8, zorder=2,
                   label="Chance (0.5)")
        ax.set_xticks(x)
        ax.set_xticklabels([COHORT_SHORT_FLAT[c] for c in cohorts],
                            fontsize=8, rotation=15, ha="right")
        ax.set_ylabel("LOCO AUC-ROC", fontsize=8.5)
        ax.set_title(f"LOCO AUC: Batch Correction Effect\n({model_lbl})", fontsize=8.5)
        ax.set_ylim(0.15, 1.05)
        ax.legend(fontsize=7.5, loc="upper right", framealpha=0.9)
        ax.grid(axis="y", linestyle=":", lw=0.5, alpha=0.6, zorder=1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.text(-0.12, 1.08, "AB"[list(axes).index(ax)],
                transform=ax.transAxes, fontsize=11, fontweight="bold",
                va="bottom")

    plt.tight_layout()
    out = f"{FIGS}/manuscript_fig5.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


def make_fig6():
    print("  Generating Figure 6 (SHAP Taxa Analysis)…")
    shap_df  = pd.read_csv(f"{TABLES}/shap_within_cohort_importance.csv")
    flip_df  = pd.read_csv(f"{TABLES}/shap_directional_flips.csv")
    over_df  = pd.read_csv(f"{TABLES}/shap_taxa_overlap.csv")

    fig_h = WIDTH_IN * 0.62
    fig = plt.figure(figsize=(WIDTH_IN, fig_h), dpi=DPI)
    fig.patch.set_facecolor("white")
    gs = gridspec.GridSpec(2, 2, figure=fig,
                           height_ratios=[1.15, 1.0],
                           hspace=0.52, wspace=0.40,
                           left=0.15, right=0.97,
                           top=0.93, bottom=0.08)
    ax_a = fig.add_subplot(gs[0, :])   # top full-width
    ax_b = fig.add_subplot(gs[1, 0])   # bottom-left
    ax_c = fig.add_subplot(gs[1, 1])   # bottom-right

    for ax, lbl in zip([ax_a, ax_b, ax_c], "ABC"):
        x_off = -0.04 if ax is ax_a else -0.08
        ax.text(x_off, 1.05, lbl, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom")

    AD_COL = "#d62728"
    CN_COL = "#1f77b4"

    TOP_PER_COHORT = 15   # taxa pool: top-15 from each cohort
    TOP_DISPLAY    = 20   # rows shown in the dot-plot

    sub_all = shap_df[(shap_df["model"] == "logreg") &
                      (shap_df["cohort"].isin(LABELED_COHORTS)) &
                      (shap_df["rank"] <= TOP_PER_COHORT)]

    taxa_max = (sub_all.groupby("taxon")["mean_abs_shap"]
                       .max()
                       .sort_values(ascending=False))
    display_taxa = taxa_max.head(TOP_DISPLAY).index.tolist()

    sub = sub_all[sub_all["taxon"].isin(display_taxa)]

    cohort_order = LABELED_COHORTS
    cohort_short = [COHORT_SHORT_FLAT[c] for c in cohort_order]
    n_coh = len(cohort_order)

    taxa_order_plot = display_taxa  # already sorted descending by max |SHAP|
    n_taxa = len(taxa_order_plot)
    y_pos  = {t: (n_taxa - 1 - i) for i, t in enumerate(taxa_order_plot)}
    x_pos  = {c: j for j, c in enumerate(cohort_order)}

    flip_taxa_lr = set(flip_df[flip_df["model"] == "logreg"]["taxon"].unique())

    for _, row in sub.iterrows():
        t = row["taxon"]
        c = row["cohort"]
        if t not in y_pos:
            continue
        yy = y_pos[t]
        xx = x_pos[c]
        size  = max(25, row["mean_abs_shap"] * 220)
        color = AD_COL if row["mean_shap_AD"] > row["mean_shap_CN"] else CN_COL
        ax_a.scatter(xx, yy, s=size, color=color, alpha=0.82,
                     zorder=3, linewidths=0.3, edgecolors="white")
        ax_a.text(xx, yy, f"{row['mean_abs_shap']:.2f}",
                  ha="center", va="center", fontsize=5.0,
                  color="white", fontweight="bold", zorder=4)

    for t, yy in y_pos.items():
        if t in flip_taxa_lr:
            ax_a.axhspan(yy - 0.46, yy + 0.46, facecolor="#fff3b0",
                         alpha=0.55, zorder=1)

    ax_a.set_xlim(-0.55, n_coh - 0.45)
    ax_a.set_ylim(-0.55, n_taxa - 0.45)
    ax_a.set_xticks(range(n_coh))
    ax_a.set_xticklabels(cohort_short, fontsize=9)
    ax_a.xaxis.set_ticks_position("top")
    ax_a.xaxis.set_label_position("top")
    ax_a.set_xlabel("Cohort", fontsize=9, labelpad=4)

    ylabels = [f"★ {t}" if t in flip_taxa_lr else f"   {t}"
               for t in sorted(y_pos, key=lambda t: y_pos[t], reverse=True)]
    ax_a.set_yticks(sorted(y_pos.values(), reverse=True))
    ax_a.set_yticklabels(ylabels, fontsize=7.5)
    ax_a.set_ylabel("Genus (top-20 by max |SHAP|)", fontsize=8.5)
    ax_a.set_title(
        "Top-20 SHAP Taxa — Logistic Regression  "
        "(dot size ∝ |SHAP|; red = AD↑, blue = CN↑; ★ = directional flip; "
        "yellow = flip rows)",
        fontsize=8, pad=4)
    ax_a.grid(axis="both", linestyle=":", lw=0.4, alpha=0.45, zorder=0)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.spines["bottom"].set_visible(False)
    legend_elements = [
        Line2D([0], [0], marker="o", color="w", markerfacecolor=AD_COL,
               markersize=8, label="AD-associated"),
        Line2D([0], [0], marker="o", color="w", markerfacecolor=CN_COL,
               markersize=8, label="CN-associated"),
        mpatches.Patch(facecolor="#fff3b0", ec="#ccaa00",
                       label="Directional flip"),
    ]
    ax_a.legend(handles=legend_elements, fontsize=7.5, loc="lower right",
                framealpha=0.92, ncol=3, columnspacing=0.8)

    flip_taxa = flip_df[(flip_df["model"] == "logreg") &
                        flip_df["in_top_n"]]["taxon"].unique()
    flip_sub = shap_df[(shap_df["model"] == "logreg") &
                       (shap_df["taxon"].isin(flip_taxa)) &
                       (shap_df["cohort"].isin(LABELED_COHORTS))]

    flip_sub = flip_sub.sort_values("rank").drop_duplicates(["taxon","cohort"])

    taxa_flip_order = (flip_sub.groupby("taxon")["mean_abs_shap"].max()
                               .sort_values(ascending=False).index.tolist())

    for _, row in flip_sub.iterrows():
        t = row["taxon"]
        c = row["cohort"]
        if t not in taxa_flip_order: continue
        yy = taxa_flip_order.index(t)
        shap_v = row["mean_shap_AD"] - row["mean_shap_CN"]
        size = max(30, abs(row["mean_abs_shap"]) * 180)
        color = AD_COL if shap_v > 0 else CN_COL
        ax_b.scatter(shap_v, yy, s=size, color=color, alpha=0.8,
                     zorder=3, linewidths=0.4, edgecolors="white")
        ax_b.text(shap_v, yy, COHORT_SHORT_FLAT[c][0],
                  ha="center", va="center", fontsize=5,
                  color="white", fontweight="bold", zorder=4)

    ax_b.axvline(0, color="black", lw=0.8, zorder=2)
    ax_b.set_yticks(range(len(taxa_flip_order)))
    ax_b.set_yticklabels([f"★ {t}" for t in taxa_flip_order], fontsize=7.5)
    ax_b.set_xlabel("Mean SHAP (AD-direction → positive)", fontsize=8)
    ax_b.set_title("Directional Flip Taxa\n(LogReg; letter = cohort initial)", fontsize=8.5)
    ax_b.set_xlim(ax_b.get_xlim()[0] * 1.15, ax_b.get_xlim()[1] * 1.15)
    ax_b.grid(axis="x", linestyle=":", lw=0.5, alpha=0.6, zorder=0)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.text(0.04, 0.98, "← CN↑ | AD↑ →", transform=ax_b.transAxes,
              ha="left", va="top", fontsize=6.5, color="#555555")

    sub_over = over_df[(over_df["model"] == "logreg") & (over_df["top_n"] == 20)]
    cohorts4 = LABELED_COHORTS
    jac_mat = np.zeros((4, 4))
    for _, row in sub_over.iterrows():
        ca = row.get("cohort_a", row.get("cohort_i", None))
        cb = row.get("cohort_b", row.get("cohort_j", None))
        ci = cohorts4.index(ca) if ca in cohorts4 else -1
        cj = cohorts4.index(cb) if cb in cohorts4 else -1
        if ci >= 0 and cj >= 0:
            jac_mat[ci, cj] = row["jaccard"]
            jac_mat[cj, ci] = row["jaccard"]
    np.fill_diagonal(jac_mat, 1.0)

    im = ax_c.imshow(jac_mat, cmap="YlOrRd", vmin=0, vmax=0.25, aspect="auto")
    for i in range(4):
        for j in range(4):
            v = jac_mat[i, j]
            tc = "white" if v > 0.18 else "black"
            ax_c.text(j, i, f"{v:.3f}", ha="center", va="center",
                      fontsize=8, color=tc, fontweight="bold")

    short_labels = [COHORT_SHORT_FLAT[c] for c in cohorts4]
    ax_c.set_xticks(range(4)); ax_c.set_yticks(range(4))
    ax_c.set_xticklabels(short_labels, fontsize=8, rotation=20, ha="right")
    ax_c.set_yticklabels(short_labels, fontsize=8)
    ax_c.set_title("Pairwise Jaccard Similarity\n(Top-20 SHAP taxa, LogReg)", fontsize=8.5)
    cb = fig.colorbar(im, ax=ax_c, fraction=0.046, pad=0.04, shrink=0.85)
    cb.set_label("Jaccard Index", fontsize=7.5)
    cb.ax.tick_params(labelsize=7)

    plt.tight_layout()
    out = f"{FIGS}/manuscript_fig6.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


def make_supp_fig_s1():
    print("  Generating Figure S1 (Zhu 2022 fecal-only AUC)…")

    WIDTH_S = 85 / 25.4

    wc = pd.read_csv(f"{TABLES}/within_cohort_auc.csv")
    sh = wc[wc["cohort"] == "shanghai2022"].set_index("model")

    lr_auc   = sh.loc["logreg", "auc_oof"]
    lr_lo    = sh.loc["logreg", "auc_ci_lo"]
    lr_hi    = sh.loc["logreg", "auc_ci_hi"]
    lgbm_auc = sh.loc["lgbm",  "auc_oof"]
    lgbm_lo  = sh.loc["lgbm",  "auc_ci_lo"]
    lgbm_hi  = sh.loc["lgbm",  "auc_ci_hi"]

    models   = ["LogReg", "LGBM"]
    aucs     = [lr_auc,   lgbm_auc]
    errs_lo  = [lr_auc  - lr_lo,   lgbm_auc - lgbm_lo]
    errs_hi  = [lr_hi   - lr_auc,  lgbm_hi  - lgbm_auc]
    bar_col  = "#2166ac"

    x = np.array([0, 1])

    fig, ax = plt.subplots(figsize=(WIDTH_S, WIDTH_S * 1.15), dpi=DPI)
    fig.patch.set_facecolor("white")

    ax.bar(x, aucs, 0.5, color=bar_col, alpha=0.88, zorder=3)
    ax.errorbar(x, aucs, yerr=[errs_lo, errs_hi],
                fmt="none", color="#333333",
                capsize=3, capthick=0.8, elinewidth=0.8, zorder=4)
    for xi, auc, hi in zip(x, aucs, errs_hi):
        ax.text(xi, auc + hi + 0.004, f"{auc:.3f}",
                ha="center", va="bottom", fontsize=7.5,
                color="#111111", fontweight="bold")

    ax.axhline(0.5, color="gray", ls="--", lw=0.7, zorder=2)
    ax.set_xticks(x)
    ax.set_xticklabels(models, fontsize=9)
    ax.set_ylabel("AUC-ROC (within-cohort)", fontsize=8.5)
    ax.set_ylim(0.88, 1.055)
    ax.set_title("Zhu 2022: Fecal-Only Within-Cohort AUC\n"
                 "(StratifiedKFold; n=60 binary: 30 AD + 30 CN)", fontsize=8)
    ax.grid(axis="y", linestyle=":", lw=0.5, alpha=0.6, zorder=1)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    plt.tight_layout()
    out = f"{FIGS}/supp_fig_s1.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


def make_supp_fig_s2():
    print("  Generating Figure S2 (top-15 SHAP taxa, LightGBM)…")

    imp = pd.read_csv(f"{TABLES}/shap_within_cohort_importance.csv")
    lgbm = imp[imp["model"] == "lgbm"].copy()

    COHORT_LABELS_S2 = {
        "zhuang2018":    "Zhuang 2018\n(China, n=86)",
        "ling2020":      "Ling 2020\n(China, n=171)",
        "shanghai2022":  "Zhu 2022\n(China, n=60)†",
        "kazakhstan2022":"Auyezbayeva 2022\n(KZ, n=84)",
    }

    cohorts = LABELED_COHORTS
    fig, axes = plt.subplots(2, 2, figsize=(WIDTH_IN, WIDTH_IN * 0.9),
                              dpi=DPI)
    fig.patch.set_facecolor("white")

    for ax, cohort in zip(axes.flat, cohorts):
        sub = (lgbm[lgbm["cohort"] == cohort]
               .sort_values("mean_abs_shap", ascending=False)
               .head(15)
               .sort_values("mean_abs_shap", ascending=True))

        colours = ["#d73027" if v >= 0 else "#4575b4" for v in sub["mean_shap"]]
        ax.barh(sub["taxon"], sub["mean_abs_shap"], color=colours, alpha=0.85)
        ax.set_xlabel("Mean |SHAP value| (CLR units)", fontsize=8)
        label = COHORT_LABELS_S2.get(cohort, cohort)
        ax.set_title(label, fontsize=8.5, fontweight="bold",
                     color=COLORS.get(cohort, "black"))
        ax.tick_params(axis="y", labelsize=7)
        ax.tick_params(axis="x", labelsize=7)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    fig.suptitle("Top-15 SHAP taxa per cohort (LightGBM)\n"
                 "Red = AD-associated; Blue = CN-associated",
                 fontsize=9, y=1.01)
    plt.tight_layout()
    out = f"{FIGS}/supp_fig_s2.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"  Saved: {out}")


if __name__ == "__main__":
    print("Generating manuscript figures\n")

    over_df = pd.read_csv(f"{TABLES}/shap_taxa_overlap.csv")
    print("  overlap cols:", over_df.columns.tolist())
    print("  overlap head:\n", over_df.head(6))

    make_fig1()
    make_fig2()
    make_fig3()
    make_fig4()
    make_fig5()
    make_fig6()
    make_supp_fig_s1()
    make_supp_fig_s2()

    print("\nAll figures saved to results/figures/")
