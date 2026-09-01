from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
import lightgbm as lgb
import shap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

PROJECT_ROOT  = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR    = PROJECT_ROOT / "results" / "tables"
FIGURES_DIR   = PROJECT_ROOT / "results" / "figures"
MODEL_DIR     = PROJECT_ROOT / "results" / "model_outputs" / "shap"

for d in [TABLES_DIR, FIGURES_DIR, MODEL_DIR]:
    d.mkdir(parents=True, exist_ok=True)

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]
RANDOM_STATE = 42
TOP_N        = 20
TOP_DISPLAY  = 15
MIN_COHORTS  = 2

LOGREG_GRID = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}
LGBM_GRID   = {
    "n_estimators": [50, 100, 200],
    "learning_rate": [0.05, 0.1],
    "max_depth":     [3, 5],
    "num_leaves":    [15, 31],
}

COHORT_COLOURS = {
    "zhuang2018":    "#E41A1C",
    "ling2020":      "#FF7F00",
    "shanghai2022":  "#984EA3",
    "kazakhstan2022":"#377EB8",
}
COHORT_LABELS = {
    "zhuang2018":    "Zhuang 2018\n(China, n=86)",
    "ling2020":      "Ling 2020\n(China, n=171)",
    "shanghai2022":  "Zhu 2022\n(China, n=60)†",
    "kazakhstan2022":"Kazakhstan 2022\n(KZ, n=84)",
}


def load_data():
    clr_df  = pd.read_csv(PROCESSED_DIR / "clr_matrix.csv", index_col=0)
    diag_df = pd.read_csv(PROCESSED_DIR / "diagnosis_labels.csv")
    genus_cols = [c for c in clr_df.columns if c != "cohort"]

    clr_reset = clr_df.reset_index().rename(columns={"sample_id": "run_id"})
    merged = clr_reset.merge(
        diag_df[["run_id", "cohort", "diagnosis", "sample_name"]],
        on="run_id", suffixes=("_clr", "")
    ).drop(columns=["cohort_clr"], errors="ignore").set_index("run_id")

    merged_binary = merged[merged["diagnosis"].isin(["AD", "CN"])].copy()

    cohort_data = {}
    for cohort in SUPERVISED_COHORTS:
        sub = merged_binary[merged_binary["cohort"] == cohort]
        X = sub[genus_cols].to_numpy()
        y = (sub["diagnosis"] == "AD").astype(int).to_numpy()

        groups = None
        if cohort == "shanghai2022":
            groups = sub["sample_name"].map(
                lambda s: "_".join(s.split("_")[1:])
            ).to_numpy()

        cohort_data[cohort] = {
            "X": X, "y": y, "groups": groups,
            "ids": list(sub.index),
            "n": len(y),
        }
    return cohort_data, genus_cols


def fit_best_logreg(X, y, groups=None):
    inner_cv = (StratifiedGroupKFold(n_splits=5) if groups is not None
                else StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE))
    base = LogisticRegression(solver="saga", l1_ratio=0.5,
                              max_iter=10000, random_state=RANDOM_STATE)
    gs = GridSearchCV(base, LOGREG_GRID, cv=inner_cv, scoring="roc_auc", n_jobs=-1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if groups is not None:
            gs.fit(X, y, groups=groups)
        else:
            gs.fit(X, y)
    return gs.best_estimator_, gs.best_params_


def fit_best_lgbm(X, y, groups=None):
    inner_cv = (StratifiedGroupKFold(n_splits=5) if groups is not None
                else StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE))
    base = lgb.LGBMClassifier(objective="binary", n_jobs=4,
                               random_state=RANDOM_STATE, verbose=-1)
    gs = GridSearchCV(base, LGBM_GRID, cv=inner_cv, scoring="roc_auc", n_jobs=-1)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if groups is not None:
            gs.fit(X, y, groups=groups)
        else:
            gs.fit(X, y)
    return gs.best_estimator_, gs.best_params_


def shap_logreg(model, X_background, X_eval, genus_cols):
    # X_background must differ from X_eval: for a linear model under
    # interventional SHAP, shap_ij ~= beta_j*(x_ij - mean_j(background)), so
    # if background==eval, sum_i shap_ij == beta_j*(sum_i x_ij - n*mean_j(X))
    # == 0 EXACTLY (any set's deviations from its own mean sum to zero). That
    # forces mean_shap to ~1e-17 noise regardless of the true feature effect,
    # silently breaking every AD/CN direction call downstream. Using a
    # different reference population (training fold / other cohorts) avoids
    # this degeneracy.
    explainer  = shap.LinearExplainer(model, X_background,
                                      feature_perturbation="interventional")
    shap_vals  = explainer.shap_values(X_eval)
    return shap_vals


def shap_lgbm(model, X, genus_cols):
    explainer = shap.TreeExplainer(model)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        shap_vals = explainer.shap_values(X, check_additivity=False)
    # version-dependent return shape, annoying
    if isinstance(shap_vals, list):
        return shap_vals[1]
    return shap_vals


def compute_importance(shap_vals, genus_cols, cohort, model_name, y):
    mean_abs = np.abs(shap_vals).mean(axis=0)
    mean_shp = shap_vals.mean(axis=0)
    mean_ad  = shap_vals[y == 1].mean(axis=0) if (y == 1).any() else np.full(len(genus_cols), np.nan)
    mean_cn  = shap_vals[y == 0].mean(axis=0) if (y == 0).any() else np.full(len(genus_cols), np.nan)

    df = pd.DataFrame({
        "taxon":         genus_cols,
        "cohort":        cohort,
        "model":         model_name,
        "mean_abs_shap": mean_abs,
        "mean_shap":     mean_shp,
        "mean_shap_AD":  mean_ad,
        "mean_shap_CN":  mean_cn,
    }).sort_values("mean_abs_shap", ascending=False).reset_index(drop=True)
    df["rank"] = df.index + 1
    return df


OUTER_FOLDS = 10  # matches 03_within_cohort_baseline.py's nested-CV outer loop


def run_within_cohort_shap(cohort_data, genus_cols, outer_folds=OUTER_FOLDS):
    # Out-of-fold SHAP: for each outer fold, fit on the other folds (with the
    # same inner-CV hyperparameter search fit_best_* already does) and explain
    # only the held-out fold. In-sample SHAP (fit and explain on all samples)
    # can reflect overfitting, especially for small cohorts (n=60-171).
    all_rows = []

    for cohort in SUPERVISED_COHORTS:
        cd = cohort_data[cohort]
        X, y, groups = cd["X"], cd["y"], cd["groups"]
        n_samples, n_genera = X.shape
        print(f"\n  {cohort}  (n={cd['n']})  — {outer_folds}-fold OOF SHAP")

        if groups is not None:
            outer_cv = StratifiedGroupKFold(n_splits=outer_folds)
            split_iter = list(outer_cv.split(X, y, groups=groups))
        else:
            outer_cv = StratifiedKFold(n_splits=outer_folds, shuffle=True,
                                       random_state=RANDOM_STATE)
            split_iter = list(outer_cv.split(X, y))

        oof_shap_lr = np.zeros((n_samples, n_genera))
        oof_shap_gb = np.zeros((n_samples, n_genera))

        for fold_idx, (train_idx, test_idx) in enumerate(split_iter):
            X_tr, X_te = X[train_idx], X[test_idx]
            y_tr, y_te = y[train_idx], y[test_idx]
            g_tr = groups[train_idx] if groups is not None else None

            logreg, lp = fit_best_logreg(X_tr, y_tr, g_tr)
            oof_shap_lr[test_idx] = shap_logreg(logreg, X_tr, X_te, genus_cols)

            lgbm_m, lbp = fit_best_lgbm(X_tr, y_tr, g_tr)
            oof_shap_gb[test_idx] = shap_lgbm(lgbm_m, X_te, genus_cols)

            print(f"    fold {fold_idx+1:2d}/{outer_folds}: "
                  f"n_test={len(test_idx)}  logreg C={lp['C']}  "
                  f"lgbm n_est={lbp['n_estimators']}, lr={lbp['learning_rate']}")

        imp_lr = compute_importance(oof_shap_lr, genus_cols, cohort, "logreg", y)
        imp_gb = compute_importance(oof_shap_gb, genus_cols, cohort, "lgbm", y)
        all_rows.append(imp_lr)
        all_rows.append(imp_gb)

        print(f"    OOF top taxon (logreg): {imp_lr.iloc[0]['taxon']} "
              f"(|SHAP|={imp_lr.iloc[0]['mean_abs_shap']:.4f}, "
              f"dir={'AD↑' if imp_lr.iloc[0]['mean_shap']>0 else 'CN↑'})")
        print(f"    OOF top taxon (lgbm):   {imp_gb.iloc[0]['taxon']} "
              f"(|SHAP|={imp_gb.iloc[0]['mean_abs_shap']:.4f}, "
              f"dir={'AD↑' if imp_gb.iloc[0]['mean_shap']>0 else 'CN↑'})")

    return pd.concat(all_rows, ignore_index=True)


def run_loco_shap(cohort_data, genus_cols):
    all_rows = []

    for held_out in SUPERVISED_COHORTS:
        train_cohorts = [c for c in SUPERVISED_COHORTS if c != held_out]
        X_train = np.vstack([cohort_data[c]["X"] for c in train_cohorts])
        y_train = np.concatenate([cohort_data[c]["y"] for c in train_cohorts])
        X_test  = cohort_data[held_out]["X"]
        y_test  = cohort_data[held_out]["y"]

        print(f"\n  LOCO held-out={held_out}  "
              f"train={'+'.join(train_cohorts)}  (train n={len(y_train)}, test n={len(X_test)})")

        logreg, lp = fit_best_logreg(X_train, y_train)
        sv_lr = shap_logreg(logreg, X_train, X_test, genus_cols)
        imp_lr = compute_importance(sv_lr, genus_cols,
                                    f"loco_test_{held_out}", "logreg", y_test)
        imp_lr["train_cohorts"] = "+".join(train_cohorts)
        all_rows.append(imp_lr)
        print(f"    logreg best C={lp['C']}. "
              f"Top taxon: {imp_lr.iloc[0]['taxon']}")

        lgbm_m, lbp = fit_best_lgbm(X_train, y_train)
        sv_gb = shap_lgbm(lgbm_m, X_test, genus_cols)
        imp_gb = compute_importance(sv_gb, genus_cols,
                                    f"loco_test_{held_out}", "lgbm", y_test)
        imp_gb["train_cohorts"] = "+".join(train_cohorts)
        all_rows.append(imp_gb)
        print(f"    lgbm n_est={lbp['n_estimators']}, lr={lbp['learning_rate']}. "
              f"Top taxon: {imp_gb.iloc[0]['taxon']}")

    return pd.concat(all_rows, ignore_index=True)


def jaccard(set_a, set_b):
    overlap = len(set_a & set_b)
    union = len(set_a | set_b)
    return overlap / union if union else 0.0


def compute_overlap(importance_df, model_name, top_n=TOP_N):
    top_sets = {}
    for cohort in SUPERVISED_COHORTS:
        sub = importance_df[
            (importance_df["cohort"] == cohort) &
            (importance_df["model"] == model_name)
        ].head(top_n)
        top_sets[cohort] = set(sub["taxon"])

    mat = pd.DataFrame(index=SUPERVISED_COHORTS, columns=SUPERVISED_COHORTS,
                       dtype=float)
    for c1 in SUPERVISED_COHORTS:
        for c2 in SUPERVISED_COHORTS:
            mat.loc[c1, c2] = jaccard(top_sets[c1], top_sets[c2])
    return mat, top_sets


N_JACCARD_PERM = 10_000


def jaccard_null_distribution(n_genera, top_n=TOP_N, n_perm=N_JACCARD_PERM,
                              rng=None):
    # Null: two independent random top-n draws from the same n_genera pool.
    # Analytical E[J] for two random size-k sets from n items (no overlap
    # constraint other than set size): E[|A∩B|] = k^2/n, so
    # E[J] ≈ E[|A∩B|] / (2k - E[|A∩B|]) = k / (2n - k) for k << n.
    if rng is None:
        rng = np.random.default_rng(RANDOM_STATE)
    null_j = np.empty(n_perm)
    for i in range(n_perm):
        set_a = set(rng.choice(n_genera, size=top_n, replace=False))
        set_b = set(rng.choice(n_genera, size=top_n, replace=False))
        null_j[i] = jaccard(set_a, set_b)
    return null_j


def compute_jaccard_null_baseline(within_imp, n_genera, top_n=TOP_N,
                                  n_perm=N_JACCARD_PERM):
    rng = np.random.default_rng(RANDOM_STATE)
    rows = []
    for model_name in ["logreg", "lgbm"]:
        mat, _ = compute_overlap(within_imp, model_name, top_n=top_n)
        off_diag = [mat.loc[c1, c2]
                    for c1 in SUPERVISED_COHORTS
                    for c2 in SUPERVISED_COHORTS if c1 != c2]
        observed_mean = float(np.mean(off_diag))

        null_j = jaccard_null_distribution(n_genera, top_n=top_n,
                                           n_perm=n_perm, rng=rng)
        mean_null = float(null_j.mean())
        ci_lo = float(np.percentile(null_j, 2.5))
        ci_hi = float(np.percentile(null_j, 97.5))
        p_value = float((null_j >= observed_mean).mean())
        analytical_e = top_n / (2 * n_genera - top_n)

        rows.append({
            "model":                model_name,
            "mean_null_jaccard":    round(mean_null, 4),
            "ci_lower_null":        round(ci_lo, 4),
            "ci_upper_null":        round(ci_hi, 4),
            "analytical_expected":  round(analytical_e, 4),
            "observed_mean_jaccard": round(observed_mean, 4),
            "p_value":              round(p_value, 5),
        })
        where = ("ABOVE" if observed_mean > ci_hi else
                 "BELOW" if observed_mean < ci_lo else "WITHIN")
        print(f"  {model_name}: observed mean Jaccard={observed_mean:.4f} | "
              f"null mean={mean_null:.4f} [{ci_lo:.4f}, {ci_hi:.4f}] "
              f"(analytical E[J]={analytical_e:.4f}) | "
              f"observed is {where} the null 95% CI | p={p_value:.5f}")

    return pd.DataFrame(rows)


def directional_analysis(importance_df, model_name, top_n=TOP_N,
                         min_cohorts=MIN_COHORTS):
    top_sets = {}
    for cohort in SUPERVISED_COHORTS:
        sub = importance_df[
            (importance_df["cohort"] == cohort) &
            (importance_df["model"] == model_name)
        ].head(top_n)
        top_sets[cohort] = set(sub["taxon"])

    all_taxa = {}
    for cohort, s in top_sets.items():
        for t in s:
            all_taxa[t] = all_taxa.get(t, 0) + 1
    shared_taxa = {t for t, cnt in all_taxa.items() if cnt >= min_cohorts}

    rows = []
    for taxon in sorted(shared_taxa):
        for cohort in SUPERVISED_COHORTS:
            sub = importance_df[
                (importance_df["cohort"] == cohort) &
                (importance_df["model"] == model_name) &
                (importance_df["taxon"] == taxon)
            ]
            if len(sub) == 0:
                continue
            rows.append({
                "taxon":              taxon,
                "cohort":             cohort,
                "mean_shap":          sub["mean_shap"].values[0],
                "mean_abs_shap":      sub["mean_abs_shap"].values[0],
                "rank_in_cohort":     sub["rank"].values[0],
                "in_top_n":           taxon in top_sets[cohort],
                "n_cohorts_in_top":   all_taxa[taxon],
            })

    return pd.DataFrame(rows)


def plot_top15_per_cohort(importance_df, model_name, out_path):
    model_label = "Logistic Regression" if model_name == "logreg" else "LightGBM"
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for ax, cohort in zip(axes, SUPERVISED_COHORTS):
        sub = importance_df[
            (importance_df["cohort"] == cohort) &
            (importance_df["model"] == model_name)
        ].head(TOP_DISPLAY).iloc[::-1]

        colours = ["#D62728" if v > 0 else "#1F77B4"
                   for v in sub["mean_shap"]]

        ax.barh(sub["taxon"], sub["mean_abs_shap"], color=colours, alpha=0.85)
        ax.set_xlabel("Mean |SHAP value| (CLR units)", fontsize=9)
        ax.set_title(f"{cohort}\n({COHORT_LABELS[cohort].replace(chr(10), ' ')})",
                     fontsize=9, fontweight="bold",
                     color=COHORT_COLOURS[cohort])
        ax.tick_params(axis="y", labelsize=8)
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", alpha=0.3)

        p_ad = mpatches.Patch(color="#D62728", label="AD↑ (positive SHAP)")
        p_cn = mpatches.Patch(color="#1F77B4", label="CN↑ (negative SHAP)")
        ax.legend(handles=[p_ad, p_cn], fontsize=7, loc="lower right")

    fig.suptitle(
        f"Top {TOP_DISPLAY} Taxa by SHAP Importance — {model_label}\n"
        f"(Red = higher CLR abundance → AD prediction, "
        f"Blue = higher CLR abundance → CN prediction)",
        fontsize=10, y=1.01
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_jaccard_heatmap(importance_df, out_path):
    models      = ["logreg", "lgbm"]
    m_labels    = {"logreg": "Logistic Regression", "lgbm": "LightGBM"}
    short_names = {c: c.replace("2022","'22").replace("2020","'20")
                     .replace("2018","'18").replace("kazakhstan","KZ")
                   for c in SUPERVISED_COHORTS}

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    for ax, model_name in zip(axes, models):
        mat, _ = compute_overlap(importance_df, model_name, top_n=TOP_N)
        labels  = [short_names[c] for c in SUPERVISED_COHORTS]

        sns.heatmap(mat.astype(float), ax=ax,
                    annot=True, fmt=".2f", cmap="YlOrRd",
                    vmin=0, vmax=1, linewidths=0.5, linecolor="white",
                    xticklabels=labels, yticklabels=labels,
                    cbar_kws={"label": "Jaccard", "shrink": 0.8})
        ax.set_title(f"{m_labels[model_name]}\n"
                     f"Jaccard overlap of top-{TOP_N} taxa",
                     fontsize=10)
        ax.tick_params(labelsize=9, rotation=0)

    fig.suptitle(
        f"Cross-Cohort Taxon Overlap — Top-{TOP_N} SHAP Taxa per Cohort\n"
        "(Jaccard = |intersection| / |union|; 1.0 = identical lists, 0.0 = no overlap)",
        fontsize=9
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_direction_dotplot(dir_df, model_name, out_path):
    if dir_df.empty:
        print(f"  No shared taxa at top-{TOP_N} across ≥{MIN_COHORTS} cohorts, skipping plot.")
        return

    model_label = "Logistic Regression" if model_name == "logreg" else "LightGBM"

    taxon_order = (
        dir_df.groupby("taxon")["mean_abs_shap"].mean()
        .sort_values(ascending=True)
        .index.tolist()
    )

    flip_taxa = set()
    for taxon in dir_df["taxon"].unique():
        signs = dir_df[dir_df["taxon"] == taxon]["mean_shap"].apply(np.sign)
        if signs.min() < 0 and signs.max() > 0:
            flip_taxa.add(taxon)

    fig_height = max(5, 0.45 * len(taxon_order))
    fig, ax = plt.subplots(figsize=(9, fig_height))

    for cohort in SUPERVISED_COHORTS:
        sub = dir_df[dir_df["cohort"] == cohort]
        y_pos = [taxon_order.index(t) for t in sub["taxon"] if t in taxon_order]
        x_val = [sub[sub["taxon"] == t]["mean_shap"].values[0]
                 for t in taxon_order if t in sub["taxon"].values]
        sizes = [sub[sub["taxon"] == t]["mean_abs_shap"].values[0] * 600
                 for t in taxon_order if t in sub["taxon"].values]

        ax.scatter(x_val, y_pos, label=cohort,
                   color=COHORT_COLOURS[cohort], s=sizes,
                   alpha=0.75, edgecolors="white", linewidths=0.5, zorder=3)

    ax.axvline(0, color="black", linewidth=0.8, linestyle="--", zorder=2)
    ax.set_yticks(range(len(taxon_order)))
    ax.set_yticklabels(
        [f"{'★ ' if t in flip_taxa else ''}{t}" for t in taxon_order],
        fontsize=8
    )
    for i, t in enumerate(taxon_order):
        if t in flip_taxa:
            ax.axhspan(i - 0.5, i + 0.5, color="#FFFACD", alpha=0.6, zorder=1)

    ax.set_xlabel("Mean SHAP value (positive → AD, negative → CN)", fontsize=10)
    ax.set_title(
        f"Cross-Cohort SHAP Direction — {model_label}\n"
        f"Taxa in top-{TOP_N} for ≥{MIN_COHORTS} cohorts "
        f"(★ = directional flip; dot size ∝ |SHAP|)",
        fontsize=10
    )
    ax.legend(title="Cohort", fontsize=8, title_fontsize=8,
              loc="lower right", framealpha=0.9)
    ax.grid(axis="x", alpha=0.25)
    ax.set_xlim(left=None)

    ax.text(0.01, 0.01, "← CN-associated", transform=ax.transAxes,
            fontsize=8, color="#1F77B4", va="bottom")
    ax.text(0.99, 0.01, "AD-associated →", transform=ax.transAxes,
            fontsize=8, color="#D62728", va="bottom", ha="right")

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}  "
          f"({len(taxon_order)} taxa shown, {len(flip_taxa)} direction flips ★)")


def plot_loco_direction(loco_df, model_name, out_path):
    if loco_df.empty:
        return

    model_label = "Logistic Regression" if model_name == "logreg" else "LightGBM"

    sub_model = loco_df[loco_df["model"] == model_name].copy()
    test_cohorts = [c for c in sub_model["cohort"].unique()
                    if c.startswith("loco_test_")]

    top_sets = {}
    for tc in test_cohorts:
        held = tc.replace("loco_test_", "")
        top_sets[held] = set(
            sub_model[sub_model["cohort"] == tc].head(15)["taxon"]
        )
    taxa_counts = {}
    for s in top_sets.values():
        for t in s:
            taxa_counts[t] = taxa_counts.get(t, 0) + 1
    shared = sorted({t for t, n in taxa_counts.items() if n >= 2},
                    key=lambda t: -taxa_counts[t])[:25]

    if not shared:
        print("  No shared LOCO taxa to plot, skipping.")
        return

    taxon_order = sorted(shared, key=lambda t:
        sub_model[sub_model["taxon"] == t]["mean_abs_shap"].mean())

    fig_h = max(5, 0.4 * len(taxon_order))
    fig, ax = plt.subplots(figsize=(9, fig_h))

    for tc in test_cohorts:
        held = tc.replace("loco_test_", "")
        sub = sub_model[sub_model["cohort"] == tc]
        y_pos, x_val, sizes = [], [], []
        for t in taxon_order:
            row = sub[sub["taxon"] == t]
            if len(row):
                y_pos.append(taxon_order.index(t))
                x_val.append(row["mean_shap"].values[0])
                sizes.append(row["mean_abs_shap"].values[0] * 600)
        ax.scatter(x_val, y_pos, label=f"LOCO→{held}",
                   color=COHORT_COLOURS.get(held, "grey"),
                   s=sizes, alpha=0.75, edgecolors="white",
                   linewidths=0.5, zorder=3)

    ax.axvline(0, color="black", lw=0.8, ls="--")
    ax.set_yticks(range(len(taxon_order)))
    ax.set_yticklabels(taxon_order, fontsize=8)
    ax.set_xlabel("Mean SHAP value on test cohort (positive → AD, negative → CN)",
                  fontsize=9)
    ax.set_title(
        f"LOCO SHAP Direction — {model_label}\n"
        f"SHAP of LOCO model evaluated on held-out test cohort",
        fontsize=10
    )
    ax.legend(title="LOCO test cohort", fontsize=8, title_fontsize=8,
              loc="lower right", framealpha=0.9)
    ax.grid(axis="x", alpha=0.25)
    ax.text(0.01, 0.01, "← CN-associated", transform=ax.transAxes,
            fontsize=8, color="#1F77B4", va="bottom")
    ax.text(0.99, 0.01, "AD-associated →", transform=ax.transAxes,
            fontsize=8, color="#D62728", va="bottom", ha="right")
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def main():
    print("08_shap_taxa_comparison.py — Phase 6")

    cohort_data, genus_cols = load_data()
    print(f"\nLoaded {len(genus_cols)} genera, {len(SUPERVISED_COHORTS)} cohorts")

    print("\nWithin-cohort SHAP")
    within_imp = run_within_cohort_shap(cohort_data, genus_cols)
    within_imp.to_csv(TABLES_DIR / "shap_within_cohort_importance.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/shap_within_cohort_importance.csv")

    print("\nLOCO SHAP (train on N-1 cohorts, SHAP on held-out)")
    loco_imp = run_loco_shap(cohort_data, genus_cols)
    loco_imp.to_csv(TABLES_DIR / "shap_loco_importance.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/shap_loco_importance.csv")

    print(f"\nOverlap analysis (top-{TOP_N} per cohort)")

    overlap_rows = []
    for model_name in ["logreg", "lgbm"]:
        for top_n in [10, 20, 50]:
            mat, top_sets = compute_overlap(within_imp, model_name, top_n=top_n)
            off_diag = [mat.loc[c1, c2]
                        for c1 in SUPERVISED_COHORTS
                        for c2 in SUPERVISED_COHORTS if c1 != c2]
            mean_j = np.mean(off_diag)
            print(f"  {model_name} top-{top_n:2d}: "
                  f"mean pairwise Jaccard = {mean_j:.3f}")
            for c1 in SUPERVISED_COHORTS:
                for c2 in SUPERVISED_COHORTS:
                    overlap_rows.append({
                        "model": model_name,
                        "top_n": top_n,
                        "cohort_a": c1,
                        "cohort_b": c2,
                        "jaccard": round(float(mat.loc[c1, c2]), 4),
                    })

    overlap_df = pd.DataFrame(overlap_rows)
    overlap_df.to_csv(TABLES_DIR / "shap_taxa_overlap.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/shap_taxa_overlap.csv")

    print(f"\nJaccard null permutation baseline (top-{TOP_N}, {N_JACCARD_PERM:,} draws)")
    null_df = compute_jaccard_null_baseline(within_imp, len(genus_cols), top_n=TOP_N,
                                            n_perm=N_JACCARD_PERM)
    null_df.to_csv(TABLES_DIR / "jaccard_null.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/jaccard_null.csv")

    print(f"\nDirectional SHAP analysis (shared top-{TOP_N} taxa)")

    dir_rows = []
    for model_name in ["logreg", "lgbm"]:
        dir_df = directional_analysis(within_imp, model_name)
        dir_df["model"] = model_name
        dir_rows.append(dir_df)

        flip_taxa = set()
        for taxon in dir_df["taxon"].unique():
            signs = dir_df[dir_df["taxon"] == taxon]["mean_shap"].apply(np.sign)
            if signs.min() < 0 and signs.max() > 0:
                flip_taxa.add(taxon)

        shared_count = dir_df["taxon"].nunique()
        print(f"\n  {model_name}: {shared_count} taxa shared across "
              f"≥{MIN_COHORTS} cohorts' top-{TOP_N}")
        print(f"  {model_name}: {len(flip_taxa)} taxa with directional flip "
              f"(AD in one cohort, CN in another):")
        for t in sorted(flip_taxa):
            cohort_signs = dir_df[dir_df["taxon"] == t][["cohort","mean_shap"]].copy()
            direction = {row["cohort"]: "AD↑" if row["mean_shap"] > 0 else "CN↑"
                         for _, row in cohort_signs.iterrows()}
            print(f"    {t}: {direction}")

    dir_all = pd.concat(dir_rows, ignore_index=True)
    dir_all.to_csv(TABLES_DIR / "shap_directional_flips.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/shap_directional_flips.csv")

    print("\nGenerating figures")

    for model_name in ["logreg", "lgbm"]:
        plot_top15_per_cohort(
            within_imp, model_name,
            FIGURES_DIR / f"shap_top15_per_cohort_{model_name}.png"
        )
    plot_jaccard_heatmap(
        within_imp,
        FIGURES_DIR / "shap_overlap_jaccard.png"
    )

    for model_name in ["logreg", "lgbm"]:
        dir_sub = dir_all[dir_all["model"] == model_name]
        plot_direction_dotplot(
            dir_sub, model_name,
            FIGURES_DIR / f"shap_direction_dotplot_{model_name}.png"
        )
        plot_loco_direction(
            loco_imp, model_name,
            FIGURES_DIR / f"shap_loco_direction_{model_name}.png"
        )

    print("\nKey findings")

    for model_name in ["logreg", "lgbm"]:
        mat, top_sets = compute_overlap(within_imp, model_name, top_n=TOP_N)
        off = [mat.loc[c1, c2]
               for c1 in SUPERVISED_COHORTS
               for c2 in SUPERVISED_COHORTS if c1 != c2]
        dir_sub = dir_all[dir_all["model"] == model_name]
        n_shared = dir_sub["taxon"].nunique()
        flip_set = {
            t for t in dir_sub["taxon"].unique()
            if (dir_sub[dir_sub["taxon"]==t]["mean_shap"].apply(np.sign).min() < 0 and
                dir_sub[dir_sub["taxon"]==t]["mean_shap"].apply(np.sign).max() > 0)
        }

        print(f"\n  {model_name}")
        print(f"    Mean pairwise Jaccard (top-{TOP_N}): {np.mean(off):.3f}")
        print(f"    Taxa shared ≥{MIN_COHORTS} cohorts: {n_shared}")
        print(f"    Direction flips: {len(flip_set)}")
        if flip_set:
            print(f"    Flipped taxa: {sorted(flip_set)}")

        print(f"\n    Top-3 taxa per cohort ({model_name}):")
        for cohort in SUPERVISED_COHORTS:
            top3 = within_imp[
                (within_imp["cohort"] == cohort) &
                (within_imp["model"] == model_name)
            ].head(3)[["taxon", "mean_abs_shap", "mean_shap"]].copy()
            top3["dir"] = top3["mean_shap"].apply(lambda x: "AD↑" if x > 0 else "CN↑")
            lines = [f"{r['taxon']} ({r['dir']}, |SHAP|={r['mean_abs_shap']:.3f})"
                     for _, r in top3.iterrows()]
            print(f"      {cohort}: {' | '.join(lines)}")

    print("\nPhase 6 complete.")
    print("Outputs:")
    for f in [
        "shap_within_cohort_importance.csv", "shap_loco_importance.csv",
        "shap_taxa_overlap.csv", "shap_directional_flips.csv", "jaccard_null.csv"
    ]:
        print(f"  {TABLES_DIR}/{f}")
    for f in [
        "shap_top15_per_cohort_logreg.png", "shap_top15_per_cohort_lgbm.png",
        "shap_overlap_jaccard.png",
        "shap_direction_dotplot_logreg.png", "shap_direction_dotplot_lgbm.png",
        "shap_loco_direction_logreg.png", "shap_loco_direction_lgbm.png",
    ]:
        print(f"  {FIGURES_DIR}/{f}")
    print("\nNext step: manuscript draft (manuscript/draft.md)")


if __name__ == "__main__":
    main()
