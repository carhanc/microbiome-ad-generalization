from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold, StratifiedGroupKFold
from sklearn.metrics import roc_auc_score
import lightgbm as lgb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
TABLES_DIR    = PROJECT_ROOT / "results" / "tables"
FIGURES_DIR   = PROJECT_ROOT / "results" / "figures"
MODEL_DIR     = PROJECT_ROOT / "results" / "model_outputs" / "corrected"

for d in [TABLES_DIR, FIGURES_DIR, MODEL_DIR]:
    d.mkdir(parents=True, exist_ok=True)

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]
RANDOM_STATE  = 42
INNER_FOLDS   = 5
N_BOOTSTRAP   = 1000

LOGREG_GRID = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}
LGBM_GRID   = {
    "n_estimators": [50, 100, 200],
    "learning_rate": [0.05, 0.1],
    "max_depth":     [3, 5],
    "num_leaves":    [15, 31],
}

CORRECTION_METHODS = {
    "combatseq": PROCESSED_DIR / "clr_matrix_combatseq.csv",
    "mmuphin":   PROCESSED_DIR / "clr_matrix_mmuphin.csv",
}


def bootstrap_auc_ci(y_true, y_score, n_boot=N_BOOTSTRAP, alpha=0.05):
    rng = np.random.default_rng(RANDOM_STATE)
    n = len(y_true)
    boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        boot.append(roc_auc_score(yt, ys))
    boot = np.array(boot)
    return float(np.percentile(boot, 100*alpha/2)), float(np.percentile(boot, 100*(1-alpha/2)))


def tune_and_fit(X_train, y_train, model_name):
    inner_cv = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    if model_name == "logreg":
        base  = LogisticRegression(solver="saga", l1_ratio=0.5,
                                   max_iter=10000, random_state=RANDOM_STATE)
        grid  = LOGREG_GRID
    else:
        base  = lgb.LGBMClassifier(objective="binary", n_jobs=4,
                                    random_state=RANDOM_STATE, verbose=-1)
        grid  = LGBM_GRID
    gs = GridSearchCV(base, grid, cv=inner_cv, scoring="roc_auc",
                      n_jobs=-1, refit=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gs.fit(X_train, y_train)
    return gs.best_estimator_, gs.best_params_


def eval_test(model, X_test, y_test):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y_score = model.predict_proba(X_test)[:, 1]
    auc = roc_auc_score(y_test, y_score)
    ci_lo, ci_hi = bootstrap_auc_ci(y_test, y_score)
    return {"auc": round(auc, 4), "ci_lo": round(ci_lo, 4),
            "ci_hi": round(ci_hi, 4), "y_score": y_score}


def load_corrected_data(clr_path):
    df = pd.read_csv(clr_path)
    id_col = df.columns[0]
    genus_cols = [c for c in df.columns if c not in [id_col, "cohort", "diagnosis"]]

    by_cohort = {}
    for cohort in SUPERVISED_COHORTS:
        sub = df[(df["cohort"] == cohort) & df["diagnosis"].isin(["AD", "CN"])].copy()
        labels = (sub["diagnosis"] == "AD").astype(int).to_numpy()
        by_cohort[cohort] = {
            "X": sub[genus_cols].to_numpy(),
            "y": labels,
            "ids":  list(sub[id_col]),
            "n":    len(labels),
            "n_AD": int(labels.sum()),
            "n_CN": int((labels == 0).sum()),
        }
    return by_cohort, genus_cols


def run_loco(cohort_data, method_label):
    results = []
    for held_out in SUPERVISED_COHORTS:
        train_cohorts = [c for c in SUPERVISED_COHORTS if c != held_out]
        X_train = np.vstack([cohort_data[c]["X"] for c in train_cohorts])
        y_train = np.concatenate([cohort_data[c]["y"] for c in train_cohorts])
        X_test  = cohort_data[held_out]["X"]
        y_test  = cohort_data[held_out]["y"]

        print(f"  {method_label} LOCO held-out={held_out} "
              f"| train n={len(y_train)}, test n={len(y_test)}")

        for model_name in ["logreg", "lgbm"]:
            model, best_p = tune_and_fit(X_train, y_train, model_name)
            m = eval_test(model, X_test, y_test)
            print(f"    {model_name} AUC={m['auc']:.4f} "
                  f"CI=[{m['ci_lo']:.4f},{m['ci_hi']:.4f}]")
            results.append({
                "correction":    method_label,
                "experiment":    "LOCO",
                "train_cohorts": "+".join(train_cohorts),
                "test_cohort":   held_out,
                "model":         model_name,
                "n_train":       len(y_train),
                "n_test":        len(y_test),
                "auc":           m["auc"],
                "ci_lo":         m["ci_lo"],
                "ci_hi":         m["ci_hi"],
            })
    return results


def run_pairwise(cohort_data, method_label):
    results = []
    for train_c in SUPERVISED_COHORTS:
        for test_c in SUPERVISED_COHORTS:
            if train_c == test_c:
                continue
            X_train = cohort_data[train_c]["X"]
            y_train = cohort_data[train_c]["y"]
            X_test  = cohort_data[test_c]["X"]
            y_test  = cohort_data[test_c]["y"]

            for model_name in ["logreg", "lgbm"]:
                model, _ = tune_and_fit(X_train, y_train, model_name)
                m = eval_test(model, X_test, y_test)
                results.append({
                    "correction":   method_label,
                    "experiment":   "pairwise",
                    "train_cohort": train_c,
                    "test_cohort":  test_c,
                    "model":        model_name,
                    "n_train":      len(y_train),
                    "n_test":       len(y_test),
                    "auc":          m["auc"],
                    "ci_lo":        m["ci_lo"],
                    "ci_hi":        m["ci_hi"],
                })
    return results


def run_within_cohort_check(cohort_data, method_label, diag_df):
    OUTER_FOLDS = 10
    results = []

    for cohort in SUPERVISED_COHORTS:
        cd   = cohort_data[cohort]
        X, y = cd["X"], cd["y"]

        groups = None
        if cohort == "shanghai2022":
            run_id_to_sname = diag_df.set_index("run_id")["sample_name"]
            snames = pd.Series(cd["ids"]).map(run_id_to_sname)
            groups = snames.map(lambda s: "_".join(s.split("_")[1:])).to_numpy()

        outer_cv = (StratifiedGroupKFold(n_splits=OUTER_FOLDS) if groups is not None
                    else StratifiedKFold(n_splits=OUTER_FOLDS, shuffle=True,
                                         random_state=RANDOM_STATE))
        inner_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)

        for model_name in ["logreg", "lgbm"]:
            oof = np.zeros(len(y))
            split_iter = (outer_cv.split(X, y, groups=groups) if groups is not None
                          else outer_cv.split(X, y))
            for train_i, test_i in split_iter:
                Xtr, Xte = X[train_i], X[test_i]
                ytr, yte = y[train_i], y[test_i]
                model, _ = tune_and_fit(Xtr, ytr, model_name)
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    oof[test_i] = model.predict_proba(Xte)[:, 1]

            auc = roc_auc_score(y, oof)
            ci_lo, ci_hi = bootstrap_auc_ci(y, oof)
            print(f"  {method_label} within {cohort}/{model_name}: "
                  f"AUC={auc:.4f} CI=[{ci_lo:.4f},{ci_hi:.4f}]")
            results.append({
                "correction": method_label,
                "cohort":     cohort,
                "model":      model_name,
                "auc_within_corrected": round(auc, 4),
                "ci_lo":      round(ci_lo, 4),
                "ci_hi":      round(ci_hi, 4),
            })
    return results


def plot_auc_comparison(comparison_df, out_path):
    models      = ["logreg", "lgbm"]
    m_labels    = {"logreg": "Logistic Regression", "lgbm": "LightGBM"}
    corrections = ["uncorrected", "combatseq", "mmuphin"]
    corr_labels = {"uncorrected": "Uncorrected",
                   "combatseq":   "ComBat-seq",
                   "mmuphin":     "MMUPHin"}
    corr_colors = {"uncorrected": "#888888",
                   "combatseq":   "#2196F3",
                   "mmuphin":     "#4CAF50"}
    cohort_labels = {
        "zhuang2018":    "Zhuang 2018\n(China, n=86)",
        "ling2020":      "Ling 2020\n(China, n=171)",
        "shanghai2022":  "Shanghai 2022\n(China, n=120)",
        "kazakhstan2022":"Kazakhstan 2022\n(KZ, n=84)",
    }

    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5), sharey=True)

    for ax, model_name in zip(axes, models):
        sub = comparison_df[comparison_df["model"] == model_name]
        x        = np.arange(len(SUPERVISED_COHORTS))
        n_bars   = len(corrections)
        width    = 0.22
        offsets  = np.linspace(-(n_bars-1)/2, (n_bars-1)/2, n_bars) * width

        for i, corr in enumerate(corrections):
            row = sub[sub["correction"] == corr]
            aucs     = [row[row["test_cohort"]==c]["loco_auc"].values[0]
                        if len(row[row["test_cohort"]==c]) else np.nan
                        for c in SUPERVISED_COHORTS]
            ci_los   = [row[row["test_cohort"]==c]["loco_auc"].values[0] -
                        row[row["test_cohort"]==c]["loco_ci_lo"].values[0]
                        if len(row[row["test_cohort"]==c]) else 0
                        for c in SUPERVISED_COHORTS]
            ci_his   = [row[row["test_cohort"]==c]["loco_ci_hi"].values[0] -
                        row[row["test_cohort"]==c]["loco_auc"].values[0]
                        if len(row[row["test_cohort"]==c]) else 0
                        for c in SUPERVISED_COHORTS]
            ax.bar(x + offsets[i], aucs, width,
                   label=corr_labels[corr],
                   color=corr_colors[corr],
                   alpha=0.85,
                   yerr=[ci_los, ci_his],
                   error_kw={"capsize": 3, "linewidth": 1, "ecolor": "black"})

        wc_auc = sub[sub["correction"] == "uncorrected"]
        for j, cohort in enumerate(SUPERVISED_COHORTS):
            wc_row = wc_auc[wc_auc["test_cohort"] == cohort]
            if len(wc_row):
                wc_val = wc_row["within_auc"].values[0]
                ax.plot(x[j], wc_val, marker="*", color="black",
                        markersize=10, zorder=5,
                        label="Within-cohort AUC" if j == 0 else "")

        ax.axhline(0.5, color="grey", lw=0.8, ls="--", label="Chance (0.5)")
        ax.set_xticks(x)
        ax.set_xticklabels([cohort_labels[c] for c in SUPERVISED_COHORTS],
                           fontsize=8, rotation=15, ha="right")
        ax.set_ylim(0.35, 1.05)
        ax.set_ylabel("LOCO AUC-ROC", fontsize=10)
        ax.set_title(m_labels[model_name], fontsize=11, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right")
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle(
        "LOCO Cross-Cohort AUC: Uncorrected vs Batch-Corrected\n"
        "(★ = within-cohort AUC; error bars = 95% bootstrap CI)",
        fontsize=11, y=1.02
    )
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def plot_pairwise_heatmap_corrected(pairwise_df, correction, model_name, out_path):
    df   = pairwise_df[(pairwise_df["correction"] == correction) &
                       (pairwise_df["model"] == model_name)]
    mat  = pd.DataFrame(np.nan, index=SUPERVISED_COHORTS, columns=SUPERVISED_COHORTS)
    for _, row in df.iterrows():
        mat.loc[row["train_cohort"], row["test_cohort"]] = row["auc"]

    label_map = {
        "zhuang2018":    "Zhuang\n2018\n(CN)",
        "ling2020":      "Ling\n2020\n(CN)",
        "shanghai2022":  "Shanghai\n2022\n(CN)",
        "kazakhstan2022":"Kazakh\n2022\n(KZ)",
    }
    short = [label_map[c] for c in SUPERVISED_COHORTS]

    fig, ax = plt.subplots(figsize=(6, 5))
    sns.heatmap(mat.astype(float), ax=ax, annot=True, fmt=".3f",
                cmap="RdYlGn", vmin=0.3, vmax=1.0, center=0.5,
                linewidths=0.5, linecolor="white",
                xticklabels=short, yticklabels=short,
                cbar_kws={"label": "AUC-ROC", "shrink": 0.8})
    corr_label = {"combatseq": "ComBat-seq", "mmuphin": "MMUPHin"}
    model_label = "LogReg" if model_name == "logreg" else "LightGBM"
    ax.set_title(f"Pairwise AUC — {corr_label[correction]} / {model_label}\n"
                 f"Row = train, Col = test", fontsize=9)
    ax.set_xlabel("Test cohort", fontsize=9)
    ax.set_ylabel("Train cohort", fontsize=9)
    ax.tick_params(labelsize=8)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {out_path}")


def main():
    print("07_corrected_generalization.py — Phase 5")

    diag_df = pd.read_csv(PROCESSED_DIR / "diagnosis_labels.csv")

    uncorrected_loco = pd.read_csv(TABLES_DIR / "loco_auc.csv")
    uncorrected_pw   = pd.read_csv(TABLES_DIR / "pairwise_auc.csv")
    within_ref       = pd.read_csv(TABLES_DIR / "within_cohort_auc.csv")
    within_ref       = within_ref[~within_ref["cohort"].str.endswith("3class")]

    all_loco     = []
    all_pairwise = []
    all_within   = []

    # correction was fit on all cohorts at once — slight leakage in LOCO, noted in paper
    for method, clr_path in CORRECTION_METHODS.items():
        if not clr_path.exists():
            print(f"\nSkipping {method}: {clr_path} not found. "
                  f"Run 06_batch_correction.R first.")
            continue

        print(f"\nMethod: {method}")

        cohort_data, genus_cols = load_corrected_data(clr_path)
        print(f"  Loaded {sum(cd['n'] for cd in cohort_data.values())} samples, "
              f"{len(genus_cols)} genera")
        for c in SUPERVISED_COHORTS:
            cd = cohort_data[c]
            print(f"  {c:20s}: n={cd['n']} ({cd['n_AD']} AD / {cd['n_CN']} CN)")

        print("\n  LOCO")
        all_loco += run_loco(cohort_data, method)

        print("\n  Pairwise")
        all_pairwise += run_pairwise(cohort_data, method)

        print("\n  Within-cohort check on corrected data")
        all_within += run_within_cohort_check(cohort_data, method, diag_df)

    loco_df     = pd.DataFrame(all_loco)
    pairwise_df = pd.DataFrame(all_pairwise)
    within_df   = pd.DataFrame(all_within)

    loco_df.to_csv(TABLES_DIR / "loco_auc_corrected.csv", index=False)
    pairwise_df.to_csv(TABLES_DIR / "pairwise_auc_corrected.csv", index=False)
    within_df.to_csv(TABLES_DIR / "within_cohort_auc_corrected.csv", index=False)

    uncorr_loco_tagged = uncorrected_loco.copy()
    uncorr_loco_tagged["correction"] = "uncorrected"

    within_ref_lookup = within_ref.set_index(["cohort", "model"])["auc_oof"]

    rows = []
    for _, unc_row in uncorr_loco_tagged.iterrows():
        cohort = unc_row["test_cohort"]
        model  = unc_row["model"]
        wc_auc = within_ref_lookup.get((cohort, model), np.nan)
        base_row = {
            "test_cohort":   cohort,
            "model":         model,
            "within_auc":    round(wc_auc, 4),
            "correction":    "uncorrected",
            "loco_auc":      unc_row["auc"],
            "loco_ci_lo":    unc_row["ci_lo"],
            "loco_ci_hi":    unc_row["ci_hi"],
            "loco_drop":     round(wc_auc - unc_row["auc"], 4),
        }
        rows.append(base_row)

    for _, corr_row in loco_df.iterrows():
        cohort = corr_row["test_cohort"]
        model  = corr_row["model"]
        wc_auc = within_ref_lookup.get((cohort, model), np.nan)
        rows.append({
            "test_cohort":   cohort,
            "model":         model,
            "within_auc":    round(wc_auc, 4),
            "correction":    corr_row["correction"],
            "loco_auc":      corr_row["auc"],
            "loco_ci_lo":    corr_row["ci_lo"],
            "loco_ci_hi":    corr_row["ci_hi"],
            "loco_drop":     round(wc_auc - corr_row["auc"], 4),
        })

    comparison_df = pd.DataFrame(rows)
    comparison_df.to_csv(TABLES_DIR / "auc_comparison_table.csv", index=False)

    print("\nResults summary")

    for model in ["logreg", "lgbm"]:
        print(f"\n  Model: {model}")
        pivot = comparison_df[comparison_df["model"] == model].pivot_table(
            index="test_cohort", columns="correction",
            values="loco_auc", aggfunc="first"
        ).round(4)
        wc_col = comparison_df[
            (comparison_df["model"] == model) &
            (comparison_df["correction"] == "uncorrected")
        ].set_index("test_cohort")["within_auc"]
        pivot.insert(0, "within_cohort", wc_col)
        col_order = (["within_cohort", "uncorrected"] +
                     [c for c in ["combatseq", "mmuphin"]
                      if c in pivot.columns])
        print(pivot[col_order].to_string())

    if not within_df.empty:
        print("\n  Within-cohort AUC on corrected data "
              f"(vs uncorrected Phase 2 baseline):")
        for model in ["logreg", "lgbm"]:
            print(f"    {model}:")
            wc_sub = within_df[within_df["model"] == model]
            for _, row in wc_sub.iterrows():
                ref = within_ref_lookup.get((row["cohort"], model), np.nan)
                delta = row["auc_within_corrected"] - ref
                direction = "▲" if delta > 0 else "▼"
                print(f"      {row['cohort']:20s} {row['correction']:10s} "
                      f"corrected={row['auc_within_corrected']:.4f}  "
                      f"uncorrected={ref:.4f}  "
                      f"Δ={delta:+.4f} {direction}")

    print("\nGenerating figures")

    plot_auc_comparison(comparison_df,
                        FIGURES_DIR / "auc_comparison_bar.png")

    pw_all = pd.concat([
        uncorrected_pw.assign(correction="uncorrected"),
        pairwise_df
    ], ignore_index=True)

    for method in [m for m in CORRECTION_METHODS if
                   (PROCESSED_DIR / f"clr_matrix_{m}.csv").exists()]:
        for model_name in ["logreg", "lgbm"]:
            plot_pairwise_heatmap_corrected(
                pw_all, method, model_name,
                FIGURES_DIR / f"pairwise_heatmap_{method}_{model_name}.png"
            )

    print("\nKey findings")

    for model in ["logreg", "lgbm"]:
        sub = comparison_df[comparison_df["model"] == model]
        print(f"\n  Model: {model}")
        for correction in ["uncorrected"] + list(CORRECTION_METHODS.keys()):
            csub = sub[sub["correction"] == correction]
            if csub.empty:
                continue
            mean_loco = csub["loco_auc"].mean()
            mean_drop = csub["loco_drop"].mean()
            print(f"    {correction:12s}  mean LOCO AUC={mean_loco:.4f}  "
                  f"mean drop={mean_drop:+.4f}")

    print("\nPhase 5 complete.")
    print("Outputs:")
    print(f"  {TABLES_DIR}/loco_auc_corrected.csv")
    print(f"  {TABLES_DIR}/pairwise_auc_corrected.csv")
    print(f"  {TABLES_DIR}/within_cohort_auc_corrected.csv")
    print(f"  {TABLES_DIR}/auc_comparison_table.csv")
    print(f"  {FIGURES_DIR}/auc_comparison_bar.png")
    print("\nNext step: python scripts/08_shap_taxa_comparison.py")


if __name__ == "__main__":
    main()
