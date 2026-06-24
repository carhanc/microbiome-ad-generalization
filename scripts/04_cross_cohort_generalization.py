from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import roc_auc_score
import lightgbm as lgb
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_ROOT / "data" / "processed"
RESULTS_DIR   = PROJECT_ROOT / "results"
TABLES_DIR    = RESULTS_DIR / "tables"
FIGURES_DIR   = RESULTS_DIR / "figures"
MODEL_DIR     = RESULTS_DIR / "model_outputs" / "cross_cohort"

for d in [TABLES_DIR, FIGURES_DIR, MODEL_DIR]:
    d.mkdir(parents=True, exist_ok=True)

SUPERVISED_COHORTS = ["zhuang2018", "ling2020", "shanghai2022", "kazakhstan2022"]

RANDOM_STATE  = 42
INNER_FOLDS   = 5
N_BOOTSTRAP   = 1000


LOGREG_GRID = {"C": [0.001, 0.01, 0.1, 1.0, 10.0]}

LGBM_GRID = {
    "n_estimators": [50, 100, 200],
    "learning_rate": [0.05, 0.1],
    "max_depth":     [3, 5],
    "num_leaves":    [15, 31],
}


def bootstrap_auc_ci(y_true: np.ndarray, y_score: np.ndarray,
                     n_boot: int = N_BOOTSTRAP,
                     alpha: float = 0.05) -> tuple[float, float]:
    rng = np.random.default_rng(RANDOM_STATE)
    n = len(y_true)
    auc_boot = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yt, ys = y_true[idx], y_score[idx]
        if len(np.unique(yt)) < 2:
            continue
        auc_boot.append(roc_auc_score(yt, ys))
    auc_boot = np.array(auc_boot)
    lo = float(np.percentile(auc_boot, 100 * alpha / 2))
    hi = float(np.percentile(auc_boot, 100 * (1 - alpha / 2)))
    return lo, hi


def tune_and_fit(X_train: np.ndarray, y_train: np.ndarray,
                 model_name: str) -> object:
    inner_cv = StratifiedKFold(n_splits=INNER_FOLDS, shuffle=True,
                               random_state=RANDOM_STATE)

    if model_name == "logreg":
        base = LogisticRegression(
            solver="saga",
            l1_ratio=0.5,
            max_iter=10000,
            random_state=RANDOM_STATE,
        )
        grid = LOGREG_GRID

    elif model_name == "lgbm":
        base = lgb.LGBMClassifier(
            objective="binary",
            n_jobs=4,
            random_state=RANDOM_STATE,
            verbose=-1,
        )
        grid = LGBM_GRID

    else:
        raise ValueError(f"Unknown model_name: {model_name}")

    gs = GridSearchCV(
        estimator=base,
        param_grid=grid,
        cv=inner_cv,
        scoring="roc_auc",
        n_jobs=-1,
        refit=True,
    )
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gs.fit(X_train, y_train)

    return gs.best_estimator_, gs.best_params_


def evaluate_on_test(model, X_test: np.ndarray,
                     y_test: np.ndarray) -> dict:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        y_score = model.predict_proba(X_test)[:, 1]

    auc = roc_auc_score(y_test, y_score)
    ci_lo, ci_hi = bootstrap_auc_ci(y_test, y_score)

    return {
        "auc":   round(auc, 4),
        "ci_lo": round(ci_lo, 4),
        "ci_hi": round(ci_hi, 4),
        "n_test": len(y_test),
        "n_AD":  int(y_test.sum()),
        "n_CN":  int((y_test == 0).sum()),
        "y_score": y_score,
    }


def load_data() -> tuple[pd.DataFrame, dict]:
    clr_df  = pd.read_csv(PROCESSED_DIR / "clr_matrix.csv", index_col=0)
    diag_df = pd.read_csv(PROCESSED_DIR / "diagnosis_labels.csv")

    genus_cols = [c for c in clr_df.columns if c != "cohort"]

    clr_reset = clr_df.reset_index().rename(columns={"sample_id": "run_id"})
    merged = clr_reset.merge(
        diag_df[["run_id", "cohort", "diagnosis"]],
        on="run_id", suffixes=("_clr", "")
    )
    if "cohort_clr" in merged.columns:
        merged = merged.drop(columns=["cohort_clr"])
    merged = merged.set_index("run_id")

    merged_binary = merged[merged["diagnosis"].isin(["AD", "CN"])].copy()
    # MCI dropped here — only AD/CN for the cross-cohort stuff

    cohort_data = {}
    for cohort in SUPERVISED_COHORTS:
        sub = merged_binary[merged_binary["cohort"] == cohort]
        X = sub[genus_cols].to_numpy()
        y = (sub["diagnosis"] == "AD").astype(int).to_numpy()
        cohort_data[cohort] = {
            "X": X,
            "y": y,
            "ids": list(sub.index),
            "n": len(y),
            "n_AD": int(y.sum()),
            "n_CN": int((y == 0).sum()),
        }

    return merged_binary, cohort_data, genus_cols


def run_loco(cohort_data: dict) -> list[dict]:
    results = []

    for held_out in SUPERVISED_COHORTS:
        train_cohorts = [c for c in SUPERVISED_COHORTS if c != held_out]
        print(f"\n  LOCO — hold out {held_out}, train on {train_cohorts}")

        X_train = np.vstack([cohort_data[c]["X"] for c in train_cohorts])
        y_train = np.concatenate([cohort_data[c]["y"] for c in train_cohorts])
        X_test  = cohort_data[held_out]["X"]
        y_test  = cohort_data[held_out]["y"]

        n_train_AD = int(y_train.sum())
        n_train_CN = int((y_train == 0).sum())
        print(f"     train: {len(y_train)} ({n_train_AD} AD / {n_train_CN} CN)")
        print(f"     test:  {len(y_test)} ({cohort_data[held_out]['n_AD']} AD / "
              f"{cohort_data[held_out]['n_CN']} CN)")

        for model_name in ["logreg", "lgbm"]:
            print(f"     {model_name}: ", end="", flush=True)
            model, best_params = tune_and_fit(X_train, y_train, model_name)
            metrics = evaluate_on_test(model, X_test, y_test)

            print(f"AUC={metrics['auc']:.4f}  CI=[{metrics['ci_lo']:.4f}, "
                  f"{metrics['ci_hi']:.4f}]  best={best_params}")

            pred_df = pd.DataFrame({
                "run_id":      cohort_data[held_out]["ids"],
                "y_true":      y_test,
                "y_score":     metrics["y_score"],
                "train_cohorts": "+".join(train_cohorts),
            })
            pred_path = MODEL_DIR / f"loco_{held_out}_{model_name}_predictions.csv"
            pred_df.to_csv(pred_path, index=False)

            if model_name == "logreg":
                _save_logreg_coefs(model, list(cohort_data[held_out]["X"].shape),
                                   model_name, f"loco_{held_out}")
            elif model_name == "lgbm":
                _save_lgbm_importances(model, f"loco_{held_out}", model_name)

            results.append({
                "experiment":    "LOCO",
                "train_cohorts": "+".join(train_cohorts),
                "test_cohort":   held_out,
                "model":         model_name,
                "n_train":       len(y_train),
                "n_train_AD":    n_train_AD,
                "n_train_CN":    n_train_CN,
                "n_test":        metrics["n_test"],
                "n_test_AD":     metrics["n_AD"],
                "n_test_CN":     metrics["n_CN"],
                "auc":           metrics["auc"],
                "ci_lo":         metrics["ci_lo"],
                "ci_hi":         metrics["ci_hi"],
                "best_params":   str(best_params),
            })

    return results


def run_pairwise(cohort_data: dict) -> list[dict]:
    results = []

    for train_cohort in SUPERVISED_COHORTS:
        for test_cohort in SUPERVISED_COHORTS:
            if train_cohort == test_cohort:
                continue

            X_train = cohort_data[train_cohort]["X"]
            y_train = cohort_data[train_cohort]["y"]
            X_test  = cohort_data[test_cohort]["X"]
            y_test  = cohort_data[test_cohort]["y"]

            print(f"\n  pairwise: {train_cohort} → {test_cohort}")
            print(f"     train: {len(y_train)} ({cohort_data[train_cohort]['n_AD']} AD / "
                  f"{cohort_data[train_cohort]['n_CN']} CN)")
            print(f"     test:  {len(y_test)} ({cohort_data[test_cohort]['n_AD']} AD / "
                  f"{cohort_data[test_cohort]['n_CN']} CN)")

            for model_name in ["logreg", "lgbm"]:
                print(f"     {model_name}: ", end="", flush=True)
                model, best_params = tune_and_fit(X_train, y_train, model_name)
                metrics = evaluate_on_test(model, X_test, y_test)

                print(f"AUC={metrics['auc']:.4f}  CI=[{metrics['ci_lo']:.4f}, "
                      f"{metrics['ci_hi']:.4f}]  best={best_params}")

                pred_df = pd.DataFrame({
                    "run_id":       cohort_data[test_cohort]["ids"],
                    "y_true":       y_test,
                    "y_score":      metrics["y_score"],
                    "train_cohort": train_cohort,
                })
                pred_path = (MODEL_DIR
                             / f"pairwise_{train_cohort}_to_{test_cohort}"
                               f"_{model_name}_predictions.csv")
                pred_df.to_csv(pred_path, index=False)

                results.append({
                    "experiment":   "pairwise",
                    "train_cohort": train_cohort,
                    "test_cohort":  test_cohort,
                    "model":        model_name,
                    "n_train":      len(y_train),
                    "n_train_AD":   cohort_data[train_cohort]["n_AD"],
                    "n_train_CN":   cohort_data[train_cohort]["n_CN"],
                    "n_test":       metrics["n_test"],
                    "n_test_AD":    metrics["n_AD"],
                    "n_test_CN":    metrics["n_CN"],
                    "auc":          metrics["auc"],
                    "ci_lo":        metrics["ci_lo"],
                    "ci_hi":        metrics["ci_hi"],
                    "best_params":  str(best_params),
                })

    return results


def _save_logreg_coefs(model, shape_info, model_name, prefix):
    coef_path = MODEL_DIR / f"{prefix}_{model_name}_coefs.csv"
    np.savetxt(coef_path, model.coef_, delimiter=",")


def _save_lgbm_importances(model, prefix, model_name):
    imp_path = MODEL_DIR / f"{prefix}_{model_name}_feature_importances.npy"
    np.save(imp_path, model.feature_importances_)


def build_auc_drop_table(loco_results: list, pairwise_results: list,
                         within_cohort_path: Path) -> pd.DataFrame:
    within_df = pd.read_csv(within_cohort_path)
    within_df = within_df[~within_df["cohort"].str.endswith("3class")]
    within_pivot = within_df.set_index(["cohort", "model"])["auc_oof"]

    loco_df     = pd.DataFrame(loco_results)
    pairwise_df = pd.DataFrame(pairwise_results)

    rows = []
    for cohort in SUPERVISED_COHORTS:
        for model in ["logreg", "lgbm"]:
            within_auc = within_pivot.get((cohort, model), np.nan)

            loco_row = loco_df[
                (loco_df["test_cohort"] == cohort) & (loco_df["model"] == model)
            ]
            loco_auc = loco_row["auc"].values[0] if len(loco_row) else np.nan
            loco_ci_lo = loco_row["ci_lo"].values[0] if len(loco_row) else np.nan
            loco_ci_hi = loco_row["ci_hi"].values[0] if len(loco_row) else np.nan

            pw_as_test = pairwise_df[
                (pairwise_df["test_cohort"] == cohort) & (pairwise_df["model"] == model)
            ]
            pw_mean = pw_as_test["auc"].mean() if len(pw_as_test) else np.nan
            pw_min  = pw_as_test["auc"].min()  if len(pw_as_test) else np.nan
            pw_max  = pw_as_test["auc"].max()  if len(pw_as_test) else np.nan

            rows.append({
                "cohort":             cohort,
                "model":              model,
                "within_auc":         round(within_auc, 4),
                "loco_auc":           round(loco_auc, 4),
                "loco_ci_lo":         round(loco_ci_lo, 4),
                "loco_ci_hi":         round(loco_ci_hi, 4),
                "loco_drop":          round(within_auc - loco_auc, 4),
                "pairwise_mean_auc":  round(pw_mean, 4),
                "pairwise_min_auc":   round(pw_min, 4),
                "pairwise_max_auc":   round(pw_max, 4),
            })

    return pd.DataFrame(rows)


def plot_pairwise_heatmap(pairwise_results: list, model_name: str,
                          out_path: Path) -> None:
    df = pd.DataFrame(pairwise_results)
    df = df[df["model"] == model_name]

    cohorts = SUPERVISED_COHORTS
    n = len(cohorts)
    matrix = pd.DataFrame(np.nan, index=cohorts, columns=cohorts)

    for _, row in df.iterrows():
        matrix.loc[row["train_cohort"], row["test_cohort"]] = row["auc"]

    wc = pd.read_csv(TABLES_DIR / "within_cohort_auc.csv")
    wc = wc[(wc["model"] == model_name) & (~wc["cohort"].str.endswith("3class"))]
    for _, row in wc.iterrows():
        matrix.loc[row["cohort"], row["cohort"]] = row["auc_oof"]

    label_map = {
        "zhuang2018":    "Zhuang\n2018\n(CN)",
        "ling2020":      "Ling\n2020\n(CN)",
        "shanghai2022":  "Shanghai\n2022\n(CN)",
        "kazakhstan2022":"Kazakh\n2022\n(KZ)",
    }
    short = [label_map[c] for c in cohorts]

    fig, ax = plt.subplots(figsize=(7, 5.5))

    sns.heatmap(
        matrix.astype(float),
        ax=ax,
        annot=True,
        fmt=".3f",
        cmap="RdYlGn",
        vmin=0.3, vmax=1.0,
        center=0.5,
        linewidths=0.5,
        linecolor="white",
        xticklabels=short,
        yticklabels=short,
        cbar_kws={"label": "AUC-ROC", "shrink": 0.8},
    )

    for i in range(n):
        ax.add_patch(plt.Rectangle((i, i), 1, 1,
                                   fill=False, edgecolor="black",
                                   lw=2.5, clip_on=False))

    model_label = "Logistic Regression" if model_name == "logreg" else "LightGBM"
    ax.set_title(
        f"Cross-Cohort Generalization — {model_label}\n"
        f"Row = train cohort, Column = test cohort | "
        f"Diagonal = within-cohort AUC",
        fontsize=10, pad=12
    )
    ax.set_xlabel("Test cohort", fontsize=10)
    ax.set_ylabel("Train cohort", fontsize=10)
    ax.tick_params(axis="x", labelrotation=0, labelsize=8)
    ax.tick_params(axis="y", labelrotation=0, labelsize=8)

    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_path}")


def plot_loco_bar(loco_results: list, drop_df: pd.DataFrame,
                  out_path: Path) -> None:
    cohorts = SUPERVISED_COHORTS
    cohort_labels = {
        "zhuang2018":    "Zhuang 2018\n(China, n=86)",
        "ling2020":      "Ling 2020\n(China, n=171)",
        "shanghai2022":  "Shanghai 2022\n(China, n=120)",
        "kazakhstan2022":"Kazakhstan 2022\n(KZ, n=84)",
    }

    models   = ["logreg", "lgbm"]
    m_labels = {"logreg": "LogReg", "lgbm": "LightGBM"}
    colors   = {"logreg": "#4C72B0", "lgbm": "#DD8452"}

    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    for ax, model_name in zip(axes, models):
        x = np.arange(len(cohorts))
        width = 0.35

        within_aucs = []
        loco_aucs   = []
        loco_ci_los = []
        loco_ci_his = []

        for cohort in cohorts:
            row = drop_df[(drop_df["cohort"] == cohort) &
                          (drop_df["model"] == model_name)].iloc[0]
            within_aucs.append(row["within_auc"])
            loco_aucs.append(row["loco_auc"])
            loco_ci_los.append(row["loco_auc"] - row["loco_ci_lo"])
            loco_ci_his.append(row["loco_ci_hi"] - row["loco_auc"])

        col = colors[model_name]
        ax.bar(x - width/2, within_aucs,
                        width, label="Within-cohort", color=col, alpha=0.85)
        ax.bar(x + width/2, loco_aucs,
                        width, label="LOCO (cross-cohort)",
                        color=col, alpha=0.45, hatch="///",
                        yerr=[loco_ci_los, loco_ci_his],
                        error_kw={"capsize": 4, "linewidth": 1.2,
                                  "ecolor": "black"})

        ax.axhline(0.5, color="grey", linewidth=0.8, linestyle="--",
                   label="Chance (AUC=0.5)")

        ax.set_xticks(x)
        ax.set_xticklabels([cohort_labels[c] for c in cohorts],
                           fontsize=8, rotation=15, ha="right")
        ax.set_ylim(0.3, 1.05)
        ax.set_ylabel("AUC-ROC", fontsize=10)
        ax.set_title(m_labels[model_name], fontsize=11, fontweight="bold")
        ax.legend(fontsize=8, loc="upper right")
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%.2f"))
        ax.grid(axis="y", alpha=0.3)

    fig.suptitle("Within-Cohort vs LOCO Cross-Cohort AUC\n"
                 "(error bars = 95% bootstrap CI of LOCO AUC)",
                 fontsize=11, y=1.02)
    plt.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {out_path}")


def main():
    print("04_cross_cohort_generalization.py")

    print("\nLoading CLR matrix and diagnosis labels...")
    merged, cohort_data, genus_cols = load_data()

    print("\nCohort summary (binary AD/CN):")
    for cohort in SUPERVISED_COHORTS:
        cd = cohort_data[cohort]
        print(f"  {cohort:20s}  n={cd['n']:3d}  "
              f"({cd['n_AD']} AD / {cd['n_CN']} CN)")

    print("\nLeave-one-cohort-out (LOCO)")
    loco_results = run_loco(cohort_data)

    loco_df = pd.DataFrame(loco_results)
    loco_path = TABLES_DIR / "loco_auc.csv"
    loco_df.drop(columns=["best_params"], errors="ignore").to_csv(loco_path, index=False)
    print(f"\n  wrote {loco_path}")

    print("\nPairwise train/test")
    pairwise_results = run_pairwise(cohort_data)

    pairwise_df = pd.DataFrame(pairwise_results)
    pairwise_path = TABLES_DIR / "pairwise_auc.csv"
    pairwise_df.drop(columns=["best_params"], errors="ignore").to_csv(
        pairwise_path, index=False)
    print(f"\n  wrote {pairwise_path}")

    print("\nAUC drop summary (within-cohort vs LOCO)")
    drop_df = build_auc_drop_table(
        loco_results, pairwise_results,
        TABLES_DIR / "within_cohort_auc.csv"
    )
    drop_path = TABLES_DIR / "auc_drop_summary.csv"
    drop_df.to_csv(drop_path, index=False)

    for model in ["logreg", "lgbm"]:
        print(f"\n  {model}")
        sub = drop_df[drop_df["model"] == model][
            ["cohort", "within_auc", "loco_auc", "loco_drop",
             "pairwise_mean_auc", "pairwise_min_auc"]
        ]
        print(sub.to_string(index=False))

    print(f"\n  wrote {drop_path}")

    print("\nFigures")
    for model_name in ["logreg", "lgbm"]:
        plot_pairwise_heatmap(
            pairwise_results, model_name,
            FIGURES_DIR / f"pairwise_heatmap_{model_name}.png"
        )

    plot_loco_bar(
        loco_results, drop_df,
        FIGURES_DIR / "loco_auc_bar.png"
    )

    print("\nSummary")
    mean_drop = drop_df["loco_drop"].mean()
    max_drop  = drop_df.loc[drop_df["loco_drop"].idxmax()]
    below_chance = drop_df[drop_df["loco_auc"] < 0.55]

    print(f"\n  mean LOCO AUC drop: {mean_drop:.3f}")
    print(f"  largest drop: {max_drop['cohort']} / {max_drop['model']} "
          f"— {max_drop['within_auc']:.3f} → {max_drop['loco_auc']:.3f} "
          f"(Δ={max_drop['loco_drop']:.3f})")

    if len(below_chance) > 0:
        print("\n  LOCO AUC < 0.55:")
        for _, row in below_chance.iterrows():
            print(f"    {row['cohort']:20s} / {row['model']:8s}  "
                  f"AUC={row['loco_auc']:.4f}")
    else:
        print("\n  no LOCO AUC below 0.55")

    print("\nDone.")
    print(f"  {loco_path}")
    print(f"  {pairwise_path}")
    print(f"  {drop_path}")
    print(f"  {FIGURES_DIR}/pairwise_heatmap_*.png")
    print(f"  {FIGURES_DIR}/loco_auc_bar.png")
    print("\nNext: python scripts/05_permanova_variance_decomp.R")


if __name__ == "__main__":
    main()
