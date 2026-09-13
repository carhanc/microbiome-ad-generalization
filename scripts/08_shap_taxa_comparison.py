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

# Reviewer 3 Round 2, concern #3: mean signed SHAP from a background-relative
# LinearExplainer is not a valid global direction statistic (its sign depends
# on the eval-vs-background feature-mean difference, not just the model
# coefficient). For logistic regression, direction is now defined from the
# fitted model's coefficient sign, tracked per outer fold, with a stability
# criterion PRE-SPECIFIED here before any results are examined:
STABILITY_MIN_FOLDS = 8  # of OUTER_FOLDS=10 folds agreeing in sign to call a
                          # direction "stable" (positive or negative); ties/
                          # minority folds and exact-zero coefficients do not
                          # count toward either direction. Not tuned post hoc.

COEF_ROWS: list[dict] = []  # per-fold logistic-regression coefficients, filled by run_within_cohort_shap


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
    # NOTE: penalty="elasticnet" is never set -- this is L2, l1_ratio is inert
    # (see manuscript Methods 2.4/Limitations; unchanged to avoid an unrequested rerun)
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

            fold_coefs = logreg.coef_.ravel()
            for taxon, coef in zip(genus_cols, fold_coefs):
                COEF_ROWS.append({
                    "cohort": cohort, "fold": fold_idx, "taxon": taxon,
                    "coefficient": float(coef),
                    "sign": "positive" if coef > 0 else ("negative" if coef < 0 else "zero"),
                    "nonzero": bool(coef != 0),
                })

            lgbm_m, lbp = fit_best_lgbm(X_tr, y_tr, g_tr)
            oof_shap_gb[test_idx] = shap_lgbm(lgbm_m, X_te, genus_cols)

            print(f"    fold {fold_idx+1:2d}/{outer_folds}: "
                  f"n_test={len(test_idx)}  logreg C={lp['C']}  "
                  f"lgbm n_est={lbp['n_estimators']}, lr={lbp['learning_rate']}")

        imp_lr = compute_importance(oof_shap_lr, genus_cols, cohort, "logreg", y)
        imp_gb = compute_importance(oof_shap_gb, genus_cols, cohort, "lgbm", y)
        all_rows.append(imp_lr)
        all_rows.append(imp_gb)

        # Console logging only (importance ranking, not a direction claim) --
        # direction is reported elsewhere from the coefficient-stability
        # table (build_coefficient_stability_table), never from mean_shap.
        print(f"    OOF top taxon (logreg): {imp_lr.iloc[0]['taxon']} "
              f"(|SHAP|={imp_lr.iloc[0]['mean_abs_shap']:.4f})")
        print(f"    OOF top taxon (lgbm):   {imp_gb.iloc[0]['taxon']} "
              f"(|SHAP|={imp_gb.iloc[0]['mean_abs_shap']:.4f})")

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


# --- Reviewer 3 Round 2, concern #4 -----------------------------------------
# The observed statistic is the MEAN OF THE SIX UNIQUE PAIRWISE JACCARDS among
# the four real cohorts' top-N sets (C(4,2)=6). The original null instead drew
# only two random sets per replicate (a single pairwise Jaccard), which does
# not match the observed statistic's sampling design -- its mean is still an
# unbiased estimate of the correct null's mean by linearity of expectation,
# but its variance (and therefore its CI and p-value) is not calibrated to a
# mean-of-six statistic. This corrected null draws FOUR independent random
# top-N sets per replicate and averages the same six pairwise Jaccards the
# observed statistic uses, exactly matching its construction.
N_JACCARD_PERM_CORRECTED = 100_000
JACCARD_N_SENSITIVITY = [10, 20, 50]


def four_set_null_replicate(n_genera, top_n, rng):
    sets = [set(rng.choice(n_genera, size=top_n, replace=False)) for _ in range(4)]
    pairwise = [jaccard(sets[a], sets[b])
                for a in range(4) for b in range(a + 1, 4)]
    assert len(pairwise) == 6
    return float(np.mean(pairwise))


def four_set_null_distribution(n_genera, top_n, n_perm, rng):
    return np.array([four_set_null_replicate(n_genera, top_n, rng) for _ in range(n_perm)])


def compute_jaccard_null_baseline_corrected(within_imp, n_genera,
                                            top_n_list=JACCARD_N_SENSITIVITY,
                                            n_perm=N_JACCARD_PERM_CORRECTED):
    rng = np.random.default_rng(RANDOM_STATE)
    rows = []
    for model_name in ["logreg", "lgbm"]:
        for top_n in top_n_list:
            mat, _ = compute_overlap(within_imp, model_name, top_n=top_n)
            off_diag_unique = [mat.loc[c1, c2]
                               for i, c1 in enumerate(SUPERVISED_COHORTS)
                               for c2 in SUPERVISED_COHORTS[i + 1:]]
            assert len(off_diag_unique) == 6, "expected C(4,2)=6 unique cohort pairs"
            observed_mean = float(np.mean(off_diag_unique))

            null_dist = four_set_null_distribution(n_genera, top_n, n_perm, rng)
            mean_null = float(null_dist.mean())
            ci_lo = float(np.percentile(null_dist, 2.5))
            ci_hi = float(np.percentile(null_dist, 97.5))
            # finite-sample empirical one-sided p (Davison & Hinkley, 1997)
            p_value = float((1 + np.sum(null_dist >= observed_mean)) / (n_perm + 1))

            rows.append({
                "model":                 model_name,
                "top_n":                 top_n,
                "n_perm":                n_perm,
                "observed_mean_jaccard_six_pairs": round(observed_mean, 4),
                "mean_null_jaccard":     round(mean_null, 4),
                "ci_lower_null":         round(ci_lo, 4),
                "ci_upper_null":         round(ci_hi, 4),
                "p_value_empirical":     round(p_value, 6),
            })
            where = ("ABOVE" if observed_mean > ci_hi else
                     "BELOW" if observed_mean < ci_lo else "WITHIN")
            print(f"  {model_name} top-{top_n}: observed mean-of-six Jaccard={observed_mean:.4f} | "
                  f"matched null mean={mean_null:.4f} [{ci_lo:.4f}, {ci_hi:.4f}] | "
                  f"observed is {where} the null 95% CI | p={p_value:.6f}")

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


def build_coefficient_stability_table(coef_rows: list[dict]) -> pd.DataFrame:
    """Per (cohort, taxon): fold-level logistic-regression coefficient sign
    stability. sign_stability uses the PRE-SPECIFIED STABILITY_MIN_FOLDS
    threshold (>=8 of 10 outer folds agreeing in sign); zero coefficients
    count toward neither direction."""
    df = pd.DataFrame(coef_rows)
    rows = []
    for (cohort, taxon), g in df.groupby(["cohort", "taxon"]):
        n_total = len(g)
        n_pos = int((g["sign"] == "positive").sum())
        n_neg = int((g["sign"] == "negative").sum())
        n_zero = int((g["sign"] == "zero").sum())
        if n_pos >= STABILITY_MIN_FOLDS:
            stability = "stable_positive"
        elif n_neg >= STABILITY_MIN_FOLDS:
            stability = "stable_negative"
        else:
            stability = "unstable"
        rows.append({
            "cohort": cohort, "taxon": taxon,
            "n_folds_total": n_total,
            "n_positive": n_pos, "n_negative": n_neg, "n_zero": n_zero,
            "median_coefficient": float(g["coefficient"].median()),
            "mean_coefficient": float(g["coefficient"].mean()),
            "sign_stability": stability,
        })
    return pd.DataFrame(rows)


def build_stable_directional_flips(stability_df: pd.DataFrame, within_imp: pd.DataFrame,
                                   top_n: int = TOP_N, min_cohorts: int = MIN_COHORTS) -> pd.DataFrame:
    """Conservative logistic-regression directional-flip definition (Reviewer
    3 Round 2, concern #3): a taxon qualifies only if (a) it appears in the
    top-N OOF mean-|SHAP| ranking in >=min_cohorts cohorts, AND (b) it has a
    STABLE positive coefficient direction in at least one of those cohorts,
    AND (c) a STABLE negative coefficient direction in at least one other.
    Unstable per-cohort directions do not count as evidence either way."""
    top_sets = {}
    for cohort in SUPERVISED_COHORTS:
        sub = within_imp[(within_imp["cohort"] == cohort) & (within_imp["model"] == "logreg")].head(top_n)
        top_sets[cohort] = set(sub["taxon"])

    taxon_cohort_count = {}
    for cohort, s in top_sets.items():
        for t in s:
            taxon_cohort_count[t] = taxon_cohort_count.get(t, 0) + 1
    candidate_taxa = {t for t, n in taxon_cohort_count.items() if n >= min_cohorts}

    rows = []
    for taxon in sorted(candidate_taxa):
        sub = stability_df[stability_df["taxon"] == taxon]
        stable_pos_cohorts = sorted(sub[sub["sign_stability"] == "stable_positive"]["cohort"])
        stable_neg_cohorts = sorted(sub[sub["sign_stability"] == "stable_negative"]["cohort"])
        is_stable_flip = len(stable_pos_cohorts) >= 1 and len(stable_neg_cohorts) >= 1
        rows.append({
            "taxon": taxon,
            "n_cohorts_in_top_n": taxon_cohort_count[taxon],
            "stable_positive_cohorts": "+".join(stable_pos_cohorts) if stable_pos_cohorts else "",
            "stable_negative_cohorts": "+".join(stable_neg_cohorts) if stable_neg_cohorts else "",
            "is_stable_directional_flip": is_stable_flip,
        })
    return pd.DataFrame(rows)


def plot_top15_per_cohort(importance_df, model_name, out_path, coef_stability=None):
    model_label = "Logistic Regression" if model_name == "logreg" else "LightGBM"
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    axes = axes.flatten()

    for ax, cohort in zip(axes, SUPERVISED_COHORTS):
        sub = importance_df[
            (importance_df["cohort"] == cohort) &
            (importance_df["model"] == model_name)
        ].head(TOP_DISPLAY).iloc[::-1]

        if model_name == "logreg" and coef_stability is not None:
            # Direction comes from the fitted model's coefficient (median
            # across outer folds), not from mean signed SHAP.
            coef_lookup = coef_stability[coef_stability["cohort"] == cohort] \
                .set_index("taxon")["median_coefficient"]
            colours = ["#D62728" if coef_lookup.get(t, 0) > 0 else "#1F77B4"
                       for t in sub["taxon"]]
            legend_handles = [
                mpatches.Patch(color="#D62728", label="Positive model coefficient"),
                mpatches.Patch(color="#1F77B4", label="Negative model coefficient"),
            ]
        else:
            # LightGBM has no single coefficient-like direction; a tree
            # ensemble's SHAP relationship can be nonlinear/nonmonotonic, so
            # no biological/predictive direction is assigned here at all
            # (Reviewer 3 Round 2, concern #3) -- bars are feature-importance
            # magnitude only.
            colours = ["#4C72B0"] * len(sub)
            legend_handles = None

        ax.barh(sub["taxon"], sub["mean_abs_shap"], color=colours, alpha=0.85)
        ax.set_xlabel("Mean |SHAP value| (CLR units)", fontsize=9)
        ax.set_title(f"{cohort}\n({COHORT_LABELS[cohort].replace(chr(10), ' ')})",
                     fontsize=9, fontweight="bold",
                     color=COHORT_COLOURS[cohort])
        ax.tick_params(axis="y", labelsize=8)
        ax.tick_params(axis="x", labelsize=8)
        ax.grid(axis="x", alpha=0.3)

        if legend_handles is not None:
            ax.legend(handles=legend_handles, fontsize=7, loc="lower right")

    if model_name == "logreg":
        subtitle = ("Red = positive fitted-model coefficient (higher CLR abundance predicts "
                    "higher AD log-odds), Blue = negative coefficient")
    else:
        subtitle = ("Bars show OOF mean |SHAP| feature importance only; LightGBM's "
                    "tree-ensemble relationships are not assigned a single AD/CN "
                    "predictive direction (see Methods)")
    fig.suptitle(
        f"Top {TOP_DISPLAY} Taxa by SHAP Importance — {model_label}\n({subtitle})",
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


# plot_direction_dotplot() was removed as part of the Round 2 publication-
# consistency patch: it plotted mean signed SHAP on its x-axis and labeled
# the axes/title "Mean SHAP value (positive -> AD, negative -> CN)",
# "Cross-Cohort SHAP Direction", "<- CN-associated" / "AD-associated ->" --
# exactly the invalid direction statistic Reviewer 3 concern #3 identified,
# even though its star/highlight marker had already been switched to the
# coefficient-stability definition. Rather than patch the axis semantics of
# an auxiliary, non-manuscript figure, this function was deleted outright:
# the manuscript's actual Figure 6 Panel B (generate_manuscript_figures.py::
# make_fig6) now shows the same eight-genus coefficient-stability screen as
# a coefficient-value grid (rows = taxa, columns = cohorts, cell = median
# fitted coefficient + fold-count stability annotation), which supersedes
# what this dot plot was for. No coefficient-like mean-SHAP-as-direction
# visualization remains anywhere in this script.


# plot_loco_direction() was removed in Round 2 (Reviewer 3, concern #3): it
# plotted mean signed SHAP from the LOCO explainer as if it were a global
# AD/CN direction, but LOCO's SHAP background (the N-1 training cohorts) and
# evaluation set (the held-out cohort) come from systematically different
# populations, so its sign is confounded with raw cross-cohort abundance
# differences rather than reflecting the model's learned direction. LOCO
# feature-IMPORTANCE (mean |SHAP|, magnitude only) remains available in
# results/tables/shap_loco_importance.csv and is not affected by this issue.


def main():
    print("08_shap_taxa_comparison.py — Phase 6 (Round 2: coefficient-stability direction, matched Jaccard null)")

    cohort_data, genus_cols = load_data()
    print(f"\nLoaded {len(genus_cols)} genera, {len(SUPERVISED_COHORTS)} cohorts")

    print("\nWithin-cohort SHAP")
    within_imp = run_within_cohort_shap(cohort_data, genus_cols)
    within_imp.to_csv(TABLES_DIR / "shap_within_cohort_importance.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/shap_within_cohort_importance.csv")

    print("\nLogistic-regression coefficient stability "
          f"(pre-specified threshold: >={STABILITY_MIN_FOLDS}/{OUTER_FOLDS} folds agreeing in sign)")
    coef_stability_df = build_coefficient_stability_table(COEF_ROWS)
    coef_stability_df.to_csv(TABLES_DIR / "logreg_coefficient_stability.csv", index=False)
    print(f"  Wrote: {TABLES_DIR}/logreg_coefficient_stability.csv "
          f"({len(coef_stability_df)} cohort x genus rows)")
    print(coef_stability_df["sign_stability"].value_counts().to_string())

    print("\nLOCO SHAP (train on N-1 cohorts, SHAP on held-out) — "
          "feature-importance ranking only; no cross-cohort AD/CN direction is assigned "
          "from LOCO SHAP (background/eval come from different cohorts, confounding sign "
          "with raw cross-cohort abundance differences; Reviewer 3 Round 2, concern #3)")
    loco_imp = run_loco_shap(cohort_data, genus_cols)
    loco_imp.to_csv(TABLES_DIR / "shap_loco_importance.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/shap_loco_importance.csv")

    print(f"\nOverlap analysis (top-N per cohort, magnitude only — unaffected by the direction issue)")

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

    print(f"\nMATCHED Jaccard null permutation baseline (Reviewer 3 Round 2, concern #4): "
          f"{N_JACCARD_PERM_CORRECTED:,} replicates, each drawing 4 random top-N sets and "
          f"averaging the same 6 unique pairwise Jaccards the observed statistic uses; "
          f"N sensitivity = {JACCARD_N_SENSITIVITY}")
    null_df_corrected = compute_jaccard_null_baseline_corrected(within_imp, len(genus_cols))
    null_df_corrected.to_csv(TABLES_DIR / "jaccard_null_corrected.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/jaccard_null_corrected.csv "
          f"(supersedes the old jaccard_null.csv, which used a mismatched 2-set null and "
          f"is retained on disk only for audit-trail provenance, not used by the manuscript)")

    print(f"\nDirectional SHAP analysis (shared top-{TOP_N} taxa) — descriptive overlap table, "
          "both models; the FLIP CLAIM itself is now reported ONLY for logistic regression, "
          "using the coefficient-stability definition, not mean-signed-SHAP")

    dir_rows = []
    for model_name in ["logreg", "lgbm"]:
        dir_df = directional_analysis(within_imp, model_name)
        dir_df["model"] = model_name
        dir_rows.append(dir_df)
        shared_count = dir_df["taxon"].nunique()
        print(f"  {model_name}: {shared_count} taxa shared across ≥{MIN_COHORTS} cohorts' top-{TOP_N} "
              f"(descriptive overlap only, no direction/flip claim{'' if model_name == 'logreg' else ' for lgbm — see below'})")

    dir_all = pd.concat(dir_rows, ignore_index=True)
    dir_all.to_csv(TABLES_DIR / "shap_directional_flips.csv", index=False)
    print(f"\n  Wrote: {TABLES_DIR}/shap_directional_flips.csv "
          f"(descriptive shared-taxa/mean-SHAP table for both models; RETAINED for provenance "
          f"but its old mean-signed-SHAP flip count is superseded by "
          f"logreg_directional_flips_stable.csv below and is not used for the LightGBM claim at all)")

    print("\nStable directional-flip taxa (logistic regression only)")
    stable_flip_df = build_stable_directional_flips(coef_stability_df, within_imp)
    stable_flip_df.to_csv(TABLES_DIR / "logreg_directional_flips_stable.csv", index=False)
    n_stable_flips = int(stable_flip_df["is_stable_directional_flip"].sum())
    stable_flip_taxa = set(stable_flip_df[stable_flip_df["is_stable_directional_flip"]]["taxon"])
    print(f"  Wrote: {TABLES_DIR}/logreg_directional_flips_stable.csv")
    print(f"  {n_stable_flips} taxa meet the pre-specified stable-flip criterion "
          f"(top-{TOP_N} in >={MIN_COHORTS} cohorts AND a stable positive coefficient in >=1 "
          f"cohort AND a stable negative coefficient in >=1 other cohort):")
    if stable_flip_taxa:
        print(f"    {sorted(stable_flip_taxa)}")

    print("\nCONFIRMATION: no LightGBM directional-flip claim is computed or reported anywhere "
          "in this script (tree ensembles have no single coefficient-like monotonic direction; "
          "SHAP relationships may be nonlinear/nonmonotonic). LightGBM retains only OOF "
          "mean-|SHAP| feature importance, top-feature overlap, and Jaccard analysis.")

    print("\nGenerating figures")

    plot_top15_per_cohort(
        within_imp, "logreg",
        FIGURES_DIR / "shap_top15_per_cohort_logreg.png",
        coef_stability=coef_stability_df,
    )
    plot_top15_per_cohort(
        within_imp, "lgbm",
        FIGURES_DIR / "shap_top15_per_cohort_lgbm.png",
    )
    plot_jaccard_heatmap(
        within_imp,
        FIGURES_DIR / "shap_overlap_jaccard.png"
    )

    # No dot-plot direction figure and no plot_loco_direction call for either
    # model -- both would assign a global AD/CN direction from mean signed
    # SHAP, which Reviewer 3 Round 2 concern #3 establishes is not valid
    # (LightGBM: no coefficient exists at all; LOCO: background/eval cohort
    # mismatch confounds sign with raw cross-cohort abundance shifts for both
    # models). The coefficient-stability screen (stable_flip_taxa) is
    # visualized in the manuscript's Figure 6 Panel B (generate_manuscript_
    # figures.py::make_fig6) as a coefficient-value grid, not as a dot plot.

    print("\nKey findings")

    for model_name in ["logreg", "lgbm"]:
        mat, top_sets = compute_overlap(within_imp, model_name, top_n=TOP_N)
        off = [mat.loc[c1, c2]
               for c1 in SUPERVISED_COHORTS
               for c2 in SUPERVISED_COHORTS if c1 != c2]
        dir_sub = dir_all[dir_all["model"] == model_name]
        n_shared = dir_sub["taxon"].nunique()

        print(f"\n  {model_name}")
        print(f"    Mean pairwise Jaccard (top-{TOP_N}): {np.mean(off):.3f}")
        print(f"    Taxa shared ≥{MIN_COHORTS} cohorts: {n_shared}")
        if model_name == "logreg":
            print(f"    Stable coefficient-based direction flips: {n_stable_flips}")
            if stable_flip_taxa:
                print(f"    Flipped taxa: {sorted(stable_flip_taxa)}")
        else:
            print("    Direction flips: not assessed for LightGBM (no coefficient-like "
                  "direction exists for a tree ensemble; see Methods)")

        print(f"\n    Top-3 taxa per cohort ({model_name}, ranked by mean |SHAP|; "
              f"no direction implied for lgbm):")
        for cohort in SUPERVISED_COHORTS:
            top3 = within_imp[
                (within_imp["cohort"] == cohort) &
                (within_imp["model"] == model_name)
            ].head(3)[["taxon", "mean_abs_shap"]].copy()
            lines = [f"{r['taxon']} (|SHAP|={r['mean_abs_shap']:.3f})"
                     for _, r in top3.iterrows()]
            print(f"      {cohort}: {' | '.join(lines)}")

    print("\nPhase 6 complete (Round 2).")
    print("Outputs:")
    for f in [
        "shap_within_cohort_importance.csv", "logreg_coefficient_stability.csv",
        "shap_loco_importance.csv", "shap_taxa_overlap.csv",
        "jaccard_null_corrected.csv", "shap_directional_flips.csv",
        "logreg_directional_flips_stable.csv",
    ]:
        print(f"  {TABLES_DIR}/{f}")
    for f in [
        "shap_top15_per_cohort_logreg.png", "shap_top15_per_cohort_lgbm.png",
        "shap_overlap_jaccard.png",
    ]:
        print(f"  {FIGURES_DIR}/{f}")
    print("\nNext step: manuscript draft (manuscript/draft.md)")


if __name__ == "__main__":
    main()
