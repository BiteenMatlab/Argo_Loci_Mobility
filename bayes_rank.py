import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import pymc as pm
import arviz as az

## run command example:
# python bayes_rank.py /path/to/file.xlsx --draws 2500 --tune 2500 --chains 4
# optional report:
# python bayes_rank.py /path/to/file.xlsx --report_path my_report.txt


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Bayesian hierarchical ranking for wide-format strain columns with replicate label as last column.\n"
            "Model: y = mu_strain + b_rep(rep within strain) + StudentT residual.\n"
            "Assumes each row is a track measurement; melts wide->long and drops NA."
        )
    )
    p.add_argument(
        "file_path",
        type=str,
        help="Path to input data file (.xlsx/.xls or .csv). Replicate column must be the last column.",
    )
    p.add_argument(
        "--draws", type=int, default=2000, help="Posterior draws per chain (default: 2000)."
    )
    p.add_argument(
        "--tune", type=int, default=2000, help="Tuning (warmup) iterations per chain (default: 2000)."
    )
    p.add_argument(
        "--chains", type=int, default=4, help="Number of MCMC chains (default: 4)."
    )
    p.add_argument(
        "--target_accept", type=float, default=0.95, help="NUTS target_accept (default: 0.95)."
    )
    p.add_argument(
        "--seed", type=int, default=109, help="Random seed (default: 109)."
    )
    # Weakly-informative prior scales (on centered log scale)
    p.add_argument("--beta_sd", type=float, default=2.0, help="Prior SD for strain means (default: 2.0).")
    p.add_argument("--sd_rep", type=float, default=0.5, help="HalfNormal scale for replicate SD (default: 0.5).")
    p.add_argument("--sd_res", type=float, default=0.5, help="HalfNormal scale for residual SD (default: 0.5).")
    p.add_argument(
        "--report_path",
        type=str,
        default=None,
        help="Path to save a text report. Default: ./<input_stem>_bayes_report.txt",
    )
    return p.parse_args()


def load_wide_table(path: Path) -> pd.DataFrame:
    if path.suffix.lower() in [".xlsx", ".xls"]:
        return pd.read_excel(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    raise ValueError("Unsupported file type. Use .xlsx/.xls or .csv")


def mode_and_prob(x: np.ndarray) -> tuple[int, float]:
    """Return (mode, probability mass at mode) for integer array."""
    vals, counts = np.unique(x, return_counts=True)
    i = np.argmax(counts)
    return int(vals[i]), float(counts[i] / x.size)


def main():
    args = parse_args()
    path = Path(args.file_path)
    if not path.exists():
        raise FileNotFoundError(f"Cannot find file: {path.resolve()}")

    wide = load_wide_table(path)

    if wide.shape[1] < 2:
        raise ValueError("Expected at least 2 columns: strain columns + replicate column (last).")

    # Replicate column is always the last column
    rep_col = wide.columns[-1]

    # Strain columns are everything except replicate col
    strain_cols = list(wide.columns[:-1])

    # Wide -> long
    long = wide.melt(
        id_vars=[rep_col],
        value_vars=strain_cols,
        var_name="strain",
        value_name="value",
    ).dropna(subset=["value"])

    # Clean labels
    long["replicate"] = long[rep_col].astype(str)
    long = long.drop(columns=[rep_col])

    # Log transform (base 10 log)
    long["value"] = pd.to_numeric(long["value"], errors="coerce")
    long = long.dropna(subset=["value"])
    long["y_raw"] = np.log10(long["value"].to_numpy())

    # Center for stable priors/sampling
    y_mean = float(long["y_raw"].mean())
    long["y"] = long["y_raw"] - y_mean

    # Grouping factor: replicate-within-strain
    long["strain"] = long["strain"].astype("category")
    long["rep_in_strain"] = (
        long["strain"].astype(str) + "__" + long["replicate"].astype(str)
    ).astype("category")

    strain_levels = list(long["strain"].cat.categories)
    S = len(strain_levels)
    strain_codes = long["strain"].cat.codes.to_numpy()

    rep_codes = long["rep_in_strain"].cat.codes.to_numpy()
    R = long["rep_in_strain"].nunique()

    y = long["y"].to_numpy()

    print(f"File: {path.name}")
    print(f"Replicate column (last): {rep_col}")
    print(f"Loaded {len(long):,} observations after melting/dropping NA.")
    print(f"Detected {S} strains and {R} replicate-within-strain groups.")
    print(f"Centering constant (mean log value): {y_mean:.6f}")

    # Bayesian model
    with pm.Model() as model:
        mu_strain = pm.Normal("mu_strain", mu=0.0, sigma=args.beta_sd, shape=S)

        sigma_rep = pm.HalfNormal("sigma_rep", sigma=args.sd_rep)
        b_rep = pm.Normal("b_rep", mu=0.0, sigma=sigma_rep, shape=R)

        sigma = pm.HalfNormal("sigma", sigma=args.sd_res)
        nu_minus1 = pm.Exponential("nu_minus1", lam=1.0)
        nu = pm.Deterministic("nu", nu_minus1 + 1.0)

        mu = mu_strain[strain_codes] + b_rep[rep_codes]

        pm.StudentT("y_obs", nu=nu, mu=mu, sigma=sigma, observed=y)

        idata = pm.sample(
            draws=args.draws,
            tune=args.tune,
            chains=args.chains,
            target_accept=args.target_accept,
            random_seed=args.seed,
            progressbar=True,
        )

    diag_df = az.summary(idata, var_names=["mu_strain", "sigma_rep", "sigma", "nu"], round_to=3)

    print("\nDiagnostics (key parameters):")
    print(diag_df)

    # Posterior draws for strain means: (samples, S)
    mu_draws = idata.posterior["mu_strain"].stack(sample=("chain", "draw")).values.T

    # Rank: 1 = most mobile = largest mu (centered log-scale)
    ranks = np.argsort(np.argsort(-mu_draws, axis=1), axis=1) + 1  # (samples, S)

    rank_mode = np.zeros(S, dtype=int)
    p_rank_mode = np.zeros(S, dtype=float)
    for j in range(S):
        rank_mode[j], p_rank_mode[j] = mode_and_prob(ranks[:, j])

    rank_median = np.median(ranks, axis=0)
    rank_lo = np.quantile(ranks, 0.025, axis=0)
    rank_hi = np.quantile(ranks, 0.975, axis=0)

    rank_table = pd.DataFrame(
        {
            "strain": strain_levels,
            "rank_mode": rank_mode,
            "P(rank_mode)": p_rank_mode,
            "rank_median": rank_median,
            "rank_2.5%": rank_lo,
            "rank_97.5%": rank_hi,
        }
    ).sort_values(["rank_mode", "P(rank_mode)"], ascending=[True, False])

    print("\nPosterior rank summary (1 = most mobile):")
    print(rank_table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    # Pairwise probability matrix P(mu_s > mu_t)
    P = np.zeros((S, S))
    for i in range(S):
        for j in range(S):
            P[i, j] = np.mean(mu_draws[:, i] > mu_draws[:, j])

    P_df = pd.DataFrame(P, index=strain_levels, columns=strain_levels)
    print("\nPairwise probability matrix P(strain row > strain col):")
    print(P_df.round(3).to_string())

    # -----------------------
    # Save report file (requested ordering)
    # -----------------------
    if args.report_path is None:
        report_path = Path.cwd() / f"{path.stem}_bayes_report.txt"
    else:
        report_path = Path(args.report_path)

    report_text = (
        "Posterior rank summary (1 = most mobile):\n"
        f"{rank_table.to_string(index=False, float_format=lambda v: f'{v:.3f}')}\n\n"
        "Pairwise probability matrix P(strain row > strain col):\n"
        f"{P_df.round(3).to_string()}\n\n"
        "----\n"
        "Diagnostics (key parameters):\n"
        f"{diag_df.to_string()}\n"
    )

    report_path.write_text(report_text, encoding="utf-8")
    print(f"\nSaved report to: {report_path.resolve()}")


if __name__ == "__main__":
    main()