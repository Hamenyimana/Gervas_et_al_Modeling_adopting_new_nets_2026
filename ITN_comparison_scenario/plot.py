
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('MacOSX')  # interactive backend on macOS; use 'Agg' for headless
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

COLOR_STD = "#708090"   # slategrey
COLOR_PBO = "#BDB76B"   # darkkhaki
COLOR_IG2 = "#ADD8E6"   # lightblue

MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
# Reporting interval = 30 days; 3 years = 36 periods from ITN deployment

RPT_DAYS = 30
ITN_ABS = 60 // RPT_DAYS  # absolute period of ITN deployment = 2 (March)
N_PERIODS = 35  # 36 months shown, starting from ITN deployment

# uncertainty (30 stochastic realizations)
N_RUNS = 30
Z95 = 1.96
def percent_reduction_with_ci(control_mean, control_sd,
                              intervention_mean, intervention_sd,
                              n=N_RUNS):
    """
    Calculate percentage reduction and approximate 95% CI using
    error propagation.

    Returns
    -------
    reduction : %
    lower95 : %
    upper95 : %
    """

    reduction = (control_mean - intervention_mean) / control_mean * 100.0

    # Standard errors
    se_control = control_sd / np.sqrt(n)
    se_intervention = intervention_sd / np.sqrt(n)

    # Error propagation
    se_reduction = 100 * np.sqrt(
        (intervention_mean / control_mean**2)**2 * se_control**2 +
        (1 / control_mean)**2 * se_intervention**2
    )
    lower = reduction - Z95 * se_reduction
    upper = reduction + Z95 * se_reduction
    return reduction, lower, upper

# Month names starting from March (ITN deployment month)
MONTHS_FROM_MAR = ['Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']

def load_incidence():
    df = pd.read_csv('incidence - 2026-06-30T183427.232.csv')
    # Keep 3 years starting from ITN deployment (March = period ITN_ABS)
    df = df[(df['Day'] >= ITN_ABS) & (df['Day'] < ITN_ABS + N_PERIODS)].copy()
    # Reindex x from 0 (= March year 1)
    df['x'] = df['Day'] - ITN_ABS
    return df
def load_prevalence():
    df = pd.read_csv('Prevalence.csv')
    df['Month'] = (df['Time'] // RPT_DAYS).astype(int)
    # Keep 3 years starting from ITN deployment
    df = df[(df['Month'] >= ITN_ABS) & (df['Month'] < ITN_ABS + N_PERIODS)]
    agg = (df.groupby(['Scenario', 'Month'])
           .agg(RDT=('RDT Prevalence', 'mean'),
                RDT_std=('RDT Prevalence_std', 'mean'))
           .reset_index())
    # Reindex x from 0
    agg['x'] = agg['Month'] - ITN_ABS
    return agg
def _month_labels_from_mar(n=N_PERIODS):
    """Month abbreviations for 36 periods starting from March."""
    return [MONTHS_FROM_MAR[i % 12] for i in range(n)]
def _bar_panel(ax, df, scenarios, colors, ylim, ylabel=None):

    width = 0.25
    for i, (scenario, color) in enumerate(zip(scenarios, colors)):

        sub = (
            df[df["Scenario"] == scenario]
            .sort_values("x")
            .reset_index(drop=True)
        )

        x = sub["x"].values + (i - 1) * width

        mean = sub["5"].clip(*ylim)

        ci95 = 1.96 * sub["5_std"] / np.sqrt(N_RUNS)

        ax.bar(
            x,
            mean,
            width=width,
            color=color,
            edgecolor="black",
            linewidth=0.4
        )
        ax.errorbar(
            x,
            mean,
            yerr=ci95,
            fmt="none",
            ecolor="black",
            elinewidth=0.3,
            capsize=0.5
        )
    ax.set_ylim(ylim)

    ax.set_ylabel(ylabel, fontsize=16)

    ax.set_yticks([0.0, 0.4, 0.8, 1.2])
    ax.tick_params(axis='y', labelsize=16)  # Increase y-axis tick font size
    _decorations(ax)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
def _line_panel(ax, df, scenarios, colors, ylim, ylabel=None):
    """
    Plot monthly RDT prevalence with 95% confidence intervals.
    """
    for scenario, color in zip(scenarios, colors):

        sub = (
            df[df["Scenario"] == scenario]
            .sort_values("x")
            .reset_index(drop=True)
        )

        if sub.empty:
            print(f"No prevalence data found for scenario: {scenario}")
            continue

        x = sub["x"].values
        y = sub["RDT"].values

        ci95 = 1.96 * sub["RDT_std"].values / np.sqrt(N_RUNS)

        ax.plot(
            x,
            y,
            lw=3.5,
            color=color,
            label=scenario
        )

        ax.fill_between(
            x,
            y - ci95,
            y + ci95,
            color=color,
            alpha=0.40
        )

    ax.set_ylim(ylim)

    ax.set_ylabel(ylabel, fontsize=16)

    ax.set_yticks([0.0, 0.2, 0.4, 0.6,0.8])
    ax.tick_params(axis='y', labelsize=16)
    _decorations(ax)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def _decorations(ax):

    labels = _month_labels_from_mar()

    ax.set_xticks(range(N_PERIODS))

    ax.set_xticklabels(
        [labels[i] if i % 2 == 0 else "" for i in range(N_PERIODS)],
        rotation=90,
        fontsize=14
    )

    ax.set_xlim(-0.5, N_PERIODS - 0.5)
def plot_figure(outfile):
    # Load data
    inc = load_incidence()
    prev = load_prevalence()

    inc_ylim = (0.0, 0.8)
    prev_ylim = (0.0, 0.6)

    fig, axes = plt.subplots(
        2,
        1,
        figsize=(16, 10),
        sharex=True
    )

    # (a) clinical incidence
    _bar_panel(
        ax=axes[0],
        df=inc,
        scenarios=["STD", "PBO", "IG2"],
        colors=[COLOR_STD, COLOR_PBO, COLOR_IG2],
        ylim=inc_ylim,

        ylabel="Clinical incidence per child"
    )
    # (b) prevalence
    _line_panel(
        ax=axes[1],
        df=prev,
        scenarios=["STD", "PBO", "IG2"],
        colors=[COLOR_STD, COLOR_PBO, COLOR_IG2],
        ylim=prev_ylim,
        ylabel="RDT prevalence"
    )

    handles = [
        mpatches.Patch(color=COLOR_STD, label="Standard ITNs"),
        mpatches.Patch(color=COLOR_PBO, label="Pyrethroid-PBO ITNs"),
        mpatches.Patch(color=COLOR_IG2, label="Interceptor G2 ITNs"),
    ]

    fig.legend(
        handles=handles,
        loc="lower center",
        ncol=3,
        frameon=False,
        fontsize=12,
        bbox_to_anchor=(0.5, -0.02)
    )
    plt.tight_layout(rect=[0, 0.05, 1, 1])

    os.makedirs("output", exist_ok=True)

    fig.savefig(
        os.path.join("output", "Figure.png"),
        dpi=300,
        bbox_inches="tight"
    )

    fig.savefig(
        os.path.join("output", "Figure.pdf"),
        bbox_inches="tight"
    )

    plt.show()

def main():
    outfile = 'Figure.png'
    plot_figure(outfile)


if __name__ == '__main__':
    main()
