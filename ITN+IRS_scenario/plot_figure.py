
"""
Plot Figure 7: impact of adding IRS to PBO-ITN campaigns.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('MacOSX')  # interactive backend on macOS; use 'Agg' for headless
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D

# constants
COLOR_PBO = '#c8a951'  # gold: PBO ITNs only
COLOR_IRS = '#2b6cb0'  # blue: PBO ITNs + IRS

MONTHS = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

# Reporting interval = 30 days; 3 years = 36 periods from ITN deployment
RPT_DAYS = 30
ITN_ABS = 60 // RPT_DAYS  # absolute period of ITN deployment = 2 (March)
N_PERIODS = 35  # 36 months shown, starting from ITN deployment

# uncertainty (30 stochastic realizations)
N_RUNS = 30
Z95 = 1.96

# percentage reduction with 95% CI
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


# Deployment positions relative to ITN deployment (x=0 = March year 1)
ITN_PERIOD = 0  # ITN at start of figure
IRS_YR2_PERIOD = (500 // RPT_DAYS) - ITN_ABS  # ~June, year 2
IRS_YR3_PERIOD = (865 // RPT_DAYS) - ITN_ABS  # ~June, year 3

# # Year-boundary positions (relative x)
YEAR_BOUNDS = [1]  # end of year 1 and year 2
# Month names starting from March (ITN deployment month)
MONTHS_FROM_MAR = ['Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec', 'Jan', 'Feb']

def find_output_dir():
    here = os.path.dirname(os.path.abspath(__file__))
    d = os.path.join(here, 'output')
    if os.path.isdir(d):
        return d
    raise FileNotFoundError(f"Output dir not found at {d}. Pass path as argument.")

def load_incidence():
    df = pd.read_csv('incidence - 2026-06-21T131630.256.csv')
    # Keep 3 years starting from ITN deployment (March = period ITN_ABS)
    df = df[(df['Day'] >= ITN_ABS) & (df['Day'] < ITN_ABS + N_PERIODS)].copy()
    # Reindex x from 0 (= March year 1)
    df['x'] = df['Day'] - ITN_ABS
    return df

def load_prevalence():
    df = pd.read_csv('ITN_distribution_example_Prevalence (63).csv')
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


def _bar_panel(ax, df, pbo_sc, irs_sc, irs_period, ylim, title, ylabel=None):
    """Grouped monthly bars with 95% confidence intervals."""

    pbo = df[df['Scenario'] == pbo_sc].sort_values('x').reset_index(drop=True)
    irs = df[df['Scenario'] == irs_sc].sort_values('x').reset_index(drop=True)
    x = pbo['x'].values
    w = 0.42

    pbo_mean = pbo['5'].clip(upper=ylim[1])
    irs_mean = irs['5'].clip(upper=ylim[1])

    # 95% CI
    pbo_ci95 = Z95 * pbo['5_std'] / np.sqrt(N_RUNS)
    irs_ci95 = Z95 * irs['5_std'] / np.sqrt(N_RUNS)

    ax.bar(
        x - w / 2,
        pbo_mean,
        width=w,
        color=COLOR_PBO,
        zorder=3
    )


    ax.errorbar(
        x - w / 2,
        pbo_mean,
        yerr=pbo_ci95,
        fmt='none',
        ecolor='0.35',
        elinewidth=0.5,
        capsize=0.6,
        zorder=4
    )

    ax.bar(
        x + w / 2,
        irs_mean,
        width=w,
        color=COLOR_IRS,
        zorder=3
    )

    ax.errorbar(
        x + w / 2,
        irs_mean,
        yerr=irs_ci95,
        fmt='none',
        ecolor='0.35',
        elinewidth=0.5,
        capsize=0.6,
        zorder=4
    )

    _decorations(ax, ylim, irs_period)

    # ax.set_ylim(ylim)
    # ax.set_title(title, fontsize=12)
    ax.set_ylim(ylim)
    # Incidence y-axis ticks
    ax.set_yticks([0.0, 0.4, 0.8])
    ax.set_title(title, fontsize=12)

    if ylabel:
        ax.set_ylabel(ylabel, fontsize=12)

def _decorations(ax, ylim, irs_period):
    """Year boundaries, month x-ticks, and deployment markers."""
    ymax = ylim[1]

    # x-axis: month abbreviations (starting from March) every 2 months
    labels = _month_labels_from_mar()
    tick_pos = list(range(N_PERIODS))
    tick_labels = [labels[i] if i % 2 == 0 else '' for i in range(N_PERIODS)]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_labels, fontsize=12, rotation=90, ha='right')
    ax.set_xlim(-0.5, N_PERIODS - 0.5)

    ax.tick_params(axis='y', labelsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def plot_figure(outfile):
    inc = load_incidence()
    prev = load_prevalence()

    inc_ylim = (0.0, 0.8)
    prev_ylim = (0.0, 0.6)
    fig, axes = plt.subplots(2, 2, figsize=(17, 6), sharex=True, sharey=False)

    _bar_panel(axes[0][0], inc, 'PBO_only', 'PBO_IRS_yr2', IRS_YR2_PERIOD, inc_ylim,
               title='',
               ylabel='Clinical incidence per child')  # (a) IRS in year 2

    _bar_panel(axes[0][1], inc, 'PBO_only', 'PBO_IRS_yr3', IRS_YR3_PERIOD, inc_ylim,
               title='')  # (c) IRS in year 3

    # percentage reduction

    # third year = x = 24-34 (or 24-35 if you later change N_PERIODS=36)
    YEAR3 = np.arange(24, inc["x"].max() + 1)

    # incidence
    inc_year3 = inc[inc["x"].isin(YEAR3)]

    def summary_year3(df, scenario, value_col, sd_col):
        sub = df[df["Scenario"] == scenario]

        mean = sub[value_col].mean()

        # pooled SD over the third year
        sd = np.sqrt(np.mean(sub[sd_col] ** 2))

        return mean, sd

    pbo_mean, pbo_sd = summary_year3(
        inc_year3, "PBO_only", "5", "5_std"
    )

    irs2_mean, irs2_sd = summary_year3(
        inc_year3, "PBO_IRS_yr2", "5", "5_std"
    )

    irs3_mean, irs3_sd = summary_year3(
        inc_year3, "PBO_IRS_yr3", "5", "5_std"
    )

    inc_red2, inc_low2, inc_up2 = percent_reduction_with_ci(
        pbo_mean, pbo_sd,
        irs2_mean, irs2_sd
    )

    inc_red3, inc_low3, inc_up3 = percent_reduction_with_ci(
        pbo_mean, pbo_sd,
        irs3_mean, irs3_sd
    )

    print("\nIncidence reduction (Year 3)")
    print(f"IRS Year2 = {inc_red2:.1f}% (95% CI {inc_low2:.1f}–{inc_up2:.1f})")
    print(f"IRS Year3 = {inc_red3:.1f}% (95% CI {inc_low3:.1f}–{inc_up3:.1f})")

    # prevalence
    prev_year3 = prev[prev["x"].isin(YEAR3)]

    pbo_mean, pbo_sd = summary_year3(
        prev_year3, "PBO_only", "RDT", "RDT_std"
    )

    irs2_mean, irs2_sd = summary_year3(
        prev_year3, "PBO_IRS_yr2", "RDT", "RDT_std"
    )

    irs3_mean, irs3_sd = summary_year3(
        prev_year3, "PBO_IRS_yr3", "RDT", "RDT_std"
    )

    prev_red2, prev_low2, prev_up2 = percent_reduction_with_ci(
        pbo_mean, pbo_sd,
        irs2_mean, irs2_sd
    )

    prev_red3, prev_low3, prev_up3 = percent_reduction_with_ci(
        pbo_mean, pbo_sd,
        irs3_mean, irs3_sd
    )

    print("\nPrevalence reduction (Year 3)")
    print(f"IRS Year2 = {prev_red2:.1f}% (95% CI {prev_low2:.1f}–{prev_up2:.1f})")
    print(f"IRS Year3 = {prev_red3:.1f}% (95% CI {prev_low3:.1f}–{prev_up3:.1f})")


    def _line_panel(ax, df, pbo_sc, irs_sc, irs_period,
                    ylim, title, ylabel=None, xlabel=False):
        """
        Monthly prevalence lines with 95% confidence intervals.
        """

        for sc, color in [(pbo_sc, COLOR_PBO),
                          (irs_sc, COLOR_IRS)]:
            sub = df[df['Scenario'] == sc].sort_values('x')

            x = sub['x'].values
            y = sub['RDT'].clip(*ylim).values

            sd = sub['RDT_std'].values

            # 95% confidence interval
            ci95 = Z95 * sd / np.sqrt(N_RUNS)

            ax.plot(
                x,
                y,
                color=color,
                lw=3,
                zorder=4
            )

            ax.fill_between(
                x,
                np.clip(y - ci95, *ylim),
                np.clip(y + ci95, *ylim),
                color=color,
                alpha=0.18,
                lw=0,
                zorder=3
            )

        _decorations(ax, ylim, irs_period)
        # ax.set_ylim(ylim)
        # ax.set_title(title, fontsize=12)
        ax.set_ylim(ylim)
        # Prevalence y-axis ticks
        ax.set_yticks([0.0, 0.2, 0.4, 0.6])
        ax.set_title(title, fontsize=12)

        if ylabel:
            ax.set_ylabel(ylabel, fontsize=12)

        if xlabel:
            ax.set_xlabel('', fontsize=12)  # Month

    _line_panel(
        axes[1][0],
        prev,
        'PBO_only',
        'PBO_IRS_yr2',
        IRS_YR2_PERIOD,
        prev_ylim,
        title='',
        ylabel='RDT Prevalence',  # Prevalence (PfHRP2, children < 5)
        xlabel=True
    )  # (b) IRS in year 2

    _line_panel(
        axes[1][1],
        prev,
        'PBO_only',
        'PBO_IRS_yr3',
        IRS_YR3_PERIOD,
        prev_ylim,
        title='',
        xlabel=True
    )  # (d) IRS in year 3

    # legend
    handles = [
        mpatches.Patch(color=COLOR_PBO, label='Pyrethroid-PBO ITNs only'),
        mpatches.Patch(color=COLOR_IRS, label='Pyrethroid-PBO ITNs (+IRS)'),
        # Line2D([0], [0], color='0.30', lw=1.2, ls='-', label='ITN deployment'),
        # Line2D([0], [0], color=COLOR_IRS, lw=1.2, ls=':', label='IRS deployment'),
    ]
    fig.legend(handles=handles, loc='lower center', ncol=4, frameon=False,
               fontsize=11, bbox_to_anchor=(0.5, -0.02))

    # bar plot: percentage reduction during year 3
    fig2, ax = plt.subplots(figsize=(7.5, 4.5))

    categories = ["Incidence", "Prevalence"]

    yr2 = [inc_red2, prev_red2]
    yr3 = [inc_red3, prev_red3]

    yr2_err = [
        [inc_red2 - inc_low2, prev_red2 - prev_low2],
        [inc_up2 - inc_red2, prev_up2 - prev_red2]
    ]

    yr3_err = [
        [inc_red3 - inc_low3, prev_red3 - prev_low3],
        [inc_up3 - inc_red3, prev_up3 - prev_red3]
    ]

    x = np.arange(len(categories))
    width = 0.26

    bars1 = ax.bar(
        x - width / 2,
        yr2,
        width,
        color="#4C78A8",
        label="IRS deployed in Year 2",
        # yerr=yr2_err,
        # capsize=4
    )

    bars2 = ax.bar(
        x + width / 2,
        yr3,
        width,
        color="#9ECAE1",
        label="IRS deployed in Year 3",
        # yerr=yr3_err,
        # capsize=4
    )
    ax.set_ylabel("Reduction (%)", fontsize=12)
    ax.set_xticks(x)
    ax.set_xticklabels(categories, fontsize=12)
    ax.set_ylim(0, max(yr2 + yr3) + 15)
    ax.set_yticks([0,20,40, 60])
    ax.legend(frameon=False)

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


    plt.tight_layout()

    fig2.savefig("output/Year3_Reduction.png", dpi=300)
    fig2.savefig("output/Year3_Reduction.pdf")

    plt.show()
    fig.tight_layout(rect=[0, 0.06, 1, 0.97])

    # save figures (high resolution)
    os.makedirs("output", exist_ok=True)

    png_file = os.path.join("output", "Figure.png")
    pdf_file = os.path.join("output", "Figure.pdf")

    fig.savefig(
        png_file,
        dpi=300,
        bbox_inches='tight'
    )

    fig.savefig(
        pdf_file,
        bbox_inches='tight'
    )

    print(f"Saved: {png_file}")
    print(f"Saved: {pdf_file}")

    # Show figure interactively
    plt.show()

    # Close after the window is closed
    plt.close(fig)

    print(f"wrote {outfile}")

    # Show figure interactively
    plt.show()

    # Close after the window is closed
    plt.close(fig)


def main():
    outfile = 'Figure.png'
    plot_figure(outfile)

if __name__ == '__main__':
    main()

