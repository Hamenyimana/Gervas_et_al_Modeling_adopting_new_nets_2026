"""
Plot Figure 8: all three resistance loci shown continuously.

Each panel (scenario x species) shows genotype fractions for loci a, b, and c
over the full 9 years, distinguished by line style:
    locus a (pyrethroid / STD)   - solid
    locus b (PBO)                - dashed
    locus c (IG2 / chlorfenapyr) - dash-dot

Genotype colour:
    homo_susceptible  dark navy
    heterozygous      light blue
    homo_resistant    orange

Background shading shows which net is active in each 3-year block.
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('MacOSX')  # interactive backend on macOS; use 'Agg' for headless
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# constants
SCENARIOS = {
    'STD-STD-STD': ['STD', 'STD', 'STD'],
    'STD-PBO-STD': ['STD', 'PBO', 'STD'],
    'STD-PBO-PBO': ['STD', 'PBO', 'PBO'],
    'STD-PBO-IG2': ['STD', 'PBO', 'IG2'],
}
NET_LOCUS = {'STD': 'a', 'PBO': 'b', 'IG2': 'c'}
SCENARIO_TITLE = {
    'STD-STD-STD': '(a) STD → STD → STD',
    'STD-PBO-STD': '(b) STD → PBO → STD',
    'STD-PBO-PBO': '(c) STD → PBO → PBO',
    'STD-PBO-IG2': '(d) STD → PBO → IG2',
}
SPECIES = ['funestus', 'arabiensis']
LOCI = ['a', 'b', 'c']
LOCUS_LABEL = {'a': 'Locus a (pyrethroid)', 'b': 'Locus b (PBO)', 'c': 'Locus c (IG2)'}
BLOCK_YEARS = 3
N_BLOCKS = 3

# uncertainty (30 stochastic realizations)
N_RUNS = 30
Z95 = 1.96

# Genotype colours
GENO_ORDER = ['homo_susceptible', 'heterozygous', 'homo_resistant']
GENO_COLOR = {
    'homo_susceptible': '#2d3561',
    'heterozygous': '#6db1d9',
    'homo_resistant': '#e07b39',
}
GENO_LABEL = {
    'homo_susceptible': r'Homo. susceptible ($x_0$:$x_0$)',
    'heterozygous': r'Heterozygous ($x_0$:$x_1$)',
    'homo_resistant': r'Homo. resistant ($x_1$:$x_1$)',
}

# Line style per locus
LOCUS_LS = {'a': '-', 'b': '--', 'c': '-.'}
LOCUS_LW = {'a': 2.6, 'b': 2.6, 'c': 2.6}
# Background shading colour per net type (very light)
NET_BG = {'STD': '#f0f0f0', 'PBO': '#e8f4ea', 'IG2': '#eee8f5'}


def load(path):
    df = pd.read_csv(path)
    df['Year'] = df['Time'] / 365.0
    return df


# final homozygous-resistant summary
def summarize_final_homozygous(df, outfile):

    final_time = df["Time"].max()

    summary = (
        df[
            (df["Time"] == final_time) &
            (df["Genotype"] == "homo_resistant")
        ]
        .copy()
    )

    summary["CI95"] = (
        1.96 *
        summary["Fraction_std"] /
        np.sqrt(N_RUNS)
    )

    summary["Percent"] = 100 * summary["Fraction"]
    summary["Lower95"] = 100 * (summary["Fraction"] - summary["CI95"])
    summary["Upper95"] = 100 * (summary["Fraction"] + summary["CI95"])

    summary["Lower95"] = summary["Lower95"].clip(lower=0)
    summary["Upper95"] = summary["Upper95"].clip(upper=100)

    print("\n")
    print("=" * 80)
    print("Final homozygous resistant genotype frequencies")
    print("=" * 80)

    for _, r in summary.iterrows():

        species = (
            "An. funestus"
            if r["Species"] == "funestus"
            else "An. arabiensis"
        )

        print(
            f"{r['Scenario']:12s} | "
            f"{species:18s} | "
            f"Locus {r['Locus']} | "
            f"{r['Percent']:.1f}% "
            f"(95% CI {r['Lower95']:.1f}–{r['Upper95']:.1f})"
        )

    summary.to_csv(outfile, index=False)


    return summary


def _plot_locus(ax, data, locus):
    """
    Draw genotype fractions for one locus with 95% confidence intervals.
    """
    ls = LOCUS_LS[locus]
    lw = LOCUS_LW[locus]
    for geno in GENO_ORDER:

        g = (
            data[(data['Locus'] == locus) &
                 (data['Genotype'] == geno)]
            .sort_values('Year')
        )

        if g.empty:
            continue

        color = GENO_COLOR[geno]

        ax.plot(
            g['Year'],
            g['Fraction'],
            color=color,
            ls=ls,
            lw=lw,
            zorder=3
        )

        # 95% confidence intervals
        if 'Fraction_std' in g.columns and g['Fraction_std'].notna().any():
            ci95 = (
                    Z95 *
                    g['Fraction_std'] /
                    np.sqrt(N_RUNS)
            )

            lo = (
                    g['Fraction'] - ci95
            ).clip(lower=0)

            hi = (
                    g['Fraction'] + ci95
            ).clip(upper=1)

            ax.fill_between(
                g['Year'],
                lo,
                hi,
                color=color,
                alpha=0.15,
                lw=0,
                zorder=2
            )


def plot_figure8(df, outfile):
    scenarios = [s for s in SCENARIOS if s in df['Scenario'].unique()]
    n_col = len(scenarios)
    n_row = len(SPECIES)
    fig, axes = plt.subplots(n_row, n_col,
                             figsize=(3.8 * n_col, 2.8 * n_row + 0.6),
                             sharex=True, sharey=True)
    axes = np.atleast_2d(axes)
    summary_file = os.path.join(
        "output",
        "Final_homozygous_resistant_summary.csv"
    )

    summarize_final_homozygous(df, summary_file)
    for r, species in enumerate(SPECIES):
        for c, scenario in enumerate(scenarios):
            ax = axes[r][c]
            nets = SCENARIOS[scenario]
            dsub = df[(df['Species'] == species) & (df['Scenario'] == scenario)]
            # background shading: one shade per 3-year block
            for block, net in enumerate(nets):
                x0, x1 = block * BLOCK_YEARS, (block + 1) * BLOCK_YEARS
                ax.axvspan(x0, x1, color=NET_BG[net], alpha=1.0, zorder=0, lw=0)
                ax.text((x0 + x1) / 2, 0.97, net, ha='center', va='top',
                        fontsize=7.5, color='0.35', transform=ax.get_xaxis_transform())

            # Block boundary lines
            for block in range(1, N_BLOCKS):
                ax.axvline(block * BLOCK_YEARS, color='0.55', lw=0.25, ls='--', zorder=1)

            # All three loci, continuously over 9 years
            for locus in LOCI:
                _plot_locus(ax, dsub, locus)

            ax.set_xlim(0, N_BLOCKS * BLOCK_YEARS)
            ax.set_ylim(0, 1)
            ax.set_xticks(range(0, N_BLOCKS * BLOCK_YEARS + 1, 3))
            ax.tick_params(labelsize=8)

            if r == 0:
                ax.set_title(SCENARIO_TITLE.get(scenario, scenario), fontsize=9.5, pad=4)
            if c == 0:
                sp_italic = r'$\it{An.\ funestus}$' if species == 'funestus' else r'$\it{An.\ arabiensis}$'
                ax.set_ylabel(sp_italic + '\nFraction of vector population', fontsize=8.5)
            if r == n_row - 1:
                ax.set_xlabel('Years', fontsize=8.5)

    # legend
    # Row 1: genotype colours (solid line)
    geno_handles = [
        Line2D([0], [0], color=GENO_COLOR[g], lw=3.0, ls='-', label=GENO_LABEL[g])
        for g in GENO_ORDER
    ]
    # Row 2: locus line styles (shown in neutral grey)
    locus_handles = [
        Line2D([0], [0], color='0.3', lw=2, ls=LOCUS_LS[l], label=LOCUS_LABEL[l])
        for l in LOCI
    ]
    # Row 3: background net colours
    bg_handles = [
        Patch(facecolor=NET_BG[n], edgecolor='0.6', lw=0.5, label=f'{n} net active')
        for n in ['STD', 'PBO', 'IG2']
    ]

    all_handles = geno_handles + locus_handles + bg_handles
    fig.legend(handles=all_handles, loc='lower center', ncol=3,
               frameon=False, fontsize=8, columnspacing=1.2,
               bbox_to_anchor=(0.5, 0.0))

    # fig.suptitle(
    #     'Figure 8: Vector resistance genotype dynamics by ITN strategy '
    #     '(resistance probability 0.75)\n'
    #     '(Lines show means from 30 stochastic simulations; '
    #     'shaded regions indicate 95% confidence intervals)',
    #     fontsize=10.5,
    #     y=1.02
    # )
    fig.tight_layout(rect=[0, 0.13, 1, 1.0])
    fig.savefig(outfile, dpi=300, bbox_inches='tight')
    fig.savefig(os.path.splitext(outfile)[0] + '.pdf', bbox_inches='tight')
    print(f"Saved: {outfile}")
    print(f"Saved: {os.path.splitext(outfile)[0] + '.pdf'}")
    plt.show()
    plt.close(fig)

def main():
    csv_file = "Genotype_fractions.csv"

    print(f"reading {csv_file}")

    df = load(csv_file)

    if 'Resistance' in df.columns:
        vals = df['Resistance'].astype(float)
        if (vals == 0.75).any():
            df = df[np.isclose(vals, 0.75)]

    os.makedirs("output", exist_ok=True)

    outfile = os.path.join("output", "Figure.png")

    plot_figure8(df, outfile)


if __name__ == '__main__':
    main()



