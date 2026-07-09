"""
VectorGeneticsAnalyzer

Reads the ReportVectorGenetics GENOME output for each species and converts the
per-genome counts into per-locus genotype fractions (homozygous susceptible,
heterozygous, homozygous resistant) over time.

The reporter writes one CSV per species named
ReportVectorGenetics_<species>_Female_GENOME.csv with columns:
    Time, NodeID, Genome, VectorPopulation

Genome strings look like 'X-a0-b0-c0:X-a1-b0-c0':
  - two gametes separated by ':'  (diploid)
  - each gamete = <gender>-<locus a>-<locus b>-<locus c>, tokens separated by '-'
The locus order (a, b, c) is fixed by the order genes are added in
helpers.set_param_fn, so token index 1 -> locus a, 2 -> locus b, 3 -> locus c
(token 0 is the gender X/Y locus).

Output (in <working_dir>/output/):
  <expt>_genotype_fractions_full.csv  - per-realization tidy data
  <expt>_genotype_fractions.csv       - mean (+ std) across realizations
"""
import os
import pandas as pd
import numpy as np
from typing import Dict, Any, Union
from idmtools.entities.ianalyzer import IAnalyzer as BaseAnalyzer

import matplotlib as mpl
from idmtools.entities.iworkflow_item import IWorkflowItem
from idmtools.entities.simulation import Simulation

from logging import getLogger

from idmtools.analysis.analyze_manager import AnalyzeManager
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform

mpl.use('Agg')

GENOTYPE_ORDER = ['homo_susceptible', 'heterozygous', 'homo_resistant']


def genome_to_genotype_fractions(df, loci=('a', 'b', 'c')):
    """Convert a GENOME-stratified report into per-locus genotype fractions over time.

    Pure pandas/numpy (no idmtools dependency) so it can be unit-tested directly.

    Args:
        df: DataFrame with columns Time, Genome, VectorPopulation. Genome strings are
            '<gam1>:<gam2>' where each gamete is '<gender>-<locusA>-<locusB>-...'.
        loci: locus letters in the order genes were added (token order in the genome).

    Returns:
        DataFrame with columns Time, Locus, Genotype, VectorPopulation, Fraction;
        every (Time, Locus, Genotype) combination present (absent ones filled with 0)
        so the three genotype fractions sum to 1 per (Time, Locus).
    """
    loci = list(loci)
    df = df[['Time', 'Genome', 'VectorPopulation']].copy()
    gametes = df['Genome'].str.split(':', expand=True)
    gam1 = gametes[0].str.split('-', expand=True)
    gam2 = gametes[1].str.split('-', expand=True)

    frames = []
    for i, locus in enumerate(loci):
        token = i + 1  # token 0 is the gender (X/Y) locus
        n_resistant = ((gam1[token] == f'{locus}1').astype(int)
                       + (gam2[token] == f'{locus}1').astype(int))
        genotype = n_resistant.map({0: 'homo_susceptible', 1: 'heterozygous', 2: 'homo_resistant'})
        frames.append(pd.DataFrame({'Time': df['Time'].values,
                                    'Locus': locus,
                                    'Genotype': genotype.values,
                                    'VectorPopulation': df['VectorPopulation'].values}))
    long = pd.concat(frames, ignore_index=True)
    summed = long.groupby(['Time', 'Locus', 'Genotype'], as_index=False)['VectorPopulation'].sum()

    # Ensure every (Time, Locus, Genotype) combination exists so the three genotype
    # fractions always sum to 1 (fill absent genotypes with 0).
    times = sorted(summed['Time'].unique())
    full_idx = pd.MultiIndex.from_product([times, loci, GENOTYPE_ORDER],
                                          names=['Time', 'Locus', 'Genotype'])
    summed = (summed.set_index(['Time', 'Locus', 'Genotype'])
                    .reindex(full_idx, fill_value=0)
                    .reset_index())
    totals = summed.groupby(['Time', 'Locus'])['VectorPopulation'].transform('sum')
    summed['Fraction'] = summed['VectorPopulation'] / totals.replace(0, np.nan)
    return summed


class VectorGeneticsAnalyzer(BaseAnalyzer):

    def __init__(self, species=('funestus', 'arabiensis'), loci=('a', 'b', 'c'),
                 expt_name='ITN_distribution_example', tags=None, title='idm'):
        self.species = list(species)
        self.loci = list(loci)
        self.expt_name = expt_name
        filenames = [f"output\\ReportVectorGenetics_{sp}_Female_GENOME.csv" for sp in self.species]
        super().__init__(filenames=filenames)
        # Tags used to group/average across stochastic realizations (Run_Number is
        # intentionally excluded so it is averaged over).
        self.tags = tags if tags is not None else [
            'Scenario', 'Resistance', 'Larval_Capacity_funestus', 'Larval_Capacity_arabiensis',
        ]
        print(title)

    def initialize(self):
        if not os.path.exists(os.path.join(self.working_dir, "output")):
            os.mkdir(os.path.join(self.working_dir, "output"))

    def map(self, data: Dict[str, Any], item: Union[IWorkflowItem, Simulation]) -> Any:
        per_species = []
        for sp, fn in zip(self.species, self.filenames):
            g = genome_to_genotype_fractions(data[fn], loci=self.loci)
            g['Species'] = sp
            per_species.append(g)
        return pd.concat(per_species, ignore_index=True)

    def reduce(self, all_data: Dict[Union[IWorkflowItem, Simulation], Any]) -> Any:
        output_dir = os.path.join(self.working_dir, "output")
        df = pd.DataFrame()
        for s, v in all_data.items():
            tmp = v.copy()
            for t in self.tags:
                tmp[t] = s.tags[t]
            df = pd.concat([df, tmp], ignore_index=True)
        df.to_csv(os.path.join(output_dir, f"{self.expt_name}_genotype_fractions_full.csv"), index=False)

        group_cols = self.tags + ['Species', 'Locus', 'Genotype', 'Time']
        agg = df.groupby(group_cols, as_index=False).agg(
            Fraction=('Fraction', 'mean'),
            Fraction_std=('Fraction', 'std'),
        )
        agg.to_csv(os.path.join(output_dir, f"{self.expt_name}_genotype_fractions.csv"), index=False)


if __name__ == '__main__':
    logger = getLogger()
    with Platform('CALCULON') as platform:
        analyzers = [VectorGeneticsAnalyzer()]
        # Set the intervention (Step 2) experiment id you want to analyze
        experiment_id = '<REPLACE_WITH_INTERVENTION_EXPERIMENT_ID>'
        manager = AnalyzeManager(partial_analyze_ok=True, ids=[(experiment_id, ItemType.EXPERIMENT)],
                                 analyzers=analyzers)
        manager.analyze()


