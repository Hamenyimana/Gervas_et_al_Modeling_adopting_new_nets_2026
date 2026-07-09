import json
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


# InsetChart.json

class InsetChartAnalyzer(BaseAnalyzer):

    def __init__(self, title='idm', tags=None):
        super().__init__(filenames=["output\\ReportMalariaFiltered_Yearly.json"])
        self.tags = tags if tags is not None else [
            'Scenario', 'Baseline', 'Start_Day', 'Coverage', 'Resistance', 'Start_Day_IRS', 'Coverage_IRS',
            'Larval_Capacity_arabiensis', 'Larval_Capacity_funestus',
        ]
        self.age_bins = ['5']
        print(title)

    def initialize(self):
        """
        Initialize our Analyzer. At the moment, this just creates our output folder
        Returns:
        """
        if not os.path.exists(os.path.join(self.working_dir, "output")):
            os.mkdir(os.path.join(self.working_dir, "output"))

    def map(self, data: Dict[str, Any], item: Union[IWorkflowItem, Simulation]) -> Any:
        """
        Extracts the Statistical Population, Data channel from InsetChart.
        Called for Each WorkItem/Simulation.
        Args:
            data: Data mapping str to content of file
            item: Item to Extract Data from(Usually a Simulation)
        Returns:
        """
        datatemp = data[self.filenames[0]]

        prevalence = datatemp['Channels']['PfHRP2 Prevalence']['Data']
        true_prevalence = datatemp['Channels']['True Prevalence']['Data']
        adult_vector = datatemp['Channels']['Adult Vectors']['Data']
        incidence = datatemp['Channels']['New Clinical Cases']['Data']
        annual_eir = sum(datatemp['Channels']['Daily EIR']['Data']) / (
                    len((datatemp['Channels']['Daily EIR']['Data'])) / 365)
        elimination = 0
        if sum(true_prevalence[4 * 365:5 * 365]) == 0:
            elimination = 1

        df_p = pd.DataFrame(list(zip([i for i in range(len(prevalence))], prevalence, true_prevalence, adult_vector)),
                            columns=['Time', 'RDT Prevalence', 'True Prevalence', 'Adult Vectors'])
        df_i = pd.DataFrame(list(zip([sum(incidence)], [elimination], [annual_eir])),
                            columns=['New Clinical Cases', 'Elimination', 'Annual EIR'])

        return df_p, df_i

    def reduce(self, all_data: Dict[Union[IWorkflowItem, Simulation], Any]) -> Any:
        """
        Create the Final Population JSON and Plot
        Args:
            all_data: Populate data from all the Simulations
        Returns:
            None
        """
        output_dir = os.path.join(self.working_dir, "output")
        df_prevalence = pd.DataFrame()
        df_incidence = pd.DataFrame()
        for s, v in all_data.items():
            dftemp_prev = v[0].copy()
            dftemp_incidence = v[1].copy()
            for t in self.tags:
                dftemp_prev[t] = [s.tags[t]] * len(v[0])
                dftemp_incidence[t] = [s.tags[t]] * len(v[1])
            df_prevalence = pd.concat([df_prevalence, dftemp_prev])
            df_incidence = pd.concat([df_incidence, dftemp_incidence])
        df_prevalence.to_csv(os.path.join(output_dir, "ITN_distribution_example_Prevalence_full.csv"))
        df_incidence.to_csv(os.path.join(output_dir, "ITN_distribution_example_incidence_full.csv"))

        groupby_tags = self.tags
        # groupby_tags.remove('Run_Number')
        df_prev_final = df_prevalence.groupby(groupby_tags + ['Time'])[['RDT Prevalence', 'True Prevalence',
                                                                        'Adult Vectors']].mean().reset_index()
        df_prev_final_std = df_prevalence.groupby(groupby_tags + ['Time'])[['RDT Prevalence', 'True Prevalence',
                                                                            'Adult Vectors']].std()
        for c in ['RDT Prevalence', 'True Prevalence', 'Adult Vectors']:
            df_prev_final[c + '_std'] = list(df_prev_final_std[c])
        df_inc_final = df_incidence.groupby(groupby_tags)[['New Clinical Cases', 'Elimination',
                                                           'Annual EIR']].mean().reset_index()
        df_inc_final['New Clinical Cases_std'] = list(df_incidence.groupby(groupby_tags)['New Clinical Cases'].std())
        df_inc_final['Annual EIR_std'] = list(
            df_incidence.groupby(groupby_tags)['Annual EIR'].std())

        df_prev_final.to_csv(os.path.join(output_dir, "ITN_distribution_example_Prevalence.csv"))
        df_inc_final.to_csv(os.path.join(output_dir, "ITN_distribution_example_incidence.csv"))


if __name__ == '__main__':
    # Set the platform where you want to run your analysis
    # In this case we are running in BELEGOST since the Work Item we are analyzing was run on COMPS
    logger = getLogger()
    with Platform('CALCULON') as platform:
        # Initialize the analyser class with the path of the output csv file
        analyzers = [InsetChartAnalyzer()]

        # Set the experiment id you want to analyze
        experiment_id = '7cb0e02a-6483-ee11-aa0c-b88303911bc1'

        # Specify the id Type, in this case an Experiment on COMPS
        manager = AnalyzeManager(partial_analyze_ok=True, ids=[(experiment_id, ItemType.EXPERIMENT)],
                                 analyzers=analyzers)
        manager.analyze()
