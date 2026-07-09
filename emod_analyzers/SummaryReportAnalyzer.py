
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
from idmtools.analysis.csv_analyzer import CSVAnalyzer
from idmtools.core import ItemType
from idmtools.core.platform_factory import Platform

mpl.use('Agg')

class SummaryReportAnalyzer(BaseAnalyzer):

    def __init__(self, title='idm', tags=None):
        super().__init__(filenames=["output\\MalariaSummaryReport_Monthly.json"])
        self.tags = tags if tags is not None else [
            'Scenario', 'Baseline', 'Start_Day', 'Coverage', 'Resistance', 'Start_Day_IRS', 'Coverage_IRS',
            'Larval_Capacity_arabiensis', 'Larval_Capacity_funestus',
        ]
        self.age_bins = ['5', '15', '125']
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
        incidence = datatemp['DataByTimeAndAgeBins']['Annual Clinical Incidence by Age Bin']
        true_incidence = np.array(incidence)
        # population = datatemp['DataByTimeAndAgeBins']['Average Population by Age Bin']
        # true_incidence = np.array(incidence) * np.array(population)
        df = pd.DataFrame()
        for i in range(np.shape(true_incidence)[0]):
            incidence_dict = dict(zip(self.age_bins,
                                      list(true_incidence[i])))
            df = pd.concat([df, pd.DataFrame(incidence_dict, index=[i])])
        df.reset_index(names=['Day'], inplace=True)

        return df

    def reduce(self, all_data: Dict[Union[IWorkflowItem, Simulation], Any]) -> Any:
        """
        Create the Final Population JSON and Plot
        Args:
            all_data: Populate data from all the Simulations
        Returns:
            None
        """
        output_dir = os.path.join(self.working_dir, "output")
        df = pd.DataFrame()
        for s, v in all_data.items():
            dftemp = v.copy()
            for t in self.tags:
                dftemp[t] = [s.tags[t]] * len(v)
            df = pd.concat([df, dftemp])
        df.to_csv(os.path.join(output_dir, "incidence_all.csv"))

        groupby_tags = self.tags
        # groupby_tags.remove('Run_Number')
        df_final = df.groupby(groupby_tags + ['Day'])[self.age_bins].mean().reset_index()
        df_final_std = df.groupby(groupby_tags + ['Day'])[self.age_bins].std()
        for c in self.age_bins:
            df_final[c + '_std'] = list(df_final_std[c])
        df_final.to_csv(os.path.join(output_dir, "incidence.csv"))


if __name__ == '__main__':
    # Set the platform where you want to run your analysis
    # In this case we are running in BELEGOST since the Work Item we are analyzing was run on COMPS
    logger = getLogger()
    with Platform('CALCULON') as platform:
        # Initialize the analyser class with the path of the output csv file
        analyzers = [SummaryReportAnalyzer()]

        # Set the experiment id you want to analyze
        experiment_id = 'f36d11e6-8cdf-eb11-a9ec-b88303911bc1'

        # Specify the id Type, in this case an Experiment on COMPS
        manager = AnalyzeManager(partial_analyze_ok=True, ids=[(experiment_id, ItemType.EXPERIMENT)],
                                 analyzers=analyzers)
        manager.analyze()



