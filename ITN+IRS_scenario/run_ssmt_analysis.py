"""
Server-side (SSMT) analysis of the intervention experiment on COMPS.
Runs SummaryReportAnalyzer, InsetChartAnalyzer, and VectorGeneticsAnalyzer.
"""


import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from emod_analyzers.SummaryReportAnalyzer import SummaryReportAnalyzer
from emod_analyzers.InsetChartAnalyzer import InsetChartAnalyzer
from emod_analyzers.VectorGeneticsAnalyzer import VectorGeneticsAnalyzer

from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis


if __name__ == "__main__":
    platform = Platform('CALCULON')
    analysis = PlatformAnalysis(
        platform=platform,
        experiment_ids=["8aa20273-596d-f111-92e9-000d3af5294c"],
        analyzers=[SummaryReportAnalyzer, InsetChartAnalyzer, VectorGeneticsAnalyzer],
        analyzers_args=[
            # SummaryReportAnalyzer: tags that match what update_camp_type returns
            {"tags": ["Scenario", "Coverage",
                      "Larval_Capacity_funestus", "Larval_Capacity_arabiensis"]},
            # InsetChartAnalyzer: same tags
            {"tags": ["Scenario", "Coverage",
                      "Larval_Capacity_funestus", "Larval_Capacity_arabiensis"]},
            # VectorGeneticsAnalyzer: two loci (a = PBO, b = IRS)
            {"species": ("funestus", "arabiensis"), "loci": ("a", "b"),
             "expt_name": "ITN_IRS",
             "tags": ["Scenario", "Coverage",
                      "Larval_Capacity_funestus", "Larval_Capacity_arabiensis"]},
        ],
        analysis_name="PBO-ITN + IRS scenario",
        extra_args=dict(partial_analyze_ok=True),
    )
    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)
