"""
Server-side (SSMT) analysis of the intervention experiment on COMPS.
Runs SummaryReportAnalyzer (monthly incidence/prevalence) and InsetChartAnalyzer.
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
        experiment_ids=["186197ec-9874-f111-92e9-000d3af5294c"],
        analyzers=[SummaryReportAnalyzer, InsetChartAnalyzer],
        analyzers_args=[
            # SummaryReportAnalyzer: tags that match what update_camp_type returns
            {"tags": ["Scenario", "Coverage",
                      "Larval_Capacity_funestus", "Larval_Capacity_arabiensis"]},
            # InsetChartAnalyzer: same tags
            {"tags": ["Scenario", "Coverage",
                      "Larval_Capacity_funestus", "Larval_Capacity_arabiensis"]},
            # # VectorGeneticsAnalyzer (uncomment to add per-locus genotype fractions):
            # {"species": ("funestus", "arabiensis"), "loci": ("a", "b", "c"),
            #  "expt_name": "ITN_comparison",
            #  "tags": ["Scenario", "Coverage",
            #           "Larval_Capacity_funestus", "Larval_Capacity_arabiensis"]},
        ],
        analysis_name="STD vs PBO vs IG2 ITN comparison",
        extra_args=dict(partial_analyze_ok=True),
    )
    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)
