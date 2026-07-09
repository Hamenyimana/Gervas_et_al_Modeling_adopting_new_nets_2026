import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from emod_analyzers.SummaryReportAnalyzer import SummaryReportAnalyzer
from emod_analyzers.InsetChartAnalyzer import InsetChartAnalyzer
from emod_analyzers.VectorGeneticsAnalyzer import VectorGeneticsAnalyzer
from idmtools.core.platform_factory import Platform
from idmtools.analysis.platform_anaylsis import PlatformAnalysis
if __name__ == "__main__":
    platform = Platform('CALCULON')
    # Analyze the Step 2 (intervention) experiment. VectorGeneticsAnalyzer produces the
    # per-locus genotype fractions; InsetChart/SummaryReport produce incidence/prevalence.
    # Add/remove analyzers from the list as needed.
    analysis = PlatformAnalysis(platform=platform,
                                experiment_ids=["94630b92-6774-f111-92e9-000d3af5294c"],
                                analyzers=[VectorGeneticsAnalyzer, InsetChartAnalyzer, SummaryReportAnalyzer],
                                analyzers_args=[{}, {}, {}],
                                analysis_name="Gervas_et_al_Modeling_adopting_new_ITNs_2026",
                                extra_args=dict(partial_analyze_ok=True)
                                )

    analysis.analyze(check_status=True)
    wi = analysis.get_work_item()
    print(wi)

