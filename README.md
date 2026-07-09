# Modeling the adoption of new insecticide-treated nets (ITNs)

EMOD (`emodpy-malaria`) simulation and analysis code for a study of next-generation
ITN strategies and their effect on malaria transmission and insecticide-resistance
dynamics in two vector species, *Anopheles funestus* and *An. arabiensis*.

The model represents insecticide resistance with independent genetic loci, one per net
chemistry, so that deploying one chemistry selects for resistance only to that chemistry.
Three self-contained scenarios reproduce the main results:

- **ITN comparison** - standard pyrethroid vs pyrethroid-PBO vs Interceptor G2 (chlorfenapyr) nets.
- **ITN transition** - four multi-year ITN deployment strategies (e.g. STD -> PBO -> IG2).
- **ITN + IRS** - pyrethroid-PBO nets alone vs nets combined with indoor residual spraying (IRS).

## Requirements

The scripts run on the IDM modeling stack. `emodpy-malaria` pulls in `idmtools`,
`emodpy`, and `emod-api`, and ships the compiled EMOD (`Eradication`) binary and schema:

```bash
pip install -r requirements.txt --index-url https://packages.idmod.org/api/pypi/pypi-production/simple
```

`numpy`, `pandas`, and `matplotlib` are available from PyPI. See `requirements.txt` for details.

## Model binary and schema

The `Eradication` binary and `schema.json` are not committed to this repository; they are
provided by `emodpy-malaria`. Stage them into a local `download/` directory before running:

```python
import emod_malaria.bootstrap as dtk
dtk.setup("download")
```

Each scenario's `manifest.py` expects `download/Eradication` and `download/schema.json`
one level above the scenario directory.

## Infrastructure

Simulations are submitted to an HPC cluster through idmtools (the COMPS / SLURM platform;
`idmtools.ini` in each scenario folder holds the platform settings). Server-side analysis
runs as an SSMT work item on COMPS. Running elsewhere requires adjusting the platform
configuration.

## Directory organization

- `input_files/` - shared demographics input (`single_node_demographics.json`) used by every scenario.
- `emod_analyzers/` - idmtools analyzers shared by all scenarios:
  - `InsetChartAnalyzer.py` - prevalence, adult vectors, EIR, and clinical incidence over time.
  - `SummaryReportAnalyzer.py` - monthly clinical incidence by age bin.
  - `VectorGeneticsAnalyzer.py` - per-locus genotype fractions (susceptible / heterozygous / resistant) from the GENOME report.
- `ITN_comparison_scenario/` - standard vs PBO vs IG2 nets (single 3-year deployment block).
- `ITN_transition_scenario/` - four ITN deployment / rotation strategies over 9 years (three 3-year blocks).
- `ITN+IRS_scenario/` - pyrethroid-PBO nets with and without IRS, comparing IRS timing in year 2 vs year 3.
- `download/` - (created locally, not committed) the EMOD binary and schema.

Each scenario directory contains the same set of files:

- `manifest.py` - paths to the EMOD binary, schema, and Singularity image (`dtk_sif.id`).
- `params.py` - experiment name and number of stochastic realizations (`nSims`).
- `helpers.py` - builders for the config, demographics, campaign, resistance loci, and reporters.
- `*_run.py` - runs the simulations (burn-in and intervention phases).
- `run_ssmt_analysis.py` - runs the analyzers on the intervention experiment (SSMT on COMPS).
- `plot.py` / `plot_figure.py` - produces the figure from the analyzer CSV output.
- `idmtools.ini`, `dtk_sif.id` - platform and Singularity-image settings.

## Workflow

Each scenario follows the same four steps. Run them from inside the scenario directory.

1. **Burn-in.** In the `*_run.py` script set `serialization = 1` and
   `serialization_experiment_id = None`, then run it. This simulates ~40 years to reach a
   stable, age-structured population and serializes the final state. Note the experiment id
   that is printed.
2. **Intervention.** Set `serialization = 0` and `serialization_experiment_id` to the burn-in
   experiment id from step 1, then run the script again to launch the intervention sweep
   (`nSims` stochastic realizations per scenario).
3. **Analysis.** Set the intervention experiment id in `run_ssmt_analysis.py` and run it to
   produce the incidence, prevalence, and genotype-fraction CSVs on COMPS.
4. **Plot.** Download the analysis CSVs into the scenario directory and run `plot.py`
   (or `plot_figure.py` for the ITN+IRS scenario) to generate the figures in `output/`.

## Notes

- The plotting scripts default to matplotlib's interactive `MacOSX` backend. On Linux or
  for headless runs, switch to the `Agg` backend at the top of the script.
- Experiment ids in the scripts are examples from our runs; replace them with the ids
  produced by your own burn-in and intervention experiments.
