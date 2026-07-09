
# idmtools
from idmtools.builders import SweepArm, ArmType, ArmSimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from itertools import product
# emodpy
from emodpy.emod_task import EMODTask
from emodpy_malaria.reporters.builtin import ReportVectorGenetics
from emodpy_malaria.reporters.builtin import MalariaSummaryReport
from emodpy_malaria.reporters.builtin import ReportVectorStats
from emodpy_malaria.reporters.builtin import ReportMalariaFiltered
from emodpy_malaria.reporters.builtin import add_report_vector_stats, add_report_vector_genetics, \
    add_malaria_summary_report
from helpers import *
import params as params
import manifest as manifest
from functools import partial
import pandas as pd


def get_serialization_paths(platform, serialization_exp_id):
    exp = Experiment.from_id(serialization_exp_id, children=False)
    exp.simulations = platform.get_children(exp.id, exp.item_type,
                                            children=["tags", "configuration", "files", "hpc_jobs"])
    sim_dict = {'Larval_Capacity_arabiensis': [], 'Larval_Capacity_funestus': [], 'Outpath': []}
    for simulation in exp.simulations:
        # if simulation.tags['Run_Number'] == 0:
        string = simulation.get_platform_object().hpc_jobs[0].working_directory.replace('internal.idm.ctr', 'mnt')
        string = string.replace('\\', '/')
        string = string.replace('IDM2', 'idm2')
        # sim_dict['Larval_Capacity'] += [float(simulation.tags['Larval_Capacity'])]
        # sim_dict['Outpath'] += [string]
        sim_dict['Larval_Capacity_arabiensis'] += [float(simulation.tags['Larval_Capacity_arabiensis'])]
        sim_dict['Larval_Capacity_funestus'] += [float(simulation.tags['Larval_Capacity_funestus'])]
        sim_dict['Outpath'] += [string]

    df = pd.DataFrame(sim_dict)
    return df


def general_sim(serialization=0, serialized_exp_id=None):
    """
    This function is designed to be a parameterized version of the sequence of things we do
    every time we run an emod experiment.
    """

    # Create a platform
    # Show how to dynamically set priority and node_group
    platform = Platform("SLURM")

    # create EMODTask
    print("Creating EMODTask (from files)...")

    task = EMODTask.from_default2(
        config_path="my_config.json",
        eradication_path=manifest.eradication_path,
        ep4_custom_cb=None,
        campaign_builder=None,
        schema_path=manifest.schema_file,
        param_custom_cb=set_param_fn,
        demog_builder=build_demog,
    )
    task.set_sif(manifest.sif_path)

    # Create simulation sweep with builder
    builder = ArmSimulationBuilder()
    # Add asset
    task.common_assets.add_asset("../download/schema.json")
    if serialized_exp_id:
        sim_duration = 9 * 365  # 9-year intervention period (three 3-year deployment blocks)
        serialized_population_path_df = get_serialization_paths(platform=platform,
                                                                serialization_exp_id=serialized_exp_id)
        # define our Sweep Arm
        arm = SweepArm(type=ArmType.cross)
        funestus = [8.9]
        arabiensis = [8.8]
        lh_prod = list(product(funestus, arabiensis))
        species = ['funestus', 'arabiensis']
        func = partial(update_serialize, species=species, serialization=serialization, sim_duration=sim_duration,
                       serialized_population_path_df=serialized_population_path_df)
        arm.add_sweep_definition(func, lh_prod)

        # Stochastic realizations
        arm.add_sweep_definition(sweep_sim_random_seed, range(params.nSims))
        baseline = [False]
        start_times = [60]      # nets first deployed ~early March (start of peak season)
        coverage = [0.8]
        start_times_IRS = [0]
        coverage_IRS = [0]
        scenarios = list(SCENARIOS.keys())
        sim_params = list(product(baseline, start_times, coverage,
                                  start_times_IRS, coverage_IRS, scenarios))
        func = partial(update_camp_type, serialize=serialization, sim_duration=sim_duration)
        arm.add_sweep_definition(func, sim_params)
        builder.add_arm(arm)

        reporter = MalariaSummaryReport()  # Create the reporter
        reporter.config(msr_config_builder, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)  # Add the reporter
        # add_malaria_summary_report(task, manifest,)


        reporter = ReportMalariaFiltered()  # Create the reporter
        reporter.config(rmf_config_builder, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)
        ###Another reporter
        reporter = ReportVectorGenetics()  # Create the reporter
        reporter.config(rvg_genome_config_builder1, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)

        reporter = ReportVectorGenetics()  # Create the reporter
        reporter.config(rvg_genome_config_builder2, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)

        reporter = ReportVectorStats()  # Create the reporter
        reporter.config(rvs_config_builder, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)  # Add the reporter

        exp_name = 'ITN_sweep_resistance_test'

    else:
        arm = SweepArm(type=ArmType.cross)
        funestus = [8.9]
        arabiensis = [8.8]
        lh_prod = list(product(funestus, arabiensis))
        species = ['funestus', 'arabiensis']
        func = partial(update_serialize, species=species, serialization=serialization, sim_duration=40 * 365,
                       serialized_population_path_df=None)
        arm.add_sweep_definition(func, lh_prod)

        arm.add_sweep_definition(sweep_sim_random_seed, [0])

        baseline = [True]
        start_times = [0]
        coverage = [0]
        start_times_IRS = [0]
        coverage_IRS = [0]
        sim_params = list(product(baseline, start_times, coverage, start_times_IRS, coverage_IRS))
        func = partial(update_camp_type, serialize=serialization, sim_duration=40 * 365)
        arm.add_sweep_definition(func, sim_params)

        reporter = ReportMalariaFiltered()  # Create the reporter
        reporter.config(rmf_config_builder, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)  # Add the reporter

        reporter = MalariaSummaryReport()  # Create the reporter
        reporter.config(msr_config_builder, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)  # Add the reporter

        reporter = ReportVectorGenetics()  # Create the reporter
        reporter.config(rvg_genome_config_builder1, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)

        reporter = ReportVectorGenetics()  # Create the reporter
        reporter.config(rvg_genome_config_builder2, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)

        reporter = ReportVectorStats()  # Create the reporter
        reporter.config(rvs_config_builder, manifest)  # Config the reporter
        task.reporters.add_reporter(reporter)  # Add the reporter


        exp_name = params.exp_name + '_serialization'
        builder.add_arm(arm)
    # create experiment from builder
    print(f"Prompting for COMPS creds if necessary...")
    experiment = Experiment.from_builder(builder, task, name=exp_name)
    # The last step is to call run() on the ExperimentManager to run the simulations.
    experiment.run(wait_until_done=True, platform=platform)

    # Check result
    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.\n")
        exit()

    print(f"Experiment {experiment.uid} succeeded.")

    # Save experiment id to file
    with open("COMPS_ID", "w") as fd:
        fd.write(experiment.uid.hex)
    print()
    print(experiment.uid.hex)


if __name__ == "__main__":
    # Step 1 (burn-in): un-comment the next two lines and run first
    # serialization = 1
    # serialization_experiment_id = None

    # Step 2 (intervention): set the burn-in experiment id below
    serialization = 0
    serialization_experiment_id = '0dafa2b5-6770-f111-92e9-000d3af5294c'

    general_sim(serialization=serialization, serialized_exp_id=serialization_experiment_id)
