
from idmtools.builders import SweepArm, ArmType, ArmSimulationBuilder
from idmtools.core.platform_factory import Platform
from idmtools.entities.experiment import Experiment
from itertools import product
from functools import partial

from emodpy.emod_task import EMODTask
from emodpy_malaria.reporters.builtin import (ReportVectorGenetics, MalariaSummaryReport,
                                              ReportVectorStats, ReportMalariaFiltered)

from helpers import (set_param_fn, build_demog, update_serialize, update_camp_type,
                     sweep_sim_random_seed, msr_config_builder, rmf_config_builder,
                     rvg_genome_funestus, rvg_genome_arabiensis, rvs_config_builder,
                     SCENARIOS)
import params
import manifest
import pandas as pd


def get_serialization_paths(platform, serialization_exp_id):
    exp = Experiment.from_id(serialization_exp_id, children=False)
    exp.simulations = platform.get_children(exp.id, exp.item_type,
                                            children=["tags", "configuration", "files", "hpc_jobs"])
    rows = {'Larval_Capacity_arabiensis': [], 'Larval_Capacity_funestus': [], 'Outpath': []}
    for sim in exp.simulations:
        wd = sim.get_platform_object().hpc_jobs[0].working_directory
        wd = wd.replace('internal.idm.ctr', 'mnt').replace('\\', '/').replace('IDM2', 'idm2')
        rows['Larval_Capacity_funestus'].append(float(sim.tags['Larval_Capacity_funestus']))
        rows['Larval_Capacity_arabiensis'].append(float(sim.tags['Larval_Capacity_arabiensis']))
        rows['Outpath'].append(wd)
    return pd.DataFrame(rows)


def general_sim(serialization=0, serialized_exp_id=None):
    platform = Platform("SLURM")

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
    task.common_assets.add_asset("../download/schema.json")

    builder = ArmSimulationBuilder()
    species = ['funestus', 'arabiensis']

    if serialized_exp_id:
        # intervention runs
        sim_duration = 3 * 365
        serialized_df = get_serialization_paths(platform, serialized_exp_id)

        arm = SweepArm(type=ArmType.cross)
        lh_prod = list(product([8.9], [8.8]))
        func = partial(update_serialize, species=species, serialization=serialization,
                       sim_duration=sim_duration, serialized_population_path_df=serialized_df)
        arm.add_sweep_definition(func, lh_prod)

        arm.add_sweep_definition(sweep_sim_random_seed, range(params.nSims))

        func = partial(update_camp_type, coverage=0.8, serialize=serialization,
                       sim_duration=sim_duration)
        arm.add_sweep_definition(func, list(SCENARIOS.keys()))
        builder.add_arm(arm)

        # reporters
        for reporter_fn in (rvg_genome_arabiensis, rvg_genome_funestus):
            r = ReportVectorGenetics()
            r.config(reporter_fn, manifest)
            task.reporters.add_reporter(r)

        r = MalariaSummaryReport()
        r.config(msr_config_builder, manifest)
        task.reporters.add_reporter(r)

        r = ReportMalariaFiltered()
        r.config(rmf_config_builder, manifest)
        task.reporters.add_reporter(r)

        r = ReportVectorStats()
        r.config(rvs_config_builder, manifest)
        task.reporters.add_reporter(r)

        exp_name = params.exp_name + '_intervention'

    else:
        # burn-in
        sim_duration = 40 * 365
        arm = SweepArm(type=ArmType.cross)
        lh_prod = list(product([8.9], [8.8]))
        func = partial(update_serialize, species=species, serialization=serialization,
                       sim_duration=sim_duration, serialized_population_path_df=None)
        arm.add_sweep_definition(func, lh_prod)
        arm.add_sweep_definition(sweep_sim_random_seed, [0])

        func = partial(update_camp_type, coverage=0, serialize=serialization,
                       sim_duration=sim_duration)
        arm.add_sweep_definition(func, ['PBO_only'])  # dummy scenario; no nets during burn-in
        builder.add_arm(arm)

        exp_name = params.exp_name + '_burnin'

    experiment = Experiment.from_builder(builder, task, name=exp_name)
    experiment.run(wait_until_done=True, platform=platform)

    if not experiment.succeeded:
        print(f"Experiment {experiment.uid} failed.")
        exit()

    print(f"Experiment {experiment.uid} succeeded.")
    with open("COMPS_ID", "w") as fd:
        fd.write(experiment.uid.hex)
    print(experiment.uid.hex)


if __name__ == "__main__":

    # Step 1 (burn-in): un-comment the next two lines and run first
    # serialization = 1
    # serialization_experiment_id = None

    # Step 2 (intervention): set the burn-in experiment id below
    serialization = 0
    serialization_experiment_id = '1c9c7e16-4375-f111-92e9-000d3af5294c'
    general_sim(serialization=serialization, serialized_exp_id=serialization_experiment_id)
