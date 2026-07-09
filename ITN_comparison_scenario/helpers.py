import os
import numpy as np
from functools import partial

from emodpy_malaria.interventions.irs import add_scheduled_irs_housing_modification
import manifest

import emod_api.demographics.Demographics as Demographics
from emodpy_malaria import malaria_config as malconf
from emodpy_malaria import vector_config as vecconf
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
from emodpy_malaria.interventions.usage_dependent_bednet import add_scheduled_usage_dependent_bednet
from emodpy_malaria.vector_config import add_insecticide_resistance
import emod_api.campaign as camp


STD_RESISTANCE = {
    'funestus': {
        'heterozygous': 0.5,
        'homozygous': 0.25
    },
    'arabiensis': {
        'heterozygous': 0.51,
        'homozygous': 0.255
    }
}

PBO_RESISTANCE = {
    'funestus': {
        'heterozygous': 0.67,
        'homozygous': 0.335
    },
    'arabiensis': {
        'heterozygous': 0.715,
        'homozygous': 0.3575
    }
}

IG2_RESISTANCE = {
    'funestus': {
        'heterozygous': 0.73,
        'homozygous': 0.365
    },
    'arabiensis': {
        'heterozygous': 0.8,
        'homozygous': 0.4
    }
}
# One net selects one locus
NET_LOCI = [('STD', 'b'), ('PBO-bednet', 'b'), ('IG2', 'c')]

# Initial resistant-allele frequencies per locus
INITIAL_RESISTANCE_FREQ = {
    'funestus': {'a': 0.012, 'b': 0.01, 'c': 0.005},
    'arabiensis': {'a': 0.011, 'b': 0.01, 'c': 0.009},
}


STD_SPEC = dict(
    killing_initial_effect=0.93, killing_decay_time_constant=1420,
    blocking_initial_effect=0.9, blocking_decay_time_constant=730,
)
PBO_SPEC = dict(
    killing_initial_effect=0.98, killing_decay_time_constant=1260,
    blocking_initial_effect=0.92, blocking_decay_time_constant=730,
)

IG2_SPEC = dict(
    killing_initial_effect=0.987, killing_decay_time_constant=1580,
    blocking_initial_effect=0.932, blocking_decay_time_constant=730,
)

ITN_START_DAY = 60

N_BLOCKS = 1  # single 3-year deployment block
BLOCK_DAYS = 3 * 365

# ITN types to compare
SCENARIOS = {
    'STD': None,
    'PBO': None,
    'IG2': None,
}

def sweep_sim_random_seed(simulation, value):
    simulation.task.config.parameters.Run_Number = value
    return {"Run_Number": value}


def sweep_sim_larval_capacity(simulation, species, value):
    for j, s in enumerate(species):
        for i, vsp in enumerate(simulation.task.config.parameters.Vector_Species_Params):
            if vsp['Name'] == s:
                simulation.task.config.parameters.Vector_Species_Params[i].Habitats[0]['Max_Larval_Capacity'] \
                    = pow(10, value[j])
    return {"Larval_Capacity": pow(10, value[j])}


def update_serialize(simulation, larval_multiplier, species, serialization=0,
                     sim_duration=40 * 365, serialized_population_path_df=None):
    if serialization:
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Time_Steps = [sim_duration]
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'TIMESTEP'
        sweep_sim_larval_capacity(simulation, species, value=larval_multiplier)
    else:
        serialized_population_path = serialized_population_path_df[
            (serialized_population_path_df['Larval_Capacity_' + species[0]] == larval_multiplier[0])
            & (serialized_population_path_df['Larval_Capacity_' + species[1]] == larval_multiplier[1])
            ]['Outpath'].values[0]
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Mask_Node_Read = 0
        simulation.task.config.parameters.Serialization_Mask_Node_Write = 0
        simulation.task.config.parameters.Serialized_Population_Path = os.path.join(
            serialized_population_path, 'output')
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'READ'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Filenames = ['state-14600.dtk']
        sweep_sim_larval_capacity(simulation, species, value=larval_multiplier)

    lh_dict = dict(zip([f'Larval_Capacity_{s}' for s in species], larval_multiplier))
    return {"Serialization": serialization, **lh_dict}


def update_camp_type(simulation, scenario, coverage=0.8, serialize=0, sim_duration=9 * 365):
    """Sweep callback: build the campaign for the given scenario."""
    build_fn = partial(build_camp, scenario=scenario, coverage=coverage,
                       serialize=serialize, sim_duration=sim_duration)
    simulation.task.create_campaign_from_callback(build_fn)
    return {"Scenario": scenario, "Coverage": coverage}


def build_camp(scenario, coverage=0.8, serialize=0, sim_duration=9 * 365):
    """Build the campaign for a scenario.

    Burn-in (serialize=True): treatment-seeking only.
    Intervention: treatment-seeking plus the selected ITN type.
    """
    camp.set_schema(manifest.schema_file)

    if serialize:
        add_treatment_seeking(camp,
                              targets=[{"trigger": "NewClinicalCase", "coverage": 0.4,
                                        "agemin": 15, "agemax": 70, "rate": 0.3},
                                       {"trigger": "NewSevereCase", "coverage": 0.8, "rate": 0.5}],
                              drug=['Artemether', 'Lumefantrine'],
                              start_day=sim_duration - 10 * 365,
                              broadcast_event_name='Received_Treatment')
        return camp

    # treatment seeking (intervention phase)
    add_treatment_seeking(camp,
                          targets=[{"trigger": "NewClinicalCase", "coverage": 0.4,
                                    "agemin": 15, "agemax": 70, "rate": 0.3},
                                   {"trigger": "NewSevereCase", "coverage": 0.8, "rate": 0.5}],
                          drug=['Artemether', 'Lumefantrine'],
                          start_day=0,
                          broadcast_event_name='Received_Treatment')



    # deploy the selected ITN type
    NET_SPECS = {
        "STD": ("STD", STD_SPEC),
        "PBO": ("PBO-bednet", PBO_SPEC),
        "IG2": ("IG2", IG2_SPEC),
    }

    if scenario not in NET_SPECS:
        raise ValueError(f"Unknown scenario: {scenario}")

    insecticide, spec = NET_SPECS[scenario]

    for block in range(N_BLOCKS):
        itn_day = ITN_START_DAY + block * BLOCK_DAYS

        add_scheduled_usage_dependent_bednet(
            camp,
            start_day=itn_day,
            demographic_coverage=coverage,
            insecticide=insecticide,
            discard_config={
                'Expiration_Period_Distribution': "EXPONENTIAL_DISTRIBUTION",
                'Expiration_Period_Exponential': 730
            },
            killing_initial_effect=spec['killing_initial_effect'],
            killing_decay_time_constant=spec['killing_decay_time_constant'],
            blocking_initial_effect=spec['blocking_initial_effect'],
            blocking_decay_time_constant=spec['blocking_decay_time_constant'],
        )
        return camp
def build_demog():
    demog = Demographics.from_file("../input_files/single_node_demographics.json")
    return demog


def set_param_fn(config):
    config = malconf.set_team_defaults(config, manifest)

    config.parameters.Base_Rainfall = 150
    config.parameters.Climate_Model = "CLIMATE_CONSTANT"
    config.parameters.Enable_Disease_Mortality = 0
    config.parameters.Enable_Vector_Species_Report = 0
    config.parameters.Enable_Initial_Prevalence = 1
    config.parameters.Enable_Vector_Aging = 1
    config.parameters.Enable_Natural_Mortality = 1
    config.parameters.Demographics_Filenames = ['single_node_demographics.json']

    # An. funestus
    vecconf.add_species(config, manifest, species_to_select=["funestus"])
    vecconf.set_species_param(config, species='funestus', parameter="Habitats",
                              value=[vecconf.configure_linear_spline(manifest, max_larval_capacity=pow(10, 6),
                                                                     capacity_distribution_number_of_years=1,
                                                                     capacity_distribution_over_time={
                                                                         "Times": [0.0, 30.417, 60.833, 91.25, 121.667,
                                                                                   152.083, 182.5, 212.917, 243.333,
                                                                                   273.75, 304.167, 334.583],
                                                                         "Values": [0.16, .11, .25, .3, .35, 0.9, .3,
                                                                                    .2, .25, .2, .1, .15]})],
                              overwrite=True)
    for locus, freq in INITIAL_RESISTANCE_FREQ['funestus'].items():
        vecconf.add_genes_and_alleles(config, manifest, 'funestus',
                                      [(f'{locus}0', 1 - freq), (f'{locus}1', freq)])
    vecconf.set_species_param(config, species='funestus', parameter="Indoor_Feeding_Fraction",
                              value=0.782, overwrite=True)
    vecconf.set_species_param(config, species='funestus', parameter="Anthropophily",
                              value=0.64, overwrite=True)
    vecconf.set_species_param(config, species='funestus', parameter="Adult_Life_Expectancy",
                              value=25, overwrite=True)

    # An. arabiensis
    vecconf.add_species(config, manifest, species_to_select=["arabiensis"])
    vecconf.set_species_param(config, species='arabiensis', parameter="Habitats",
                              value=[vecconf.configure_linear_spline(manifest, max_larval_capacity=pow(10, 6),
                                                                     capacity_distribution_number_of_years=1,
                                                                     capacity_distribution_over_time={
                                                                         "Times": [0.0, 30.417, 60.833, 91.25, 121.667,
                                                                                   152.083, 182.5, 212.917, 243.333,
                                                                                   273.75, 304.167, 334.583],
                                                                         "Values": [.19, .22, .25, .3, .8, .18, .12,
                                                                                    .05, .042, .02, .11, .09]})],
                              overwrite=True)
    for locus, freq in INITIAL_RESISTANCE_FREQ['arabiensis'].items():
        vecconf.add_genes_and_alleles(config, manifest, 'arabiensis',
                                      [(f'{locus}0', 1 - freq), (f'{locus}1', freq)])
    vecconf.set_species_param(config, species='arabiensis', parameter="Indoor_Feeding_Fraction",
                              value=0.631, overwrite=True)
    vecconf.set_species_param(config, species='arabiensis', parameter="Anthropophily",
                              value=0.44, overwrite=True)
    vecconf.set_species_param(config, species='arabiensis', parameter="Adult_Life_Expectancy",
                              value=20, overwrite=True)

    # standard ITN resistance (locus a)
    for sp in ['funestus', 'arabiensis']:
        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='STD',
            species=sp,
            allele_combo=[['a1', 'a1']],
            killing=STD_RESISTANCE[sp]['homozygous']
        )

        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='STD',
            species=sp,
            allele_combo=[['a0', 'a1']],
            killing=STD_RESISTANCE[sp]['heterozygous']
        )

    # PBO resistance (locus b)
    for sp in ['funestus', 'arabiensis']:
        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='PBO-bednet',
            species=sp,
            allele_combo=[['b1', 'b1']],
            killing=PBO_RESISTANCE[sp]['homozygous']
        )

        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='PBO-bednet',
            species=sp,
            allele_combo=[['b0', 'b1']],
            killing=PBO_RESISTANCE[sp]['heterozygous']
        )

    # IG2 resistance (locus c)
    for sp in ['funestus', 'arabiensis']:
        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='IG2',
            species=sp,
            allele_combo=[['c1', 'c1']],
            killing=IG2_RESISTANCE[sp]['homozygous']
        )

        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='IG2',
            species=sp,
            allele_combo=[['c0', 'c1']],
            killing=IG2_RESISTANCE[sp]['heterozygous']
        )

    return config

# reporter config builders
def rvg_genome_funestus(params):
    params.Include_Vector_State_Columns = False
    params.Include_Death_By_State_Columns = False
    params.Species = 'funestus'
    params.Stratify_By = 'GENOME'
    params.Combine_Similar_Genomes = True
    params.Start_Day = 0
    params.End_Day = 9 * 365
    return params


def rvg_genome_arabiensis(params):
    params.Include_Vector_State_Columns = False
    params.Include_Death_By_State_Columns = False
    params.Species = 'arabiensis'
    params.Stratify_By = 'GENOME'
    params.Combine_Similar_Genomes = True
    params.Start_Day = 0
    params.End_Day = 9 * 365
    return params


def msr_config_builder(params):
    params.Age_Bins = [5]
    params.Reporting_Interval = 30  # monthly
    params.Filename_Suffix = 'Monthly'
    params.Max_Number_Reports = int(np.ceil(9 * 365 / 30))
    params.End_Day = 9 * 365
    return params


def rmf_config_builder(params):
    params.Start_Day = 1
    params.End_Day = 3 * 365
    params.Min_Age_Years = 0
    params.Max_Age_Years = 5
    params.Filename_Suffix = 'Yearly'
    return params

def rvs_config_builder(params):
    params.Species_List = ['funestus', 'arabiensis']
    params.Stratify_By_Species = True
    params.Include_Death_By_State_Columns = False
    params.Include_Wolbachia_Columns = False
    params.Include_Gestation_Columns = False
    params.Include_Microsporidia_Columns = False
    return params
