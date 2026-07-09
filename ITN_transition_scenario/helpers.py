
import os
import numpy as np
import pandas as pd
import datetime
from functools import partial

from emodpy_malaria.interventions.irs import add_scheduled_irs_housing_modification

import manifest as manifest

import emod_api.demographics.Demographics as Demographics
from emodpy_malaria import malaria_config as malconf
from emodpy_malaria import vector_config as vecconf
from emodpy_malaria.interventions.treatment_seeking import add_treatment_seeking
from emodpy_malaria.interventions.usage_dependent_bednet import (
    add_scheduled_usage_dependent_bednet,
    _get_seasonal_times_and_values,
    _get_age_times_and_values,
)
from emodpy_malaria.interventions import common as malaria_common
from emodpy_malaria.vector_config import add_insecticide_resistance
from emodpy_malaria.reporters.builtin import add_report_malaria_filtered
from emodpy_malaria.reporters.builtin import add_malaria_summary_report
from emod_api import schema_to_class as s2c
from emod_api.interventions import utils as waning_utils
import emod_api.campaign as camp

# Three-locus insecticide-resistance model. Each net type selects its own
# independent locus, so deploying one chemistry does not select for resistance
# to the others:
#     locus a (a0/a1) -> pyrethroid (standard ITN)
#     locus b (b0/b1) -> PBO-net
#     locus c (c0/c1) -> IG2 (chlorfenapyr)

# Resistance mortality modifiers
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

IRS_RESISTANCE = {
    'heterozygous': 0.83,
    'homozygous': 0.415
}
# (insecticide name, resistance-locus letter): one net selects one locus
NET_LOCI = [('pyrethroid', 'a'), ('PBO-synergist', 'b'), ('chlorfenapyr', 'c')]

INITIAL_RESISTANCE_FREQ = {
    'funestus': {'a': 0.045, 'b': 0.016, 'c': 0.01},
    'arabiensis': {'a': 0.04, 'b': 0.015, 'c': 0.009},
}

# Per-component efficacy parameters (each entry is one insecticide in a net).
# 'blocking' and 'block_decay' are omitted for components with no blocking effect.
_PYRETHROID_STD = dict(name='pyrethroid', killing=0.93, kill_decay=1420,
                       blocking=0.9, block_decay=730)

_PYRETHROID_PBO = dict(name='pyrethroid', killing=0.98, kill_decay=1420,
                       blocking=0.92, block_decay=730)
_PBO_SYNERGIST = dict(name='PBO-synergist', killing=0.71, kill_decay=1260)

_PYRETHROID_IG2 = dict(name='pyrethroid', killing=0.987, kill_decay=1420,
                       blocking=0.932, block_decay=730)
_CHLORFENAPYR = dict(name='chlorfenapyr', killing=0.81, kill_decay=1580)

# Net name to the list of component dicts passed to add_multi_insecticide_bednet()
# (STD uses the single-insecticide wrapper; PBO and IG2 use the multi-insecticide one)
NET_COMPONENTS = {
    'STD': None,                              # handled by add_scheduled_usage_dependent_bednet
    'PBO': [_PYRETHROID_PBO, _PBO_SYNERGIST],
    'IG2': [_PYRETHROID_IG2, _CHLORFENAPYR],
}

# STD uses the single-insecticide path
_STD_SPEC = dict(insecticide='pyrethroid',
                 killing_initial_effect=0.93, killing_decay_time_constant=1420,
                 blocking_initial_effect=0.9, blocking_decay_time_constant=730)

# Four ITN-deployment strategies (one net per 3-year block).
SCENARIOS = {
    'STD-STD-STD': ['STD', 'STD', 'STD'],  # repeated standard pyrethroid-only
    'STD-PBO-STD': ['STD', 'PBO', 'STD'],  # rotation standard/PBO
    'STD-PBO-PBO': ['STD', 'PBO', 'PBO'],  # transition to PBO
    'STD-PBO-IG2': ['STD', 'PBO', 'IG2'],  # transition standard -> PBO -> IG2
}

_DISCARD_CONFIG = {
    'Expiration_Period_Distribution': 'EXPONENTIAL_DISTRIBUTION',
    'Expiration_Period_Exponential': 730,
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

def update_serialize(simulation, larval_multiplier, species, serialization=0, sim_duration=40 * 365,
                     serialized_population_path_df=None):
    if serialization:
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Time_Steps = [sim_duration]
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'TIMESTEP'
        sweep_sim_larval_capacity(simulation, species, value=larval_multiplier)

    else:
        serialized_population_path = serialized_population_path_df[
            (serialized_population_path_df['Larval_Capacity_' + species[0]] == larval_multiplier[0])
            & (serialized_population_path_df['Larval_Capacity_' + species[1]] == larval_multiplier[1])][
            'Outpath'].values[0]
        simulation.task.config.parameters.Simulation_Duration = sim_duration
        simulation.task.config.parameters.Serialization_Mask_Node_Read = 0
        simulation.task.config.parameters.Serialization_Mask_Node_Write = 0
        simulation.task.config.parameters.Serialized_Population_Path = os.path.join(serialized_population_path,
                                                                                    'output')
        simulation.task.config.parameters.Serialized_Population_Reading_Type = 'READ'
        simulation.task.config.parameters.Serialized_Population_Writing_Type = 'NONE'
        simulation.task.config.parameters.Serialized_Population_Filenames = ['state-14600.dtk']
        sweep_sim_larval_capacity(simulation, species, value=larval_multiplier)

    lh_dict = dict(zip([f'Larval_Capacity_{s}' for s in species], larval_multiplier))
    return_dict = {"Serialization": serialization}
    return_dict.update(lh_dict)
    return return_dict


def add_multi_insecticide_bednet(campaign, start_day, coverage, components, discard_config):
    """Add a MultiInsecticideUsageDependentBednet event to the campaign.

    components: list of dicts, each with:
        name        - insecticide name (must match a name in config Insecticides)
        killing     - Initial_Effect for Killing_Config
        kill_decay  - Decay_Time_Constant (days) for Killing_Config
        blocking    - Initial_Effect for Blocking_Config (omit or 0 for no blocking)
        block_decay - Decay_Time_Constant (days) for Blocking_Config (used if blocking > 0)
    discard_config: dict mapping Expiration_Period_* parameter names to values
    """
    schema_path = campaign.schema_path
    intervention = s2c.get_class_with_defaults("MultiInsecticideUsageDependentBednet", schema_path)

    insecticides = []
    for comp in components:
        item = s2c.get_class_with_defaults("idmType:InsecticideWaningEffect", schema_path)
        item.Insecticide_Name = comp['name']
        item.Killing_Config = waning_utils.get_waning_from_parameters(
            schema_path, initial=comp['killing'],
            box_duration=0, decay_time_constant=comp['kill_decay'])
        blocking = comp.get('blocking', 0)
        if blocking > 0:
            item.Blocking_Config = waning_utils.get_waning_from_parameters(
                schema_path, initial=blocking,
                box_duration=0, decay_time_constant=comp['block_decay'])
        else:
            item.Blocking_Config = waning_utils.get_waning_from_parameters(
                schema_path, initial=0, box_duration=-1)
        item.Repelling_Config = waning_utils.get_waning_from_parameters(
            schema_path, initial=0, box_duration=-1)
        insecticides.append(item)

    intervention.Insecticides = insecticides

    intervention.Usage_Config_List = [
        _get_seasonal_times_and_values(campaign, None),
        _get_age_times_and_values(campaign, None),
    ]

    intervention.Received_Event = campaign.get_send_trigger("Bednet_Got_New_One", old=True)
    intervention.Using_Event    = campaign.get_send_trigger("Bednet_Using", old=True)
    intervention.Discard_Event  = campaign.get_send_trigger("Bednet_Discarded", old=True)

    for param, val in discard_config.items():
        setattr(intervention, param, val)

    malaria_common.add_campaign_event(campaign=campaign,
                                      start_day=start_day,
                                      demographic_coverage=coverage,
                                      individual_intervention=intervention)

def update_camp_type(simulation, sim_params, serialize=0, sim_duration=40 * 365):

    build_camp_partial = partial(build_camp, sim_params=sim_params, serialize=serialize, sim_duration=sim_duration)

    simulation.task.create_campaign_from_callback(build_camp_partial)

    scenario = sim_params[6] if len(sim_params) > 6 else 'STD-STD-STD'
    return {"Baseline": sim_params[0], 'Start_Day': sim_params[1], 'Coverage': sim_params[2],
            'Start_Day_IRS': sim_params[3], 'Coverage_IRS': sim_params[4],
            'Scenario': scenario}
def build_camp(sim_params, serialize=0, sim_duration=40 * 365):
    """Build the campaign for a scenario using emod_api.

    Intervention (serialize == 0, not baseline): one net per 3-year block following
    the scenario's net sequence (e.g. STD -> PBO -> IG2). Genotype-specific resistance
    is applied separately via the per-locus Killing_Modifiers in set_param_fn.
    """
    camp.set_schema(manifest.schema_file)

    if not serialize:
        add_treatment_seeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.4, "agemin": 15,
                                              "agemax": 70, "rate": 0.3},
                                             {"trigger": "NewSevereCase", "coverage": 0.8, "rate": 0.5}],
                              drug=['Artemether', 'Lumefantrine'],
                              start_day=0,
                              broadcast_event_name='Received_Treatment'
                              )

        # Deploy the scenario's net sequence: one net per 3-year block.
        if not sim_params[0]:  # not a no-intervention baseline
            scenario = sim_params[6] if len(sim_params) > 6 else 'STD-STD-STD'
            net_sequence = SCENARIOS[scenario]
            start_day = sim_params[1]
            coverage = sim_params[2]
            n_blocks = int(np.ceil(sim_duration / (3 * 365)))
            for m in range(n_blocks):
                if m >= len(net_sequence):
                    break
                net_type = net_sequence[m]
                deploy_day = start_day + 3 * m * 365
                components = NET_COMPONENTS[net_type]
                if components is None:
                    # STD: single-insecticide path
                    add_scheduled_usage_dependent_bednet(
                        camp,
                        start_day=deploy_day,
                        demographic_coverage=coverage,
                        insecticide=_STD_SPEC['insecticide'],
                        discard_config=_DISCARD_CONFIG,
                        blocking_initial_effect=_STD_SPEC['blocking_initial_effect'],
                        blocking_decay_time_constant=_STD_SPEC['blocking_decay_time_constant'],
                        killing_initial_effect=_STD_SPEC['killing_initial_effect'],
                        killing_decay_time_constant=_STD_SPEC['killing_decay_time_constant'],
                    )
                else:
                    # PBO / IG2: dual-insecticide path
                    add_multi_insecticide_bednet(
                        camp,
                        start_day=deploy_day,
                        coverage=coverage,
                        components=components,
                        discard_config=_DISCARD_CONFIG,
                    )

    else:
        add_treatment_seeking(camp, targets=[{"trigger": "NewClinicalCase", "coverage": 0.4, "agemin": 15,
                                              "agemax": 70, "rate": 0.3},
                                             {"trigger": "NewSevereCase", "coverage": 0.8, "rate": 0.5}],
                              drug=['Artemether', 'Lumefantrine'],
                              start_day=sim_duration - 10 * 365,
                              broadcast_event_name='Received_Treatment'
                              )
    return camp


def build_demog():
    """
    Build a demographics input file for the DTK using emod_api.
    """
    demog = Demographics.from_file(
        "../input_files/single_node_demographics.json")
    return demog
def set_param_fn(config):
    """
    Callback passed to emod-api.config to set parameters.

    Sets up the three independent resistance loci (a/b/c), the per-net resistance
    Killing_Modifiers, and the two vector species (An. funestus, An. arabiensis).
    """
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
    # species must be added before insecticide resistance is configured for them
    vecconf.add_species(config, manifest, species_to_select=["funestus"])
    vecconf.set_species_param(config, species='funestus', parameter="Habitats",
                              value=[vecconf.configure_linear_spline(manifest, max_larval_capacity=pow(10, 6),
                                                                     capacity_distribution_number_of_years=1,
                                                                     capacity_distribution_over_time={
                                                                         "Times": [0.0, 30.417, 60.833, 91.25, 121.667,
                                                                                   152.083, 182.5, 212.917, 243.333,
                                                                                   273.75, 304.167, 334.583],
                                                                         "Values": [0.16, .11, .25, .3, .35, 0.7, .3,
                                                                                    .2, .25, .2, .1, .15]})],
                              overwrite=True)

    # three independent resistance loci, added in order a, b, c; this fixes the
    # genome-string locus order to gender-a-b-c, which the analyzer relies on
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
                                                                         "Values": [.29, .32, .35, .5, 1., .48, .42,
                                                                                    .35, .092, .08, .31, .19]})],

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


    # per-net, per-locus insecticide resistance
    # added after both species exist in the config (required by emodpy_malaria);
    # each net keys its Killing_Modifier off its own locus
    # pyrethroid resistance (locus a)
    for sp in ['funestus', 'arabiensis']:

        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='pyrethroid',
            species=sp,
            allele_combo=[['a1', 'a1']],
            killing=STD_RESISTANCE[sp]['homozygous']
        )

        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='pyrethroid',
            species=sp,
            allele_combo=[['a0', 'a1']],
            killing=STD_RESISTANCE[sp]['heterozygous']
        )

    # PBO resistance (locus b)
    for sp in ['funestus', 'arabiensis']:
        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='PBO-synergist',
            species=sp,
            allele_combo=[['b1', 'b1']],
            killing=PBO_RESISTANCE[sp]['homozygous']
        )

        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='PBO-synergist',
            species=sp,
            allele_combo=[['b0', 'b1']],
            killing=PBO_RESISTANCE[sp]['heterozygous']
        )

    # IG2 resistance (locus c)
    for sp in ['funestus', 'arabiensis']:
        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='chlorfenapyr',
            species=sp,
            allele_combo=[['c1', 'c1']],
            killing=IG2_RESISTANCE[sp]['homozygous']
        )

        add_insecticide_resistance(
            config,
            manifest,
            insecticide_name='chlorfenapyr',
            species=sp,
            allele_combo=[['c0', 'c1']],
            killing=IG2_RESISTANCE[sp]['heterozygous']
        )

    return config

def rvg_genome_config_builder2(params):
    params.Include_Vector_State_Columns = False
    params.Include_Death_By_State_Columns = False
    params.Species = 'funestus'
    params.Stratify_By = 'GENOME'
    params.Combine_Similar_Genomes = True
    params.Start_Day = 0
    params.End_Day = 9 * 365
    return params


def rvg_genome_config_builder1(params):
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
    params.Reporting_Interval = 1
    params.Filename_Suffix = 'Monthly'
    params.Max_Number_Reports = 9 * 365
    params.End_Day = 9 * 365
    return params


def rmf_config_builder(params):
    params.Start_Day = 1
    params.End_Day = 9 * 365
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
