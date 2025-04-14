
from dingo import MetabolicNetwork, PolytopeSampler
import numpy as np
import pickle
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import load_json_model, read_sbml_model, load_matlab_model


ec_model = read_sbml_model("/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/models/e_coli_core.xml")
reaction_list = ec_model.reactions

ec_model_loopless_fva_FRD7 = ec_model.copy()
ec_model_loopless_fva_FRD7_SUCDi = ec_model.copy()

# for 100 opt. it returns 0.0-0.0 for FRD7 and 5.064376-5.064376 for SUCDi
default_fva = flux_variability_analysis(ec_model, reaction_list=reaction_list, fraction_of_optimum=1, loopless=False)
loopless_fva = flux_variability_analysis(ec_model, reaction_list=reaction_list, fraction_of_optimum=1, loopless=True)


# Sample with default constraints
ec_dingo_model = MetabolicNetwork.from_cobra_model(ec_model)
ec_dingo_model.set_opt_percentage(70)
sampler_default = PolytopeSampler(ec_dingo_model)
steady_states_default = sampler_default.generate_steady_states()


# Sample with loopless constraints on FRD7
ec_model_loopless_fva_FRD7.reactions.FRD7.lower_bound = 0
ec_model_loopless_fva_FRD7.reactions.FRD7.upper_bound = 0

ec_loopless_fva_FRD7_dingo_model = MetabolicNetwork.from_cobra_model(ec_model_loopless_fva_FRD7)
ec_loopless_fva_FRD7_dingo_model.set_opt_percentage(70)
sampler_loopless_FRD7 = PolytopeSampler(ec_loopless_fva_FRD7_dingo_model)
steady_states_loopless_FRD7 = sampler_loopless_FRD7.generate_steady_states()


# Sample with loopless constraints on FRD7 and SUCDi
ec_model_loopless_fva_FRD7_SUCDi.reactions.FRD7.lower_bound = 0
ec_model_loopless_fva_FRD7_SUCDi.reactions.FRD7.upper_bound = 0
ec_model_loopless_fva_FRD7_SUCDi.reactions.SUCDi.lower_bound = 5.064376
ec_model_loopless_fva_FRD7_SUCDi.reactions.SUCDi.upper_bound = 5.064376

ec_loopless_fva_FRD7_SUCDi_dingo_model = MetabolicNetwork.from_cobra_model(ec_model_loopless_fva_FRD7_SUCDi)
ec_loopless_fva_FRD7_SUCDi_dingo_model.set_opt_percentage(70)
sampler_loopless_FRD7_SUCDi = PolytopeSampler(ec_loopless_fva_FRD7_SUCDi_dingo_model)
steady_states_loopless_FRD7_SUCDi = sampler_loopless_FRD7_SUCDi.generate_steady_states()


# Save samples to pickle
file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/dingo_samples/samples_dingo_ec_default_opt_70.pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(steady_states_default, hopsy_samples_file)
    
    
file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/dingo_samples/samples_dingo_ec_loopless_FRD7_opt_70.pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(steady_states_loopless_FRD7, hopsy_samples_file)
    
    
file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/dingo_samples/samples_dingo_ec_loopless_FRD7_SUCDi_opt_70.pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(steady_states_loopless_FRD7_SUCDi, hopsy_samples_file)