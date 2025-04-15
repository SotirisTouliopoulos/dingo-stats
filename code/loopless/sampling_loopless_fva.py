
"""
This script performs sampling of the E. coli core model using Dingo and optGP with different 
constraints based on flux variability analysis (FVA), aiming to remove loops.

It generates samples for three different cases:
1. Default constraints
2. Constraints based on loopless FVA of the FRD7 reaction
3. Constraints based on loopless FVA of both the FRD7 and SUCDi reactions

The script saves the generated samples to pickle files for each case. 

Run this file with: python3 dingo_sampling_loopless_fva.py <optimal_percentage>
"""



# Load libraries
import sys
from dingo import MetabolicNetwork, PolytopeSampler
import numpy as np
import pickle
from cobra.flux_analysis import flux_variability_analysis
from cobra.io import load_json_model, read_sbml_model, load_matlab_model
import cobra



# Load model and create deep copies
ec_model = read_sbml_model("/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/models/e_coli_core.xml")
reaction_list = ec_model.reactions

default_fba_solution = ec_model.optimize()

ec_model_loopless_fva_FRD7 = ec_model.copy()
ec_model_loopless_fva_FRD7_SUCDi = ec_model.copy()



# Define optimal percentage
optimal_percentage = float(sys.argv[1])
print(f"Optimal percentage: {optimal_percentage}")
fraction_of_optimum = float(optimal_percentage / 100)
print(f"Fraction of optimum: {fraction_of_optimum}")
minimum_objective_value = default_fba_solution.objective_value * fraction_of_optimum

# Perform flux variability analysis
default_fva = flux_variability_analysis(ec_model, reaction_list=reaction_list, fraction_of_optimum=fraction_of_optimum, loopless=False)
fva_dict = default_fva.apply(lambda row: (row['minimum'], row['maximum']), axis=1).to_dict()

# Perform Loopless flux variability analysis
loopless_fva = flux_variability_analysis(ec_model, reaction_list=reaction_list, fraction_of_optimum=fraction_of_optimum, loopless=True)
fva_loopless_dict = loopless_fva.apply(lambda row: (row['minimum'], row['maximum']), axis=1).to_dict()



# Check the results from flux variability analysis
print("Default FRD7 flux range", fva_dict['FRD7'], "Loopless FRD7 flux range", fva_loopless_dict['FRD7'])
print("Default SUCDi flux range", fva_dict['SUCDi'], "Loopless SUCDi flux range", fva_loopless_dict['SUCDi'])



# Constraint FRD7 with the flux variability analysis results
ec_model_loopless_fva_FRD7.reactions.FRD7.bounds = fva_loopless_dict['FRD7']

# Constraint FRD7 and SUCDi with the flux variability analysis results
ec_model_loopless_fva_FRD7_SUCDi.reactions.FRD7.bounds = fva_loopless_dict['FRD7']
ec_model_loopless_fva_FRD7_SUCDi.reactions.SUCDi.bounds = fva_loopless_dict['SUCDi']



# Dingo sampling with default constraints
ec_dingo_model = MetabolicNetwork.from_cobra_model(ec_model)
ec_dingo_model.set_opt_percentage(optimal_percentage)
sampler_dingo_default = PolytopeSampler(ec_dingo_model)
samples_dingo_default = sampler_dingo_default.generate_steady_states()

# Dingo sampling with FRD7 constraints
ec_loopless_fva_FRD7_dingo_model = MetabolicNetwork.from_cobra_model(ec_model_loopless_fva_FRD7)
ec_loopless_fva_FRD7_dingo_model.set_opt_percentage(optimal_percentage)
sampler_dingo_loopless_FRD7 = PolytopeSampler(ec_loopless_fva_FRD7_dingo_model)
samples_dingo_loopless_FRD7 = sampler_dingo_loopless_FRD7.generate_steady_states()

# Dingo sampling with FRD7 and SUCDi constraints
ec_loopless_fva_FRD7_SUCDi_dingo_model = MetabolicNetwork.from_cobra_model(ec_model_loopless_fva_FRD7_SUCDi)
ec_loopless_fva_FRD7_SUCDi_dingo_model.set_opt_percentage(optimal_percentage)
sampler_dingo_loopless_FRD7_SUCDi = PolytopeSampler(ec_loopless_fva_FRD7_SUCDi_dingo_model)
samples_dingo_loopless_FRD7_SUCDi = sampler_dingo_loopless_FRD7_SUCDi.generate_steady_states()



# optGP sampling with default constraints
ec_model.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM").lower_bound = minimum_objective_value
sampler_optgp_default = cobra.sampling.optgp.OptGPSampler(ec_model)
samples_optgp_default = sampler_optgp_default.sample(n=3000)
samples_optgp_default = samples_optgp_default.to_numpy().T

# optGP sampling with FRD7 constraints
ec_model_loopless_fva_FRD7.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM").lower_bound = minimum_objective_value
sampler_optgp_loopless_fva_FRD7 = cobra.sampling.optgp.OptGPSampler(ec_model_loopless_fva_FRD7)
samples_optgp_loopless_fva_FRD7 = sampler_optgp_loopless_fva_FRD7.sample(n=3000)
samples_optgp_loopless_fva_FRD7 = samples_optgp_loopless_fva_FRD7.to_numpy().T

# optGP sampling with FRD7 and SUCDi constraints
ec_model_loopless_fva_FRD7_SUCDi.reactions.get_by_id("BIOMASS_Ecoli_core_w_GAM").lower_bound = minimum_objective_value
sampler_optgp_loopless_fva_FRD7_SUCDi = cobra.sampling.optgp.OptGPSampler(ec_model_loopless_fva_FRD7_SUCDi)
samples_optgp_loopless_fva_FRD7_SUCDi = sampler_optgp_loopless_fva_FRD7_SUCDi.sample(n=3000)
samples_optgp_loopless_fva_FRD7_SUCDi = samples_optgp_loopless_fva_FRD7_SUCDi.to_numpy().T



# Save dingo samples to pickle
file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/dingo_samples/samples_dingo_ec_default_opt_" + str(optimal_percentage) + ".pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(samples_dingo_default, hopsy_samples_file)
    
file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/dingo_samples/samples_dingo_ec_loopless_FRD7_opt_" + str(optimal_percentage) + ".pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(samples_dingo_loopless_FRD7, hopsy_samples_file)

file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/dingo_samples/samples_dingo_ec_loopless_FRD7_SUCDi_opt_" + str(optimal_percentage) + ".pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(samples_dingo_loopless_FRD7_SUCDi, hopsy_samples_file)
    


# Save optGP samples to pickle
file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/optgp_samples/samples_optgp_ec_default_opt_" + str(optimal_percentage) + ".pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(samples_optgp_default, hopsy_samples_file)
    
file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/optgp_samples/samples_optgp_ec_loopless_FRD7_opt_" + str(optimal_percentage) + ".pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(samples_optgp_loopless_fva_FRD7, hopsy_samples_file)

file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/ext_data/optgp_samples/samples_optgp_ec_loopless_FRD7_SUCDi_opt_" + str(optimal_percentage) + ".pckl"
with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(samples_optgp_loopless_fva_FRD7_SUCDi, hopsy_samples_file)