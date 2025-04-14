
import cobra
from dingo import MetabolicNetwork, PolytopeSampler
import numpy as np
import pickle
from cobra.io import load_json_model, read_sbml_model, load_matlab_model



"""
#roseburia_dingo_model = MetabolicNetwork.from_mat("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_withMucins.mat")
#roseburia_dingo_model = MetabolicNetwork.from_mat("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat")
#roseburia_dingo_model = MetabolicNetwork.from_sbml("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/e_coli_core.xml")

#roseburia_cobra_model = read_sbml_model("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/e_coli_core.xml")


for i in (99, 99, 99):
#for i in (20, 40, 60, 80):
    
    roseburia_dingo_model = MetabolicNetwork.from_mat("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat")
    #roseburia_dingo_model = MetabolicNetwork.from_sbml("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/e_coli_core.xml")

    reactions = roseburia_dingo_model.reactions

    roseburia_dingo_model.set_opt_percentage(i)
    
    sampler = PolytopeSampler(roseburia_dingo_model)
    steady_states = sampler.generate_steady_states()

    file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Steady_States_Samples/RI_no_mucins/RI_steady_states_opt_" + str(i) + "_rep_002.pckl"
    #file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/EC_steady_states_opt_" + str(i) + ".pckl"

    with open(file_name_samples, "wb") as hopsy_samples_file: 
        pickle.dump(steady_states, hopsy_samples_file)
        
    del steady_states
""" 



"""
modelfile = "/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat"

roseburia_dingo_model = MetabolicNetwork.from_mat(modelfile)
roseburia_dingo_model.set_opt_percentage(100)

sampler = PolytopeSampler(roseburia_dingo_model)
steady_states = sampler.generate_steady_states()
print(np.max(steady_states))

#file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Steady_States_Samples/RI_no_mucins/RI_steady_states_opt_" + str(i) + "_rep_003.pckl"
file_name_samples = "/home/touliopoulos/100_rep_3.pckl"

with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(steady_states, hopsy_samples_file)
"""



from cobra.sampling import sample
roseburia_cobra_model = load_matlab_model("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat")
samples_df = sample(roseburia_cobra_model, 20000)
samples_np = samples_df.to_numpy(dtype=float)
samples_np = samples_np.T

print(samples_np.shape)
print(np.max(samples_np[2944,]))
print(np.min(samples_np[2944,]))
print(np.max(samples_np))

file_name_samples = "/home/touliopoulos/100_cobrapy_rep_3.pckl"

with open(file_name_samples, "wb") as hopsy_samples_file: 
    pickle.dump(samples_np, hopsy_samples_file)




