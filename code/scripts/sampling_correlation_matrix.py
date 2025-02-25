# %% this is an interactive chunk 

from dingo import MetabolicNetwork, PolytopeSampler
import numpy as np
import pickle
from cobra.io import load_json_model, read_sbml_model, load_matlab_model

print("libraries loaded.")

# %% 

#roseburia_dingo_model = MetabolicNetwork.from_mat("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_withMucins.mat")
#roseburia_dingo_model = MetabolicNetwork.from_mat("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat")
#roseburia_dingo_model = MetabolicNetwork.from_sbml("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/e_coli_core.xml")

#roseburia_cobra_model = read_sbml_model("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/e_coli_core.xml")

#reactions = roseburia_cobra_model.reactions
#for reaction in reactions:
#    print(reaction)



for i in (0.01,1,20,40,60,80,100):
#for i in (20, 40, 60, 80):
    roseburia_dingo_model = MetabolicNetwork.from_mat("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/Roseburia_intestinalis_L1_82_fixed.mat")
    #roseburia_dingo_model = MetabolicNetwork.from_sbml("/home/touliopoulos/project/erasmus/erasmus_2025_project/data/Reconstructions_data/e_coli_core.xml")

    reactions = roseburia_dingo_model.reactions

    roseburia_dingo_model.set_opt_percentage(i)
    
    sampler = PolytopeSampler(roseburia_dingo_model)
    steady_states = sampler.generate_steady_states()

    #corr_matrix = np.corrcoef(steady_states, rowvar=True)

    file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/RI_steady_states_opt_" + str(i) + "_002.pckl"
    #file_name_samples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/EC_steady_states_opt_" + str(i) + ".pckl"

    with open(file_name_samples, "wb") as hopsy_samples_file: 
        pickle.dump(steady_states, hopsy_samples_file)
        
    #file_name_correlations = "/home/touliopoulos/project/erasmus/erasmus_2025_project/RI_corr_matrix_opt_" + str(i) + ".pckl"
    #file_name_correlations = "/home/touliopoulos/project/erasmus/erasmus_2025_project/EC_corr_matrix_opt_" + str(i) + ".pckl"

    #with open(file_name_correlations, "wb") as hopsy_correlations_file: 
    #    pickle.dump(corr_matrix, hopsy_correlations_file)
    
        
    #del steady_states, corr_matrix
    del steady_states






#for row in corr_matrix:
#    print(" ".join(map(str, row)))
    

# To save samples
#with open("/home/touliopoulos/project/erasmus/erasmus_2025_project/RI_steady_states_opt_1.pckl", "wb") as hopsy_samples_file: 
#    pickle.dump(steady_states, hopsy_samples_file)


#with open("/home/touliopoulos/project/erasmus/erasmus_2025_project/RI_corr_matrix_opt_1.pckl", "wb") as hopsy_correlations_file: 
#    pickle.dump(corr_matrix, hopsy_correlations_file)



# To load samples
#with open("my_samples.pckl", "rb") as f:
#    obj = pickle.load(f)