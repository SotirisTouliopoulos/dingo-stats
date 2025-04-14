
import pickle
import glob
import os
import numpy as np


# list all files starting with "RI_steady_states_opt" in the directory

"""
directory = 'data/Steady_States_Samples/RI_no_mucins'

files = glob.glob(os.path.join(directory, "RI_steady_states_opt*"))

minimum_samples = float("inf")
for file in files:
    with open(file, "rb") as absolute_steady_states:
        absolute_steady_states = abs(pickle.load(absolute_steady_states))
        print(file, absolute_steady_states.shape[1])
        if absolute_steady_states.shape[1] < minimum_samples:
            minimum_samples = absolute_steady_states.shape[1]
            
    del absolute_steady_states

# 23200
print("Minimum samples: ", minimum_samples)
"""


"""
directory = 'data/Steady_States_Samples/RI_no_mucins'

files = glob.glob(os.path.join(directory, "RI_steady_states_opt_100*"))

minimum_samples = 23200
for file in files:
    with open(file, "rb") as absolute_steady_states:
        absolute_steady_states = abs(pickle.load(absolute_steady_states))
        samples = absolute_steady_states.shape[1]
        step = int(samples / minimum_samples)
        
        subsample_indices = np.arange(0, samples, step)
        print(subsample_indices)
        
        subsamples_steady_states = absolute_steady_states[:, subsample_indices]
        print(subsamples_steady_states.shape)
        
    del absolute_steady_states
"""




# rm data/Steady_States_Samples/RI_no_mucins/*subsampled.pckl
directory = 'data/Steady_States_Samples/RI_no_mucins'

files = glob.glob(os.path.join(directory, "RI_steady_states_opt_80_*"))

minimum_samples = 23200
for file in files:
    with open(file, "rb") as absolute_steady_states:
        absolute_steady_states = abs(pickle.load(absolute_steady_states))
        samples = absolute_steady_states.shape[1]
        
        step = samples / minimum_samples
        
        subsample_indices = np.arange(0, samples, step).astype(int)
        
        # force exactly the same number of samples
        # this may take duplicate samples
        subsample_indices = subsample_indices[:minimum_samples]
        
        subsamples_steady_states = absolute_steady_states[:, subsample_indices]
        print(file, samples, subsamples_steady_states.shape[1])
        
        # extract subsamples to new file
        #file_name_subsamples = "/home/touliopoulos/project/erasmus/erasmus_2025_project/data" + file[4:-5] + "_subsampled.pckl"
        #with open(file_name_subsamples, "wb") as hopsy_samples_file: 
        #    pickle.dump(subsamples_steady_states, hopsy_samples_file)
        
    del absolute_steady_states, subsamples_steady_states




