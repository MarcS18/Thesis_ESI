Energy_Analysis_CrossMiner.py: Calculation of activation energies from compiutational output files.

TS_input_gen.py: Generation of computational input files.

TS_mol_gen_xyz: Generation of molSimplify input files for generation of structures. Uses .xyz/.mol files only.

deprotonation_index_correction.py: Correction of coordinating atom indexes between ligand from CrossMiner and ligand in complex due to deprotonation of the ligand during complex generation.

get_descriptors.py: Extraction of descriptors from computational output files for machine learning models. Workflow Only.

get_descriptors_flexible.py: Extraction of descriptors from computational output files for machine learning models. Flexible for use with different computational methods. Using in DFT comparisons.

get_descriptors_nots_pip.py: Extraction of descriptors from computational outfiles without using the transition state structures. Flexible for use with different computational methods. Workflow and comparisons.

ts_check.py: Automatic checking for correct transition states.
