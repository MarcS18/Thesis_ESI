Energy_Analysis_CrossMiner.py: Calculation of activation energies from compiutational output files. Workflow Only.

Energy_Analysis_CrossMiner_flexible.py: Calculation of activation energies from compiutational output files. Flexible for use with different DFT methods. Used in comparisons.

TS_input_gen.py: Generation of computational input files. Uses only ORCA.

TS_mols_gen_xyz.py: Generation of molSimplify input files for generation of structures. Uses .xyz/.mol files only.

deprotonation_index_correction.py: Correction of coordinating atom indexes between ligand from CrossMiner and ligand in complex due to deprotonation of the ligand during complex generation.

get_descriptors.py: Extraction of descriptors from computational output files for machine learning models. Workflow Only.

get_descriptors_flexible_expanded.py: Extraction of descriptors from computational output files for machine learning models. Expanded to explain strain on the nucleophile. Flexible for use with different computational methods. Using in DFT comparisons.

get_descriptors_nots_pyr_expanded.py: Extraction of descriptors from computational outfiles without using the transition state structures. Flexible for use with different computational methods. Workflow and comparisons.

ts_check.py: Automatic checking for correct transition states.
