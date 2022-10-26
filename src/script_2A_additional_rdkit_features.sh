# activate venv
source activate chemprop

# navigate into the scripts folder to use features generator - path may change for you
cd ../../chemprop/scripts

# save features for all the drugs in screen
python save_features.py --data_path ../../nontoxic_stat_phase_killing_abx/out/data_prep_for_ml_fullset.csv --features_generator rdkit_2d --save_path ../../nontoxic_stat_phase_killing_abx/out/data_prep_for_ml_fullset.npz --smiles_column SMILES

# save features for all 800K compounds - not included for space reasons
python save_features.py --data_path ../../nontoxic_stat_phase_killing_abx/data/broad800k.csv --features_generator rdkit_2d --save_path ../../nontoxic_stat_phase_killing_abx/data/broad.npz --smiles_column smiles
