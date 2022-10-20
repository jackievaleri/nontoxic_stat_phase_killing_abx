# activate venv
source activate chemprop

# navigate into the scripts folder to use features generator
cd ../../../beta_lactam/models/chemprop/scripts

# save features for all the drugs in screen
#python save_features.py --data_path ../../../../ian_stat_phase_ML/data/round3finalval08162021_FULL.csv --features_generator rdkit_2d --save_path ../../../../ian_stat_phase_ML/data/round3finalval08162021_FULL_unnorm_fts_for_pub.npz --smiles_column SMILES

# and also save the validated compounds from the 800k library
python save_features.py --data_path ../../../../ian_stat_phase_ML/toxicity_05102022_felixtox/out/train_val_hits.csv --features_generator rdkit_2d --save_path ../../../../ian_stat_phase_ML/data/round3finalval08162021_hits_plus_16val_hits_FULL_unnorm_fts_for_pub.npz --smiles_column SMILES
