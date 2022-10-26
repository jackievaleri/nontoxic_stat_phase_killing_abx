# conda activate chemprop
# then bash [name of this file].sh
chemprop_interpret --data_path data/training_plus_ml_curation_hits.csv --checkpoint_dir ../felix_cytotox_models/final_tox_hepg2/ --property_id 1 --smiles_column SMILES --max_atoms 20 --min_atoms 8 --prop_delta 0.1 --features_generator rdkit_2d_normalized --no_features_scaling > out/3B_toxicity_rationales_results_from_HepG2_model.txt

# to use interpret, needed to downgrade rdkit (as per https://github.com/chemprop/chemprop/issues/178)
#conda activate chemprop
#conda install -c conda-forge rdkit=2019.09.1