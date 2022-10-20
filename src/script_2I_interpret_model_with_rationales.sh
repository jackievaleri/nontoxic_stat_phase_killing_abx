# conda activate chemprop
# then bash script_2I_interpret_model_with_rationales.sh
chemprop_interpret --data_path data/training_plus_ml_curation_hits.csv --checkpoint_dir models/primary_plus_ML_val_models/32FINAL/ --property_id 1 --smiles_column SMILES --max_atoms 20 --min_atoms 8 --prop_delta 0.1 --features_generator rdkit_2d_normalized --no_features_scaling > out/2I_interpretation_results_inhibit.txt
chemprop_interpret --data_path data/training_plus_ml_curation_hits.csv --checkpoint_dir models/primary_plus_ML_val_models/32FINAL/ --property_id 2 --smiles_column SMILES --max_atoms 20 --min_atoms 8 --prop_delta 0.1 --features_generator rdkit_2d_normalized --no_features_scaling > out/2I_interpretation_results_killing.txt

# to use interpret, needed to downgrade rdkit (as per https://github.com/chemprop/chemprop/issues/178)
#conda activate chemprop
#conda install -c conda-forge rdkit=2019.09.1