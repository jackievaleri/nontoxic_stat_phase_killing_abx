cd ../../chemprop/;


# ffn
mkdir ../../nontoxic_stat_phase_killing_abx/models/final_hyperopt_shallow_models_for_pub/FINALFFN26; python train.py --ffn_num_layers 3 --hidden_size 1500 --dropout 0.3 --save_dir ../../nontoxic_stat_phase_killing_abx/models/final_hyperopt_shallow_models_for_pub/FINALFFN26 --data_path ../../nontoxic_stat_phase_killing_abx/out/data_for_sklearn.csv --num_folds 20 --dataset_type classification --features_generator morgan --no_features_scaling --split_type scaffold_balanced --split_sizes 0.8 0.1 0.1 --smiles_columns SMILES --target_columns hit_inh hit_kill --depth 0 --features_only --metric auc --extra_metrics prc-auc

# rfc
mkdir ../../nontoxic_stat_phase_killing_abx/models/final_hyperopt_shallow_models_for_pub/FINALRFC76; python sklearn_train.py --num_bits 4096 --radius 3 --class_weight balanced --num_trees 2000 --save_dir ../../nontoxic_stat_phase_killing_abx/models/final_hyperopt_shallow_models_for_pub/FINALRFC76 --data_path ../../nontoxic_stat_phase_killing_abx/out/data_for_sklearn.csv --num_folds 20 --dataset_type classification --features_path ../../nontoxic_stat_phase_killing_abx/out/data_prep_for_ml_fullset.npz  --no_features_scaling --split_type scaffold_balanced --split_sizes 0.8 0.1 0.1 --smiles_columns SMILES --target_columns hit_inh hit_kill --model_type random_forest --metric auc --extra_metrics prc-auc

# svm inh
mkdir ../../nontoxic_stat_phase_killing_abx/models/final_hyperopt_shallow_models_for_pub/FINALSVMINH63; python sklearn_train.py --num_bits 4096 --radius 2 --class_weight balanced --save_dir ../../nontoxic_stat_phase_killing_abx/models/final_hyperopt_shallow_models_for_pub/FINALSVMINH63 --data_path ../../nontoxic_stat_phase_killing_abx/out/data_for_sklearn.csv --num_folds 20 --dataset_type classification --features_path ../../nontoxic_stat_phase_killing_abx/data/out.npz  --no_features_scaling --split_type scaffold_balanced --split_sizes 0.8 0.1 0.1 --smiles_columns SMILES --target_columns hit_inh --model_type svm --metric auc --extra_metrics prc-auc

# svm kill
mkdir ../../nontoxic_stat_phase_killing_abx/models/final_hyperopt_shallow_models_for_pub/FINALSVMKILL70; python sklearn_train.py --num_bits 4096 --radius 3 --class_weight balanced --save_dir ../../nontoxic_stat_phase_killing_abx/models/final_hyperopt_shallow_models_for_pub/FINALSVMKILL70 --data_path ../../nontoxic_stat_phase_killing_abx/out/data_for_sklearn.csv --num_folds 20 --dataset_type classification --features_path ../../nontoxic_stat_phase_killing_abx/out/data_prep_for_ml_fullset.npz  --no_features_scaling --split_type scaffold_balanced --split_sizes 0.8 0.1 0.1 --smiles_columns SMILES --target_columns hit_kill --model_type svm --metric auc --extra_metrics prc-auc
