import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import subprocess

import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdmolops
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.ML.Cluster import Butina
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from sklearn.cluster import AgglomerativeClustering
from rdkit.Chem.Draw import rdMolDraw2D

# shut off warnings
from rdkit import RDLogger                                                                                                                                                               
RDLogger.DisableLog('rdApp.*') 

# get properties from this package: https://github.com/ikmckenz/adme-pred-py
from adme_pred import ADME # pip install ADME_predict is NOT the right one
# directly downloaded adme-pred-py into my src folder

## Tanimoto Similarity Scoring ## 

# function specifically to add info to a pandas dataframe for easy display to experimental team
def print_drug_dict(drugdict, smiles_list, col_list, names_list, merging_df):

    ret = pd.DataFrame(columns = col_list)
    for key, item in drugdict.items():
        tans, indexes = item
        for tan, index in zip(tans, indexes):
            try: 
                if names_list == None:
                    ret = ret.append(pd.DataFrame([[key, tan, smiles_list[index], index]], columns = col_list))
                else:
                    ret = ret.append(pd.DataFrame([[key, tan, smiles_list[index], names_list[index]]], columns = col_list))
            except Exception as e:
                print(e)
                continue
    ret = ret.merge(merging_df, left_on = 'SMILES', right_on = 'SMILES')
    return(ret)

# for every molecule, get similarity to closest antibiotic
def get_lowest_tanimoto_from_drug_set(new_set, abx_fps, smiles_list, col_list, merging_df, names_list = None):
    mols = [Chem.MolFromSmiles(x) for x in new_set]
    
    best_similarity = {}
    index = 0
   
    for m in mols:
        try:
            mol = Chem.RDKFingerprint(m)
        except Exception as e:
            index = index + 1
            continue
        set_index = 0
        curr_highest_sim = 0
        curr_highest_drug = None
        for abx in abx_fps:
            sim = DataStructs.FingerprintSimilarity(mol,abx)
            if sim >= curr_highest_sim: # use greater or equal to in order to handle weird mols (for example, CH4) that have no tan sim to anything
                curr_highest_sim = sim
                curr_highest_drug = set_index
            set_index = set_index + 1
        if curr_highest_drug == None:
            display(m)
        info_to_include = [[curr_highest_sim], [curr_highest_drug]]
        best_similarity[new_set[index]] = info_to_include
        
        index = index + 1
    
    ret = print_drug_dict(best_similarity, smiles_list, col_list, names_list, merging_df)
    return(ret)

def compute_tanimoto_against_dataset(smis, merging_df, df, dataset_name, smi_col='SMILES', name_col = 'Name'):

    df = df[[name_col, smi_col]]
    not_nans = [type(smi) != float for smi in list(df[smi_col])]
    df = df[not_nans]

    smiles = list(df[smi_col])
    names = list(df[name_col])
    mols = [Chem.MolFromSmiles(x) for x in smiles]
    
    new_smis = []
    new_names = []
    new_mols = []
    for smi, na, mo in zip(smiles, names, mols):
        try: 
            fp = Chem.RDKFingerprint(mo)
            new_smis.append(smi)
            new_names.append(na)
            new_mols.append(fp)
        except Exception as e:
            print(e)
            continue
    col_list = ['SMILES', 'tanimoto similarity to closest ' + dataset_name, 'closest ' + dataset_name + ' smiles', 'closest ' + dataset_name + ' name']
    gen_mols = get_lowest_tanimoto_from_drug_set(smis, new_mols, new_smis, col_list, merging_df, new_names)
 
    return(gen_mols)

## Drug Likeness Filtering ##

def keep_valid_molecules(df, smiles_column):
    smis = list(df[smiles_column])
    mols = [Chem.MolFromSmiles(smi) for smi in smis]
    keep_indices = [m is not None for m in mols]
    df = df[keep_indices]
    print('length of df with valid mols: ', len(df))
    mols = [m for i,m in enumerate(mols) if keep_indices[i]]
    smis = list(df[smiles_column])
    for m, smi in zip(mols, smis):
        m.SetProp('SMILES', smi)
    return(df, mols)

def filter_for_clogp(df, mols, smiles_column, thresh = 3):
    smis = list(df[smiles_column])
    logps = [ADME(smi)._logp() for smi in smis]
    keep_indices = [logp < thresh for logp in logps]
    df = df[keep_indices]
    print('length of df with logp < ' + str(thresh) + ': ', len(df))
    mols = [m for i,m in enumerate(mols) if keep_indices[i]]
    return(df, mols)

def filter_pains(df, mols, smiles_column, thresh = 0):
    # initialize filter
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)

    def search_for_pains_or_brenk(mol):
        entry = catalog.GetMatches(mol)  # Get # matching
        return(len(entry) <= thresh)

    keep_indices = [search_for_pains_or_brenk(m) for m in mols]
    df = df[keep_indices]
    mols = [m for i,m in enumerate(mols) if keep_indices[i]]
    print('length of all preds with less than or equal to ' + str(thresh) + ' PAINS alerts: ', len(df))
    return(df, mols)


def filter_brenk(df, mols, smiles_column, thresh = 1):
    # initialize filter
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
    catalog = FilterCatalog(params)

    def search_for_pains_or_brenk(mol):
        entry = catalog.GetMatches(mol)  # Get # matching
        return(len(entry) <= thresh)

    keep_indices = [search_for_pains_or_brenk(m) for m in mols]
    df = df[keep_indices]
    mols = [m for i,m in enumerate(mols) if keep_indices[i]]
    print('length of all preds with less than or equal to ' + str(thresh) + ' Brenk alerts: ', len(df))
    return(df, mols)

def filter_for_druglikeness_egan(df, smiles_column):
    smis = list(df[smiles_column])
    druglikeness = [ADME(smi).druglikeness_egan() for smi in smis]
    df = df[druglikeness]
    print('length of df satisfying druglikeness Egan: ', len(df))
    return(df)

def filter_for_druglikeness_ghose(df, smiles_column):
    smis = list(df[smiles_column])
    druglikeness = [ADME(smi).druglikeness_ghose() for smi in smis]
    df = df[druglikeness]
    print('length of df satisfying druglikeness Ghose: ', len(df))
    return(df)

def filter_for_druglikeness_lipinski(df, smiles_column):
    smis = list(df[smiles_column])
    druglikeness = [ADME(smi).druglikeness_lipinski() for smi in smis]
    df = df[druglikeness]
    print('length of df satisfying druglikeness Lipinski: ', len(df))
    return(df)

def filter_for_druglikeness_muegge(df, smiles_column):
    smis = list(df[smiles_column])
    druglikeness = [ADME(smi).druglikeness_muegge() for smi in smis]
    df = df[druglikeness]
    print('length of df satisfying druglikeness Muegge: ', len(df))
    return(df)

## Clustering Functions ##

# code adapted from https://www.macinchem.org/reviews/clustering/clustering.php
def clusterFps(fps,num_clusters):

    tan_array = [DataStructs.BulkTanimotoSimilarity(i, fps) for i in fps]
    tan_array = np.array(tan_array)
    clusterer= AgglomerativeClustering(n_clusters = num_clusters, compute_full_tree = True).fit(tan_array)
    final_clusters = {}
    for ix, m in enumerate(clusterer.labels_):
        if m in final_clusters:
            curr_list = final_clusters[m]
            curr_list.append(ix)
            final_clusters[m] = curr_list
        else:
            final_clusters[m] = [ix]
    
    return clusterer.labels_, final_clusters

def determine_optimal_clustering_number(df, max_num_clusters):

    df['row_num'] = list(range(len(df)))
    df = df.drop_duplicates(subset='SMILES')
    smis = list(df['SMILES'])
    mols = [Chem.MolFromSmiles(mol) for mol in smis]
    murcks = [MurckoScaffold.GetScaffoldForMol(mol) for mol in mols]
    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in murcks]
    
    max_dists = []
    avg_dists = []
    print(max_num_clusters)
    for number_of_clusters in range(1, max_num_clusters):
        raw_cluster_labels, final_clusters=clusterFps(fps,num_clusters=number_of_clusters)
        max_dist = []
        avg_dist = []
        for cluster_key in final_clusters:
            cluster_mols = final_clusters[cluster_key]
            cluster_mols = [mols[i] for i in cluster_mols]
            
            # get similarities
            #cluster_murcks = [MurckoScaffold.GetScaffoldForMol(mol) for mol in cluster_mols]
            cluster_fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in cluster_mols]
            tan_array = [DataStructs.BulkTanimotoSimilarity(i, cluster_fps) for i in cluster_fps]
            flattened_tan_array = [item for sublist in tan_array for item in sublist]
            avg_dist.append(np.mean(flattened_tan_array))
            max_dist.append(np.min(flattened_tan_array))
        max_dists.append(np.average(max_dist))
        avg_dists.append(np.average(avg_dist))
    plt.scatter(list(range(1,max_num_clusters)), max_dists)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Average of Minimum Similarity Within Cluster')
    plt.show()
    plt.scatter(list(range(1,max_num_clusters)), avg_dists)
    plt.xlabel('Number of Clusters')
    plt.ylabel('Average of Mean Similarity Within Cluster')
    plt.show()

    return(df, max_dists)

def extract_legends_and_plot(df, name, folder, num_clusters, name_col):

    df['row_num'] = list(range(len(df)))
    df = df.drop_duplicates(subset='SMILES')
    smis = list(df['SMILES'])
    legends = []
    row_num = 0
    
    subImgSize= (500,500)
    mols = [Chem.MolFromSmiles(mol) for mol in smis]
    
    for smi in smis:
        try:
            mol = Chem.MolFromSmiles(smi)
            row = df.iloc[row_num,:]
            actual_row_num = str(row.loc['row_num'])
            actualrowname = str(row.loc[name_col])
            killsco = str(row.loc['hit_kill'])
            inhsco = str(row.loc['hit_inh'])
            tanabx = str(row.loc['tanimoto similarity to closest abx'])
            tants = str(row.loc['tanimoto similarity to closest train set'])
            legend = str(actual_row_num) + ', ' + actualrowname + '\n' + 'kill: ' + str(np.round(float(killsco),3)) + ', inh: ' + str(np.round(float(inhsco),3)) + '\n tanabx: ' + str(np.round(float(tanabx),3)) + '\n tan ts: ' + str(np.round(float(tants),3))
        except Exception as e:
            print(e)
            actual_row_num = str(row.loc['row_num'])
            legend = 'row: ' + actual_row_num
        mols[row_num].SetProp('legend', legend)
        legends.append(legend)
        row_num = row_num + 1 
    
    # code with help from the OG greg landrum: https://gist.github.com/greglandrum/d5f12058682f6b336905450e278d3399
    molsPerRow = 4
    subImgSize= (500,500)
    
    murcks = [MurckoScaffold.GetScaffoldForMol(mol) for mol in mols]
    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2,1024) for x in murcks]
    raw_cluster_labels, final_clusters=clusterFps(fps,num_clusters=num_clusters)
    img_list = []
    name_index = 0
    
    for cluster_key in final_clusters:
        cluster_mols = final_clusters[cluster_key]
        cluster_mols = [mols[i] for i in cluster_mols]

        nRows = len(cluster_mols) // molsPerRow
        if len(cluster_mols) % molsPerRow:
            nRows += 1
        fullSize = (molsPerRow * subImgSize[0], nRows * subImgSize[1])
        d2d = rdMolDraw2D.MolDraw2DCairo(fullSize[0],fullSize[1], subImgSize[0], subImgSize[1])
        d2d.drawOptions().legendFontSize=100
        #d2d.drawOptions().useBWAtomPalette()
        d2d.DrawMolecules(cluster_mols,legends=[mol.GetProp('legend') for mol in cluster_mols])
        d2d.FinishDrawing()
        new_name = folder + str(name_index) + '_' + name
        open(new_name,'wb+').write(d2d.GetDrawingText())
        img_list.append(new_name)
        name_index = name_index + 1
    df['cluster'] = [str(i) for i in raw_cluster_labels] # so it gets interpreted by plotly in a good way
    return(df, mols)
        
