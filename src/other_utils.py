import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

def pca_tsne_mols(fp_list, fp_labels, colors_for_vis, num_components = 2048):
    # PCA first
    pca = PCA(n_components=2)
    crds = pca.fit_transform(fp_list)

    print('variance explained by pc1+pc2: ' + str(np.sum(pca.explained_variance_ratio_)))

    crds_df = pd.DataFrame(crds,columns=["PC_1","PC_2"])
    crds_df['label'] = fp_labels
    crds_df.head()

    plt.figure(figsize=(5,5), dpi = 300)
    ax = sns.scatterplot(data=crds_df,x="PC_1",y="PC_2",hue="label", alpha = 0.7, s = 10, palette=colors_for_vis)
    
    # TSNE next
    pca = PCA(n_components=num_components)
    crds = pca.fit_transform(fp_list)

    crds_embedded = TSNE(n_components=2).fit_transform(crds)

    tsne_df = pd.DataFrame(crds_embedded,columns=["X","Y"])
    tsne_df['label'] = fp_labels

    plt.figure(figsize=(5,5), dpi = 300)
    ax = sns.scatterplot(data=tsne_df,x="X",y="Y",hue="label", alpha = 0.7,  s = 10, palette=colors_for_vis)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.show()
    
    return(tsne_df)