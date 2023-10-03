import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from sklearn.decomposition import PCA
from sklearn.manifold import TSNE

# function to run PCA and t-SNE on the principal components
def pca_tsne_mols(fp_list, fp_labels, colors_for_vis, num_components = 2048):
    # PCA first, using 2 components
    pca = PCA(n_components=2)
    crds = pca.fit_transform(fp_list)

    # print variance explained by first two principal components
    print('variance explained by pc1+pc2: ' + str(np.sum(pca.explained_variance_ratio_)))

    # get the information for both PC1 and PC2
    crds_df = pd.DataFrame(crds,columns=["PC_1","PC_2"])
    crds_df['label'] = fp_labels
    crds_df.head()

    # plot the PCA done on the first two principal components
    plt.figure(figsize=(5,5), dpi = 300)
    ax = sns.scatterplot(data=crds_df,x="PC_1",y="PC_2",hue="label", alpha = 0.7, s = 10, palette=colors_for_vis)
    
    # TSNE next - first do PCA and calculate all principal components (not just first two)
    pca = PCA(n_components=num_components)
    crds = pca.fit_transform(fp_list)

    # compute t-SNE on principal components
    crds_embedded = TSNE(n_components=2).fit_transform(crds)

    # get the information for t-SNE 1 and t-SNE 2
    tsne_df = pd.DataFrame(crds_embedded,columns=["X","Y"])
    tsne_df['label'] = fp_labels

    # plot the t-SNE
    plt.figure(figsize=(5,5), dpi = 300)
    ax = sns.scatterplot(data=tsne_df,x="X",y="Y",hue="label", alpha = 0.7,  s = 10, palette=colors_for_vis)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.show()
    
    # return the t-SNE df
    return(tsne_df)