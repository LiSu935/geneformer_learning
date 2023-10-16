# this is only for plotting the umap and calculate matrics after computing the embedding

#input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/"
#output_dir = input_dir
#prefix = "NSCLC_subsetted"

#input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/ms/ms_test_ori/"
#output_dir = input_dir
#prefix = "ms_test"

import getopt
import argparse

parser = argparse.ArgumentParser(description='Plotting the umap and compute the matrics.')

parser.add_argument('--input_dir', type=str, default='/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/', help='Directory of input: default(%(default)s).')
parser.add_argument('--prefix', type=str, default='NSCLC_subsetted', help=' prefix of output files: default(%(default)s).')

args = parser.parse_args()

input_dir = args.input_dir
prefix = args.prefix
output_dir = input_dir




import pandas as pd
import scanpy as sc 
from pathlib import Path
import anndata
import seaborn as sns
import glob

sc.set_figure_params(format="png")

from matplotlib import pyplot as plt

def plot_umap(embs_df, emb_dims, label, output_file, kwargs_dict):
    only_embs_df = embs_df.iloc[:,:emb_dims]
    only_embs_df.index = pd.RangeIndex(0, only_embs_df.shape[0], name=None).astype(str)
    only_embs_df.columns = pd.RangeIndex(0, only_embs_df.shape[1], name=None).astype(str)
    vars_dict = {"embs": only_embs_df.columns}
    obs_dict = {"cell_id": list(only_embs_df.index),
                "celltype": label}
    adata = anndata.AnnData(X=only_embs_df, obs=obs_dict, var=vars_dict)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    #sc.tl.louvain(adata)
    sns.set(rc={'figure.figsize':(10,10)}, font_scale=2.3)
    sns.set_style("white")
    default_kwargs_dict = {"palette":"Set2", "size":200}
    if kwargs_dict is not None:
        default_kwargs_dict.update(kwargs_dict)
    with plt.rc_context():  # Use this to set figure params like size and dpi
        sc.pl.umap(adata, color="celltype", save=output_file, **default_kwargs_dict, show=False)
        plt.savefig(output_file, bbox_inches="tight")
    return(adata)


    

label = "celltype"
output_directory = output_dir
output_prefix = prefix+"_geneformer_out"
output_prefix_label = "_" + output_prefix + f"_umap_{label}"
output_file = (Path(output_directory) / output_prefix_label).with_suffix(".png")
# data_directory.glob("*.{}".format(file_format))
label_list = list(sc.read_h5ad(glob.glob(output_directory+"*.{}".format("h5ad"))[0]).obs["celltype"])
embs = pd.read_csv(glob.glob(output_dir+prefix+"_geneformer_out"+"*.csv")[0], header=0, index_col=0)

adata = plot_umap(embs_df=embs, emb_dims=embs.shape[1], label=label_list, output_file=output_file, kwargs_dict=None)

# the following need to check carefully.
from sklearn.metrics import silhouette_score,silhouette_samples,davies_bouldin_score,calinski_harabasz_score

sih_score = silhouette_score(X=adata.obsp['distances'], labels=label_list, metric="precomputed")
print(f"silhouette_score is {sih_score}")
print(sih_score)

sih_sample = silhouette_samples(X=adata.obsp['distances'], labels=label_list, metric="precomputed")
print(len(sih_sample))

db_score = davies_bouldin_score(X=adata.obsp['distances'].toarray(), labels=label_list)
print(f"davies_bouldin_score is {db_score}")

ch_score = calinski_harabasz_score(X=adata.obsp['distances'].toarray(), labels=label_list)
print(f"calinski_harabasz_score is {ch_score}")
