# this is only for plotting the umap and calculate matrics after computing the embedding

input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/"
output_dir = input_dir
prefix = "NSCLC_subsetted"

import pandas as pd
import scanpy as sc 
from pathlib import Path
import anndata
import seaborn as sns
import glob

sc.set_figure_params(format="png")

def plot_umap(embs_df, emb_dims, label, output_file, kwargs_dict):
    only_embs_df = embs_df.iloc[:,:emb_dims]
    only_embs_df.index = pd.RangeIndex(0, only_embs_df.shape[0], name=None).astype(str)
    only_embs_df.columns = pd.RangeIndex(0, only_embs_df.shape[1], name=None).astype(str)
    vars_dict = {"embs": only_embs_df.columns}
    obs_dict = {"cell_id": list(only_embs_df.index),
                "cell_type": label}
    adata = anndata.AnnData(X=only_embs_df, obs=obs_dict, var=vars_dict)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sns.set(rc={'figure.figsize':(10,10)}, font_scale=2.3)
    sns.set_style("white")
    default_kwargs_dict = {"palette":"Set2", "size":200}
    if kwargs_dict is not None:
      default_kwargs_dict.update(kwargs_dict)
    sc.pl.umap(adata, color="cell_type", save=output_file, **default_kwargs_dict)


label = "cell_type"
output_directory = output_dir
output_prefix = prefix+"_geneformer_out"
output_prefix_label = "_" + output_prefix + f"_umap_{label}"
output_file = (Path(output_directory) / output_prefix_label).with_suffix(".pdf")
# data_directory.glob("*.{}".format(file_format))
label_list = list(sc.read_h5ad(output_directory.glob("*.{}".format("h5ad"))[0]).obs["cell_type"])
embs = pd.read_csv(output_dir+prefix+"_geneformer_out"+".csv", header=0, index_col=0)

plot_umap(embs_df=embs, emb_dims=embs.shape[1], label=label_list, output_file=output_file, kwargs_dict=None)
