#input_dir = "/home/lsxgf/tem/NSCLC_subsetted_741/"
#output_dir = input_dir
#prefix = "NSCLC_subsetted_741"

import getopt
import argparse

parser = argparse.ArgumentParser(description='run_geneformer with raw data, Plotting the umap and compute the matrics.')

parser.add_argument('--input_dir', type=str, default='/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/', help='Directory of input: default(%(default)s).')
parser.add_argument('--prefix', type=str, default='NSCLC_subsetted', help=' prefix of output files: default(%(default)s).')

args = parser.parse_args()

input_dir = args.input_dir
prefix = args.prefix
output_dir = input_dir


from geneformer import TranscriptomeTokenizer

# if previously wrote the 'ensembl_id', the following might not be necessary.
import scanpy as sc 
#tem = sc.read_h5ad("/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/NSCLC_subsetted_raw.h5ad")
# read the ensembl_id from the saved id mapping file:
# import pandas as pd
# id_df = pd.read_csv("NSCLC_subsetted_raw_geneName_IDMapping.csv", header=0)
#tem.var['ensembl_id'] = id_df['ensembl_id']
# tem.var['ensembl_id'] = list(id_df['ensembl_id'])
#tem.__dict__['_raw'].__dict__['_var'] = tem.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})
#tem.write_h5ad("/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/NSCLC_subsetted_raw.h5ad", compression='gzip')

#tem = sc.read_h5ad("NSCLC_subsetted_raw.h5ad")
#import pandas as pd
#id_df = pd.read_csv("NSCLC_subsetted_raw_geneName_IDMapping.csv", header=0)
#tem.var['ensembl_id'] = list(id_df['ensembl_id'])
#tem.var[['_index', 'ensembl_id', 'features']].head()
#tem.var_names = tem.var['features']
#del tem.var['_index']
#tem.var[[ 'ensembl_id', 'features']].head()
#tem.write_h5ad("NSCLC_subsetted_raw.h5ad", compression='gzip')




#input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/"
#output_dir = input_dir
#prefix = "NSCLC_subsetted"

tk = TranscriptomeTokenizer(nproc=15)

tk.tokenize_data(input_dir, 
                 output_dir, 
                 prefix, 
                 file_format="h5ad")

import sys
sys.path.append('/home/lsxgf/Geneformer/geneformer/')
from emb_extractor import EmbExtractor
# initiate EmbExtractor
embex = EmbExtractor(model_type="Pretrained",
                     num_classes=0,
                     emb_mode="cell",
                     filter_data=None,
                     max_ncells=None,
                     emb_layer=-1,
                     emb_label=None,#["cell_type"],
                     labels_to_plot=None, #["cell_type"],
                     forward_batch_size=200,
                     nproc=15)

# extracts embedding from input data
# example dataset: https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/tree/main/example_input_files/cell_classification/disease_classification/human_dcm_hcm_nf.dataset
embs = embex.extract_embs("/mnt/pixstor/dbllab/suli/tools_related/geneformer/geneformer-12L-30M",
                          input_dir+prefix+".dataset/",
                          output_dir,
                          prefix+"_geneformer_out")


# plot UMAP of cell embeddings
# note: scanpy umap necessarily saves figs to figures directory
#embex.plot_embs(embs=embs, 
#                plot_style="umap",
#                output_directory="path/to/output_directory/",  
#                output_prefix="emb_plot")

import pandas as pd
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
    with plt.rc_context({"figure.figsize": (8, 8), "figure.dpi": (300)}):  # Use this to set figure params like size and dpi
        sc.pl.umap(adata, color="celltype", **default_kwargs_dict, show=False)
        plt.savefig(output_file, bbox_inches="tight")
    return(adata)


label = "celltype"
output_directory = output_dir
output_prefix = prefix+"_geneformer_out"
output_prefix_label = "_" + output_prefix + f"_umap_{label}"
output_file = (Path(output_directory) / output_prefix_label).with_suffix(".png")
# data_directory.glob("*.{}".format(file_format))
tem = sc.read_h5ad(glob.glob(output_directory+"*.{}".format("h5ad"))[0])
if "celltype" in tem.obs.keys():
    label = "celltype"
elif "cell_type" in tem.obs.keys():
    label = "cell_type"
else:
    print("cell type info is not in the original h5ad file!! Please check!!!")
label_list = list(tem.obs[label])

#embs = pd.read_csv(glob.glob(output_dir+prefix+"_geneformer_out"+"*.csv")[0], header=0, index_col=0)

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
