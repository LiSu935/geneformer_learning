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



input_dir = "/home/lsxgf/tem/NSCLC_subsetted_741/"
output_dir = input_dir
prefix = "NSCLC_subsetted_741"


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

from pathlib import Path
import anndata
import seaborn as sns

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


import pandas as pd
import scanpy as sc 

label = "cell_type"
output_directory = output_dir
output_prefix = prefix+"_geneformer_out"
output_prefix_label = "_" + output_prefix + f"_umap_{label}"
output_file = (Path(output_directory) / output_prefix_label).with_suffix(".pdf")
label_list = list(sc.read_h5ad("/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/NSCLC_subsetted_raw.h5ad").obs["cell_type"])
embs = pd.read_csv(output_dir+prefix+"_geneformer_out"+".csv", header=0, index_col=0)

plot_umap(embs_df=embs, emb_dims=embs.shape[1], label=label_list, output_file=output_file, kwargs_dict=None)


