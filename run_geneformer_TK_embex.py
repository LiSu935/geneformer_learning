# Hellbender
# source activate py-geneformer 

import scanpy as sc 
import os

# if previously wrote the 'ensembl_id', the following might not be necessary.
def rename_ad_geneName(h5ad_fileName="NSCLC_subsetted_raw.h5ad", oldID="index_column"):
  adata = sc.read_h5ad(input_dir+h5ad_fileName)
  adata.var['ensembl_id'] = adata.var[oldID]
  #adata.__dict__['_raw'].__dict__['_var'] = adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'ensembl_id'})
  adata.write_h5ad(output_dir+h5ad_fileName)


input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/ms/"
output_dir = input_dir+"rename_ad/"

if not os.path.exists(output_dir):
   os.makedirs(output_dir)
  

for file in os.listdir(input_dir):
  if file.endswith(".h5ad"):
    rename_ad_geneName(h5ad_fileName=file, oldID="index_column")
  

from geneformer import TranscriptomeTokenizer
tk = TranscriptomeTokenizer(nproc=15)


tk.tokenize_data("/mnt/pixstor/dbllab/suli/Alg_development/use_scgpt/", 
                 "/mnt/pixstor/dbllab/suli/Alg_development/use_scgpt/", 
                 "NSCLC_subsetted", 
                 file_format="h5ad")


from geneformer import EmbExtractor
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
                          "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/NSCLC_subsetted.dataset/",
                          "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/",
                          "NSCLC_subsetted_geneformer_out")


#File "/mnt/pixstor/data/lsxgf/miniconda/envs/py-geneformer/lib/python3.8/site-packages/huggingface_hub/utils/_validators.py", line 164, in validate_repo_id
#    raise HFValidationError(
#huggingface_hub.utils._validators.HFValidationError: Repo id must use alphanumeric chars or '-', '_', '.', '--' and '..' are forbidden, '-' and '.' cannot start or end the name, max length is 96: '../geneformer-12L-30M'.


# what pretrained model to use? model_type?
# Feel necessary to check the size of "only_embs_df.index" https://huggingface.co/ctheodoris/Geneformer/blob/main/geneformer/emb_extractor.py line 149
