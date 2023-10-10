#from geneformer import TranscriptomeTokenizer
import os 
os.chdir("/home/lsxgf/Geneformer/geneformer/")
from tokenizer import TranscriptomeTokenizer

# if previously wrote the 'ensembl_id', the following might not be necessary.
import scanpy as sc 
#tem = sc.read_h5ad("/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/NSCLC_subsetted_raw.h5ad")
#tem.var['ensembl_id'] = tem.var['features']
#tem.__dict__['_raw'].__dict__['_var'] = tem.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'ensembl_id'})
#tem.write_h5ad("/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/NSCLC_subsetted/NSCLC_subsetted_raw.h5ad")

input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/ms/ms_test_ori/"
output_dir = input_dir
prefix = "ms_test"

tk = TranscriptomeTokenizer(nproc=15)

tk.tokenize_data(input_dir, 
                 output_dir, 
                 prefix, 
                 ensembl_id_key="index_column",
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
                          input_dir+prefix+".dataset/",
                          output_dir,
                          prefix+"_geneformer_out")


# plot UMAP of cell embeddings
# note: scanpy umap necessarily saves figs to figures directory
#embex.plot_embs(embs=embs, 
#                plot_style="umap",
#                output_directory="path/to/output_directory/",  
#                output_prefix="emb_plot")

