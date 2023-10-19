from geneformer import EmbExtractor

import os 
os.getcwd()

input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/cell_type_train_data/"
output_dir = input_dir
prefix = "cell_type_train_data"


embex = EmbExtractor(model_type="Pretrained",
                     num_classes=0,
                     emb_mode="cell",
                     filter_data=None,
                     max_ncells=None,
                     emb_layer=-1,
                     emb_label=["cell_type"],
                     labels_to_plot=["cell_type"],
                     forward_batch_size=200,
                     nproc=15)

# extracts embedding from input data
# example dataset: https://huggingface.co/datasets/ctheodoris/Genecorpus-30M/tree/main/example_input_files/cell_classification/disease_classification/human_dcm_hcm_nf.dataset
embs = embex.extract_embs("/mnt/pixstor/dbllab/suli/tools_related/geneformer/geneformer-12L-30M",
                          input_dir+prefix+".dataset/",
                          output_dir,
                          prefix+"_geneformer_out")



embex.plot_embs(embs=embs, 
                plot_style="umap",
                output_directory=output_dir, 
                max_ncells_to_plot=None,
                output_prefix="emb_plot",
               kwargs_dict = {"show ":False)
