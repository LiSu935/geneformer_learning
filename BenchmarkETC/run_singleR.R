# conda env: r-mofa-singler 
# conda env on HB: r-singlecell-r4.1-py3.10
# this is for use part of the data as reference and do the supervised learning:

suppressPackageStartupMessages({
    
  library(Seurat)
  library(SeuratDisk)
  library(SingleR)
  library(optparse)
})
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = "/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/TrainData/", 
              help = "input path", metavar = "character"),
  make_option(c("-o", "--output_dir"), type = "character", default = '/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/TestData/',
              help = "output path", metavar = "character"),
  make_option(c("-r","--ref_file"),type = "character",default = "hPancreas_train.h5ad",
              help = "ref_file", metavar = "character"),
  make_option(c("-q","--querry_file"), type = "character", default = "hPancreas_test.h5ad",
              help = "querry_file", metavar = "character"),
  make_option(c("-l","--ref_label_file"), type = "character", default = "ms_train_label.csv",
              help = "ref_label_file", metavar = "character")
)



opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/Benchmark_data_top2000/COVID/"
output_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_singleR/results/Benchmark_data_top2000/COVID/"
train_file_input = "COVID_train_df.csv"
test_file_input = "COVID_test_df.csv"
ref_label_file_input = "COVID_train_label.csv"


       
input_dir = opt$input_dir
output_dir = opt$output_dir
train_file_input = opt$ref_file
test_file_input = opt$querry_file
ref_label_file_input = opt$ref_label_file

dir.create(output_dir, recursive = TRUE)



train_singler = function(train_file="seu_obj_few25_s400.h5ad", test_file='seu_obj_ForTest.h5ad', ref_label_file) {
      
       train_pre = strsplit(train_file, split='.',fixed=T)[[1]][1]
       train_ext = strsplit(train_file, split='.',fixed=T)[[1]][2]
       test_pre = strsplit(test_file, split='.',fixed=T)[[1]][1]
       test_ext = strsplit(test_file, split='.',fixed=T)[[1]][2]  

       if (train_ext == 'h5ad') {
         Convert(paste0(input_dir, train_file), dest = "h5seurat", overwrite = TRUE)
         s400_ob <- LoadH5Seurat(paste0(input_dir, train_pre, ".h5seurat"), meta.data = FALSE, misc = FALSE)
         s400_ob_mat = s400_ob@assays$RNA@data
       } else if (train_ext == 'csv') {
         train_mat = read.csv(paste0(input_dir, train_file), row.names = 1, header= TRUE)#row.names = NULL, header= FALSE)
         s400_ob_mat = as.matrix(train_mat)
       }
       
      
       if (test_ext == 'h5ad') {
         Convert(paste0(input_dir, test_file), dest = "h5seurat", overwrite = TRUE)
         test_ob <- LoadH5Seurat(paste0(input_dir, test_pre, ".h5seurat"), meta.data = FALSE, misc = FALSE)
         test_ob_mat = test_ob@assays$RNA@data
       } else if (test_ext == 'csv') {
         test_mat = read.csv(paste0(input_dir, test_file), row.names = 1, header= TRUE)#row.names = NULL, header= FALSE)
         test_ob_mat = as.matrix(test_mat)
       }

       
       # load label for the s400_ob_mat:
       s400_ob_mat_label = read.csv(paste0(input_dir,ref_label_file), header= FALSE)
       
       # run singleR and output the labels
       
       pred.test <- SingleR(test=test_ob_mat, ref=s400_ob_mat, labels=s400_ob_mat_label$V1, de.method="wilcox", de.n=50)
       table(pred.test$labels)
       write.table(pred.test['labels'], file=paste(output_dir ,test_pre,"_", train_pre, "_singler_labels",".csv",sep=""), row.names = TRUE, col.names = NA, sep=",", quote = FALSE)
       
}

train_singler(train_file=train_file_input,test_file=test_file_input, ref_label_file = ref_label_file_input)



# output the test df from h5ad, since above loading has issue.
#import scanpy as sc
#test = sc.read_h5ad("/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/TestData/seu_obj_ForTest.h5ad")
#import pandas as pd 
#test_df = pd.DataFrame(data=test.layers['logcounts'].toarray(), index=test.obs_names, columns=test.var_names)
#test_df_t = test_df.T
#test_df_t.to_csv("/storage/htc/joshilab/Su_Li/Alg_development/scbert/data/data/TestData/test_df.csv")

# 20240120 
# for MS data:
#import pandas as pd
#import scanpy as sc

#input_dir = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/Benchmark_data_top2000/ms/"
#train = sc.read_h5ad(input_dir+"ms_train.h5ad")
#train_df = pd.DataFrame(data=train.X.toarray(), index=train.obs_names, columns=train.var_names)
#train_df_t = train_df.T
#train_df_t.to_csv(input_dir+"ms_train_df.csv")

#test = sc.read_h5ad(input_dir+"ms_test.h5ad")
#test_df = pd.DataFrame(data=test.X.toarray(), index=test.obs_names, columns=test.var_names)
#test_df_t = test_df.T
#test_df_t.to_csv(input_dir+"ms_test_df.csv")

train_file=train_file_input
test_file=test_file_input
