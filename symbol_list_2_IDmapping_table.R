# on HB, conda activate r-singlecell-r4.1-py3.10
# to use: 
# SCRIPT='symbol_list_2_IDmapping_table.R'
# Rscript ${SCRIPT} -i PATH_OF_input_dir

library(Seurat)
library(SeuratDisk)
library(tools)
library("AnnotationDbi")
library("org.Hs.eg.db")
library(dplyr)

library(optparse)
option_list <- list(
  make_option(c("-i", "--input_dir"), type = "character", default = "/group/xulab/Su_Li/spaVG/dataset/Simulation_20220110/Scenario_1/", 
              help = "input path", metavar = "character")
 
)
opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# input parameter
input_dir <- opt$input_dir

gen_id_mapping_df = function(ORI_INPUT_DIR){
  for (gene_name_File_INPUT in list.files(path = ORI_INPUT_DIR, pattern = "genesymbol_list\\.txt$")) {
    # Get the file name without extension
    prefix <- tools::file_path_sans_ext(basename(gene_name_File_INPUT))
    
    data <- read.table(paste0(ORI_INPUT_DIR, gene_name_File_INPUT), header = FALSE, sep = "\t")
    
    # Rename the column in the dataframe
    colnames(data) <- c("symbol")
    
    
    ENSEMBL_id <- mapIds(org.Hs.eg.db, keys=data$symbol, column="ENSEMBL", keytype="SYMBOL", multiVals="first")
    head(ENSEMBL_id)
    
    new_dataframe <- data.frame(symbol = data$symbol, ensembl_id = ENSEMBL_id)
    #rownames(new_dataframe) = new_dataframe$ensembl_id
    rows_to_remove <- which(is.na(new_dataframe$ensembl_id) | new_dataframe$ensembl_id == "")
    
    # Remove rows with empty or NA row names
    if (length(rows_to_remove) != 0) {
      subset_dataframe <- new_dataframe[-rows_to_remove, ]
    } else {
      subset_dataframe <- new_dataframe
    }
        
    # since multiple gene name refers to same ensembl id, have to 
    df_unique <- subset_dataframe[!duplicated(subset_dataframe$ensembl_id), ] 
    
    write.table(df_unique, file=paste0(ORI_INPUT_DIR ,prefix,"_ensembl.txt"), row.names = FALSE, col.names = TRUE, sep="\t")
  }
}

#ORI_INPUT_DIR = "/mnt/pixstor/dbllab/suli/Alg_development/use_geneformer/data/Benchmark_data_top2000/COVID/"
gen_id_mapping_df(ORI_INPUT_DIR=input_dir)

