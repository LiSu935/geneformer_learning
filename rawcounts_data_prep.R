# conda activate r-mofa

library(Seurat)
library(SeuratDisk)
library(dplyr)


setwd('/storage/htc/joshilab/Su_Li/Alg_development/scRNA_TCRBCR_surfaceProtein/data_collection/anti-PD1_human_Tcell_NSCLC/')
raw = readRDS("GSE179994_all.Tcell.rawCounts.rds")

# Read the first column as a list from your table (replace 'your_table.csv' with your actual file path)
first_column_list <- read.csv('GSE179994_all.scTCR.tsv', header = TRUE,sep='\t')$CellName


intersection <- intersect(first_column_list, colnames(raw))


# Assuming your_matrix is a matrix with column names
subsetted_raw <- raw[, intersection]

library("AnnotationDbi")
library("org.Hs.eg.db")

columns(org.Hs.eg.db)
# [1] "ACCNUM"       "ALIAS"        "ENSEMBL"      "ENSEMBLPROT"  "ENSEMBLTRANS"
# [6] "ENTREZID"     "ENZYME"       "EVIDENCE"     "EVIDENCEALL"  "GENENAME"    
#[11] "GENETYPE"     "GO"           "GOALL"        "IPI"          "MAP"         
#[16] "OMIM"         "ONTOLOGY"     "ONTOLOGYALL"  "PATH"         "PFAM"        
#[21] "PMID"         "PROSITE"      "REFSEQ"       "SYMBOL"       "UCSCKG"      
#[26] "UNIPROT"  

ENSEMBL_id <- mapIds(org.Hs.eg.db, keys=rownames(subsetted_raw), column="ENSEMBL", keytype="SYMBOL", multiVals="first")
head(ENSEMBL_id)

new_dataframe <- data.frame(gene_name = rownames(subsetted_raw), ensembl_id = ENSEMBL_id)

rownames(subsetted_raw) = ENSEMBL_id
#grep("", rownames(subsetted_raw))
rows_to_remove <- which(is.na(rownames(subsetted_raw)) | rownames(subsetted_raw) == "")

# Remove rows with empty or NA row names
cleaned_raw <- subsetted_raw[-rows_to_remove, ]

print(dim(cleaned_raw))

unique_raw <- cleaned_raw[!duplicated(rownames(cleaned_raw)), ]
print(dim(unique_raw))

# do the same thing for the new_dataframe:
id_rows_to_remove <- which(is.na(new_dataframe$ensembl_id) | new_dataframe$ensembl_id == "")

# Remove rows with empty or NA row names
id_dataframe <- new_dataframe[-id_rows_to_remove, ]

ID_dataframe <- id_dataframe[!duplicated(id_dataframe$ensembl_id), ]
print(dim(ID_dataframe))


unique_raw_ensID_df = data.frame(ensembl_id = rownames(unique_raw))

unique_raw_ID_df = left_join(unique_raw_ensID_df, ID_dataframe, by = "ensembl_id")
print(dim(unique_raw_ID_df))



seurat_object <- CreateSeuratObject(counts = unique_raw, 
                                    project = "NSCLC",
                                    min.cells = 1, min.features = 1)

# 17745 features across 77030 samples within 1 assay 

seurat_object$n_counts = seurat_object$nCount_RNA

meta_df=seurat_object[[]]
meta_df$cellid = rownames(meta_df)

label_df_tem = read.csv('GSE179994_Tcell.metadata.tsv', header = TRUE,sep='\t')[,c("cellid","cluster")]
colnames(label_df_tem) = c("cellid" ,"cell_type")

library(dplyr)

meta_df_tem = left_join(meta_df, label_df_tem, by='cellid')
seurat_object$cell_type = meta_df_tem$cell_type

# final gene name ID table:
seurat_object_genename_df = data.frame(ensembl_id = rownames(seurat_object))
seurat_object_genename_df$gene_name = rownames(seurat_object_genename_df) 

# change the feature name back to gene name:
rownames(seurat_object@assays$RNA@counts)<-seurat_object_genename_df$gene_name
rownames(seurat_object@assays$RNA@data)<-seurat_object_genename_df$gene_name

SaveH5Seurat(seurat_object, filename = "NSCLC_subsetted_raw.h5Seurat",overwrite = TRUE)
Convert("NSCLC_subsetted_raw.h5Seurat", dest = "h5ad", overwrite = TRUE)

# save ID mapping table for later use in scanpy
write.table(seurat_object_genename_df, file=paste0("NSCLC_subsetted_raw","_geneName_IDMapping.csv",sep=""), row.names = F, col.names = T, sep=",",quote = FALSE)

