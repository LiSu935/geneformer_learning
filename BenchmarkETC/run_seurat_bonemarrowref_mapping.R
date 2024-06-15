# seurat reference mapping
# HB: r-singlecell-4.3-py3.10

library(Azimuth)
library(Seurat)
library(SeuratData)
#bonemarrowref <- LoadData("bonemarrowref", "azimuth")
#InstallData("bmcite")
bm <- LoadData(ds = "bmcite")

bm <- RunUMAP(bm, nn.name = "weighted.nn", reduction.name = "wnn.umap", 
              reduction.key = "wnnUMAP_", return.model = TRUE)
DimPlot(bm, group.by = "celltype.l2", reduction = "wnn.umap") 

bm <- ScaleData(bm, assay = 'RNA')
bm <- RunSPCA(bm, assay = 'RNA', graph = 'wsnn')

bm <- FindNeighbors(
  object = bm,
  reduction = "spca",
  dims = 1:50,
  graph.name = "spca.annoy.neighbors", 
  k.param = 50,
  cache.index = TRUE,
  return.neighbor = TRUE,
  l2.norm = TRUE
)

# query datasets:
input_dir = "/cluster/pixstor/xudong-lab/suli/Alg_others/scBERT_jobs_Fei/data/immune_data/"

#[1] "GSM4138872_scRNA_BMMC_D1T1.rds" "GSM4138873_scRNA_BMMC_D1T2.rds"
#[3] "GSM4138874_scRNA_CD34_D2T1.rds" "GSM4138875_scRNA_CD34_D3T1.rds"
#[5] "GSM4138876_scRNA_PBMC_D4T1.rds" "GSM4138877_scRNA_PBMC_D4T2.rds"

D1T1 = readRDS(paste0(input_dir, "GSM4138872_scRNA_BMMC_D1T1.rds"))
D1T2 = readRDS(paste0(input_dir, "GSM4138873_scRNA_BMMC_D1T2.rds"))
D2T1 = readRDS(paste0(input_dir, "GSM4138874_scRNA_CD34_D2T1.rds"))
D3T1 = readRDS(paste0(input_dir, "GSM4138875_scRNA_CD34_D3T1.rds"))

D1T1_ob <- CreateSeuratObject(counts = D1T1, project = "D1T1")
D1T2_ob <- CreateSeuratObject(counts = D1T2, project = "D1T2")
D2T1_ob <- CreateSeuratObject(counts = D2T1, project = "D2T1")
D3T1_ob <- CreateSeuratObject(counts = D3T1, project = "D3T1")

query_ob <- merge(D1T1_ob, y = c(D1T2_ob, D2T1_ob, D3T1_ob), add.cell.ids = c("D1T1","D1T2", "D2T1", "D3T1"), project = "query")
query_list <- SplitObject(query_ob, split.by = "orig.ident")

query_list <- lapply(X = query_list, FUN = NormalizeData, verbose = FALSE)


# Mapping:
# https://satijalab.org/seurat/articles/multimodal_reference_mapping.html#example-2-mapping-human-bone-marrow-cells

anchors <- list()
for (i in 1:length(query_list)) {
  anchors[[i]] <- FindTransferAnchors(
    reference = bm,
    query = query_list[[i]],
    k.filter = NA,
    reference.reduction = "spca", 
    reference.neighbors = "spca.annoy.neighbors", 
    dims = 1:50
  )
}

for (i in 1:length(query_list)) {
  query_list[[i]] <- MapQuery(
    anchorset = anchors[[i]], 
    query = query_list[[i]],
    reference = bm, 
    refdata = list(
      celltype = "celltype.l2", 
      predicted_ADT = "ADT"),
    reference.reduction = "spca",
    reduction.model = "wnn.umap"
  )
}

query_merged <- merge(query_list[[1]], query_list[2:length(query_list)], merge.dr = "ref.umap")


saveRDS(query_merged, paste0(input_dir,"query_merged_bmcite_ref_mapping.rds"))
preds_df = as.data.frame(query_merged[["predicted.celltype"]])
write.table(preds_df, file=paste(input_dir ,"query_merged_bmcite_predictions.txt",sep=""), row.names = TRUE, col.names = NA, sep="\t")


# Unsupervised integration and clustering:
library(Seurat)
library(SeuratData)

# query datasets:
input_dir = "/cluster/pixstor/xudong-lab/suli/Alg_others/scBERT_jobs_Fei/data/immune_data/"

#[1] "GSM4138872_scRNA_BMMC_D1T1.rds" "GSM4138873_scRNA_BMMC_D1T2.rds"
#[3] "GSM4138874_scRNA_CD34_D2T1.rds" "GSM4138875_scRNA_CD34_D3T1.rds"
#[5] "GSM4138876_scRNA_PBMC_D4T1.rds" "GSM4138877_scRNA_PBMC_D4T2.rds"

D1T1 = readRDS(paste0(input_dir, "GSM4138872_scRNA_BMMC_D1T1.rds"))
D1T2 = readRDS(paste0(input_dir, "GSM4138873_scRNA_BMMC_D1T2.rds"))
D2T1 = readRDS(paste0(input_dir, "GSM4138874_scRNA_CD34_D2T1.rds"))
D3T1 = readRDS(paste0(input_dir, "GSM4138875_scRNA_CD34_D3T1.rds"))

D1T1_ob <- CreateSeuratObject(counts = D1T1, project = "D1T1")
D1T1_ob[['sample']] = "D1T1"
D1T2_ob <- CreateSeuratObject(counts = D1T2, project = "D1T2")
D1T2_ob[['sample']] = "D1T2"
D2T1_ob <- CreateSeuratObject(counts = D2T1, project = "D2T1")
D2T1_ob[['sample']] = "D2T1"
D3T1_ob <- CreateSeuratObject(counts = D3T1, project = "D3T1")
D3T1_ob[['sample']] = "D3T1"

query_ob <- merge(D1T1_ob, y = c(D1T2_ob, D2T1_ob, D3T1_ob), add.cell.ids = c("D1T1","D1T2", "D2T1", "D3T1"), project = "query")

perform_inte = function(ifnb) {
  ifnb <- NormalizeData(ifnb)
  ifnb <- FindVariableFeatures(ifnb)
  ifnb <- ScaleData(ifnb)
  ifnb <- RunPCA(ifnb)
  
  ifnb <- FindNeighbors(ifnb, dims = 1:30, reduction = "pca")
  ifnb <- FindClusters(ifnb, resolution = 2, cluster.name = "unintegrated_clusters")
  
  ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  
  ifnb <- IntegrateLayers(object = ifnb, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
      verbose = FALSE)
  
  # re-join layers after integration
  ifnb[["RNA"]] <- JoinLayers(ifnb[["RNA"]])
  
  ifnb <- FindNeighbors(ifnb, reduction = "integrated.cca", dims = 1:30)
  ifnb <- FindClusters(ifnb, resolution = 1)
  
  ifnb <- RunUMAP(ifnb, dims = 1:30, reduction = "integrated.cca")

}

query_ob = perform_inte(query_ob)
saveRDS(query_ob, paste0(input_dir,"query_merged_unsupervised_clustering.rds"))
