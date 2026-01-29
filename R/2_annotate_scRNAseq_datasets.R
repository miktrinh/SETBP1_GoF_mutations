# Combine and annotate scRNAseq datasets: SGS + Wang21 + MDS

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(Seurat)
library(zeallot)
source('R/utils/sc_utils.R')

# Global variables -------------------------------------------------------------
outDir = 'Results/05_combined_scRNAseq_annotation'
sc_datasets_fp = c(sgs_fp = 'Results/03_SETBP1_annotation/2505/SETBP1_SGS_clean_noMTCells_annot_2505.RDS',
                   wang21_fp = 'Results/04_public_scRNAseq/2601/Wang21__rhoLimNone__clean_noMTCells.RDS',
                   mds_l061 = 'Results/04_public_scRNAseq/2601/L061/L061__clean_noMTCells.RDS',
                   mds_l067 = 'Results/04_public_scRNAseq/2601/L067/L067__clean_noMTCells.RDS'
)

checkmate::assert_true(all(file.exists(sc_datasets_fp)))

if(!dir.exists(outDir)){
  dir.create(outDir,recursive = T)
}

# Import datasets --------------------------------------------------------------
c(sgs_srat, wang21_srat, l061_srat, l067_srat) %<-% 
  lapply(sc_datasets_fp,function(x){readRDS(x)})

sgs_pbmc = merge_seurat_objects(sgs_srat,wang21_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
leuk_srat = merge_seurat_objects(l061_srat,l067_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
combined_srat = merge_seurat_objects(sgs_pbmc,leuk_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))

c(sgs_pbmc, leuk_srat, combined_srat) %<-% 
  lapply(c(sgs_pbmc, leuk_srat, combined_srat),function(x){standard_clustering(x)})

# Add metadata to each dataset -------------------------------------------------
sample_metadata = readxl::read_excel('Results/Supplementary Table 1.xlsx',skip = 2) %>% 
  janitor::clean_names() %>% 
  dplyr::mutate(batch = dplyr::case_when(donor_id == 'BJ9' ~ 'batch_3',
                                         donor_id == 'GOSH084' & sample_id %in% c('NB11528275','NB11528276','NB11528277','NB11528278') ~ 'batch_1',
                                         donor_id == 'GOSH084' & sample_id %in% c('Leuk13234209','Leuk13234210','Leuk13234211') ~ 'batch_2',
                                         donor_id %in% c('BJ111','BJ112','BJ113','BJ114') ~ 'batch_4',
                                         donor_id == 'L061' ~ 'batch_5',
                                         donor_id == 'L067' ~ 'batch_6',
                                         .default = 'external'
                                         ))
table(sample_metadata$donor_id,sample_metadata$batch)

# add sample-level metadata
c(sgs_pbmc, leuk_srat, combined_srat) %<-% 
  lapply(c(sgs_pbmc, leuk_srat, combined_srat),function(srat){
    sample_mdat_columns_to_keep = unique(c('sample_id',setdiff(colnames(sample_metadata),c(colnames(srat@meta.data),'number_of_cells_post_qc'))))
    
    srat$sample_id = gsub('^SB\\.','',srat$orig.ident)
    checkmate::assert_true(all(srat$cellID == rownames(srat@meta.data)))
    checkmate::assert_true(all(srat$sample_id %in% sample_metadata$sample_id))
    
    srat@meta.data = srat@meta.data %>% 
      dplyr::left_join(sample_metadata[,sample_mdat_columns_to_keep], by = 'sample_id')
    rownames(srat@meta.data) = srat@meta.data$cellID
    srat
  })

# Cell-type annotation ---------------------------------------------------------
## helper functions ------------------------------------------------------------


res <- score_clusters_markers(
  sgs_pbmc,
  pbmc_markers_short
)

res$scores          # cluster × cell-type
res$predicted       # best hit per cluster
res$ambiguity       # confidence proxy

seurat_obj$marker_celltype <- res$predicted[as.character(seurat_obj$seurat_clusters)]
seurat_obj$marker_ambiguity <- res$ambiguity[as.character(seurat_obj$seurat_clusters)]
seurat_obj <- AddModuleScore(
  seurat_obj,
  features = marker_list,
  name = "ModuleScore"
)

score_cols <- grep("^ModuleScore", colnames(seurat_obj@meta.data), value = TRUE)
names(score_cols) <- names(marker_list)

DimPlot(seurat_obj,group.by = 'seurat_clusters',label = T,label.box = T)
DimPlot(seurat_obj,group.by = 'marker_celltype',label = T,label.box = T) + NoLegend()
FeaturePlot(seurat_obj,'marker_ambiguity')

# Sub-cluster by batch to get better resolution
table(seurat_obj$batch)
seurat_list <- SplitObject(
  seurat_obj,
  split.by = "batch"
)
seurat_list <- lapply(seurat_list, standard_clustering)


for (b in names(seurat_list)) {
  seurat_obj$subcluster_batch[Cells(seurat_list[[b]])] <-
    paste0(b, "_", seurat_list[[b]]$seurat_clusters)
}