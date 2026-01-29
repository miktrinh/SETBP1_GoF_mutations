##============================================================================##
##                    scRNAseq Cell Type Annotation                           ##
##       Combine and annotate scRNAseq datasets: SGS + Wang21 + MDS           ##
##============================================================================##

# Libraries --------------------------------------------------------------------
library(Matrix)
library(tidyverse)
library(Seurat)
library(ComplexHeatmap)
library(patchwork)
library(zeallot)
source('R/utils/sc_utils.R')
source("R/utils/logisticRegressionCellTypist.R")
source('R/utils/io_utils.R')

# Global parameters ------------------------------------------------------------

# Annotation Pipeline Wrapper --------------------------------------------------

#' Complete annotation pipeline for a single dataset
#' @param seurat_obj Seurat object (must have seurat_clusters)
#' @param dataset_name Name of the dataset
#' @param output_dir Output directory
#' @param marker_list Marker gene list (default: HEME_MARKERS)
#' @param marker_list_short Short marker list for visualization
#' @param run_celltypist Whether to run CellTypist
#' @param celltypist_model CellTypist model to use
#' @param params Annotation parameters
#' @return Annotated Seurat object
annotate_dataset <- function(
  seurat_obj,
  dataset_name,
  output_dir,
  marker_list = HEME_MARKERS,
  run_celltypist = TRUE,
  celltypist_model = NULL, # for default model
  params = ANNOT_PARAMS(),
  force_rerun = FALSE
) {
  
  message(sprintf("\n========================================"))
  message(sprintf("Annotating dataset: %s", dataset_name))
  message(sprintf("Number of cells: %d", ncol(seurat_obj)))
  message(sprintf("========================================\n"))
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  output_prefix <- file.path(output_dir, dataset_name)
  
  # Check if annotation already exists
  final_output <- paste0(output_prefix, "_annotated.RDS")
  if (file.exists(final_output) && !force_rerun) {
    message(sprintf("Annotation already exists: %s", final_output))
    message("Loading existing annotation. Set force_rerun=TRUE to regenerate.")
    return(readRDS(final_output))
  }
  
  # Check for required clustering
  if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    message("Running standard clustering...")
    seurat_obj <- standard_clustering(seurat_obj)
  }
  
  # 1. CellTypist Annotation ---------------------------------------------------
  if (run_celltypist) {
    ct_result_file <- paste0(output_prefix, "_celltypist_result.RDS")
    
    if (file.exists(ct_result_file)) {
      message("Loading existing CellTypist results...")
      ct_result <- readRDS(ct_result_file)
    } else {
      ct_result <- run_celltypist_annotation(
        seurat_obj,
        model = celltypist_model,
        params = params
      )
      saveRDS(ct_result, ct_result_file)
    }
    
    # Add to Seurat object
    cells_match <- match(colnames(seurat_obj), rownames(ct_result$logits))
    
    seurat_obj$celltypist_label <- ct_result$predicted_labels[cells_match]
    seurat_obj$celltypist_softmax_label <- ct_result$softmax_labels[cells_match]
    seurat_obj$celltypist_max_prob <- ct_result$max_prob[cells_match]
    seurat_obj$celltypist_entropy <- ct_result$entropy[cells_match]
    
  } else {
    ct_result <- NULL
  }
  
  # 2. Marker-based Scoring ----------------------------------------------------
  marker_result <- score_clusters_markers(
    seurat_obj,
    marker_list = marker_list,
    params = params
  )
  
  # Map cluster predictions to cells
  seurat_obj$marker_celltype <- marker_result$predicted[
    as.character(seurat_obj$seurat_clusters)
    ]
  seurat_obj$marker_ambiguity <- marker_result$ambiguity[
    as.character(seurat_obj$seurat_clusters)
    ]
  
  # Add module scores
  seurat_obj <- add_marker_module_scores(seurat_obj, marker_list)
  
  # 3. Consensus Annotation ----------------------------------------------------
  if (run_celltypist) {
    seurat_obj$consensus_celltype <- integrate_annotations(
      seurat_obj,
      celltypist_col = "celltypist_softmax_label",
      marker_col = "marker_celltype",
      ambiguity_col = "marker_ambiguity",
      ambiguity_threshold = params$ambiguity_threshold
    )
  } else {
    seurat_obj$consensus_celltype <- paste0("marker:", seurat_obj$marker_celltype)
  }
  
  # 4. Generate QC Plots -------------------------------------------------------
  if (run_celltypist) {
    plot_annotation_qc(
      seurat_obj,
      output_prefix,
      celltypist_result = ct_result,
      marker_result = marker_result
    )
  }
  
  # 5. Save Results ------------------------------------------------------------
  # Save annotated Seurat object
  saveRDS(
    seurat_obj,
    paste0(output_prefix, "_annotated.RDS")
  )
  
  # Save metadata
  write.csv(
    seurat_obj@meta.data,
    paste0(output_prefix, "_metadata_annotated.csv"),
    row.names = TRUE
  )
  
  # Save marker scoring results
  write.csv(
    marker_result$scores,
    paste0(output_prefix, "_marker_scores_by_cluster.csv"),
    row.names = TRUE
  )
  
  # Summary statistics
  summary_stats <- list(
    n_cells = ncol(seurat_obj),
    n_clusters = length(unique(seurat_obj$seurat_clusters)),
    celltypes_celltypist = if(run_celltypist) table(seurat_obj$celltypist_softmax_label) else NULL,
    celltypes_marker = table(seurat_obj$marker_celltype),
    celltypes_consensus = table(seurat_obj$consensus_celltype)
  )
  
  saveRDS(summary_stats, paste0(output_prefix, "_summary_stats.RDS"))
  
  message(sprintf("\n========================================"))
  message(sprintf("Annotation complete for: %s", dataset_name))
  message(sprintf("Output saved to: %s", output_dir))
  message(sprintf("========================================\n"))
  
  return(seurat_obj)
}

##=============================================================================##
##                    Dataset-Specific Annotation Wrappers                    ##
##=============================================================================##

#' Split SETBP1 dataset by batch and save individual RDS files
#' @param input_rds Path to combined SETBP1 Seurat object
#' @param output_dir Directory to save split objects
#' @param metadata_path Path to supplementary table with batch info
#' @return List of split Seurat objects
split_setbp1_by_batch <- function(
  input_rds = 'Results/02_SETBP1_scQC/2505/SETBP1_SGS_clean_noMTCells.RDS',
  output_dir = 'Results/02_SETBP1_scQC/2505/by_batch',
  metadata_path = 'Results/Supplementary Table 1.xlsx'
) {
  message("\n================ Splitting SETBP1 by Batch ====================\n")
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Load Seurat object
  srat <- readRDS(input_rds)
  
  # Load and add batch metadata
  sample_metadata <- readxl::read_excel(metadata_path, skip = 2) %>% 
    janitor::clean_names() %>% 
    dplyr::mutate(batch = dplyr::case_when(
      donor_id == 'BJ9' ~ 'batch_3',
      donor_id == 'GOSH084' & sample_id %in% c('NB11528275','NB11528276','NB11528277','NB11528278') ~ 'batch_1',
      donor_id == 'GOSH084' & sample_id %in% c('Leuk13234209','Leuk13234210','Leuk13234211') ~ 'batch_2',
      donor_id %in% c('BJ111','BJ112','BJ113','BJ114') ~ 'batch_4',
      .default = 'external'
    ))
  
  # Add batch info to Seurat object
  srat$sample_id <- gsub('^SB\\.', '', srat$orig.ident)
  
  if (!"batch" %in% colnames(srat@meta.data)) {
    srat@meta.data <- srat@meta.data %>%
      left_join(
        sample_metadata[, c('sample_id', 'batch')],
        by = 'sample_id'
      )
    rownames(srat@meta.data) <- srat@meta.data$cellID
  }
  
  # Check batch distribution
  message("Batch distribution:")
  print(table(srat$batch))
  
  # Split by batch
  batches <- unique(srat$batch[!is.na(srat$batch)])
  batches <- batches[batches != 'external']  # Exclude external if present
  
  srat_list <- list()
  
  for (b in batches) {
    message(sprintf("\nProcessing %s...", b))
    
    batch_srat <- subset(srat, subset = batch == b)
    batch_srat <- standard_clustering(batch_srat)
    
    message(sprintf("  Cells: %d", ncol(batch_srat)))
    message(sprintf("  Samples: %d", length(unique(batch_srat$sample_id))))
    
    # Save individual batch
    batch_file <- file.path(output_dir, paste0("SETBP1_", b, "_clean.RDS"))
    saveRDS(batch_srat, batch_file)
    message(sprintf("  Saved to: %s", batch_file))
    
    srat_list[[b]] <- batch_srat
  }
  
  message("\n================ SETBP1 Splitting Complete ====================\n")
  
  return(srat_list)
}

#' Annotate SETBP1 batches
annotate_setbp1_batches <- function(
  batch_dir = 'Results/02_SETBP1_scQC/2505/by_batch',
  output_base_dir = 'Results/03_scRNAseq_annotation',
  run_suffix = NULL,
  force_rerun = FALSE
) {
  message("\n================ Annotating SETBP1 Batches ====================\n")
  
  # Create versioned output directory
  version_dir <- create_versioned_output_dir(
    base_dir = output_base_dir,
    run_suffix = run_suffix,
    create_current_symlink = TRUE
  )
  
  # Get all batch RDS files
  batch_files <- list.files(batch_dir, pattern = "SETBP1_batch.*_clean\\.RDS$", full.names = TRUE)
  batch_names <- gsub(".*SETBP1_|_clean\\.RDS", "", basename(batch_files))
  
  if (length(batch_files) == 0) {
    stop("No batch files found. Run split_setbp1_by_batch() first.")
  }
  
  message(sprintf("Found %d batches to annotate", length(batch_files)))
  
  results <- list()
  
  for (i in seq_along(batch_files)) {
    batch_name <- batch_names[i]
    batch_file <- batch_files[i]
    
    message(sprintf("\n---- Processing %s (%d/%d) ----", 
                    batch_name, i, length(batch_files)))
    
    srat <- readRDS(batch_file)
    
    output_dir <- file.path(version_dir, batch_name)
    
    srat_annotated <- tryCatch({
      annotate_dataset(
        seurat_obj = srat,
        dataset_name = paste0("SETBP1_", batch_name),
        output_dir = output_dir,
        marker_list = HEME_MARKERS,
        marker_list_short = HEME_MARKERS,
        run_celltypist = TRUE,
        celltypist_model = NULL
      )
    }, error = function(e) {
      message(sprintf("Error annotating %s: %s", batch_name, e$message))
      NULL
    })
    
    results[[batch_name]] <- srat_annotated
  }
  
  # Log this run
  log_annotation_run(
    base_dir = output_base_dir,
    version_dir = version_dir,
    datasets = paste("SETBP1 batches:", paste(batch_names, collapse = ", ")),
    params = ANNOT_PARAMS
  )
  
  message("\n==================== SETBP1 Batch Annotation Complete ====================\n")
  
  return(results)
}

#' Annotate L061 SETBP1 leukemia
annotate_l061 <- function(
  output_base_dir = 'Results/04_public_scRNAseq/2601/L061/annotation',
  run_suffix = NULL,
  force_rerun = FALSE
) {
  message("\n==================== Annotating L061 ====================\n")
  
  # Create versioned output directory
  version_dir <- create_versioned_output_dir(
    base_dir = output_base_dir,
    run_suffix = run_suffix,
    create_current_symlink = TRUE
  )
  
  srat <- readRDS('Results/04_public_scRNAseq/2601/L061/L061__clean_noMTCells.RDS')
  
  # Add metadata
  projMan <- readxl::read_excel('~/lustre_mt22/projectManifest.xlsx', sheet = 'SETBP1')
  projMan <- projMan[!is.na(projMan$assay) & grepl('scRNA', projMan$assay) & 
                       projMan$DonorID == 'L061', ]
  
  metadata <- projMan[, c('DonorID', 'externalID', 'Tissue', 'Sex', 'assay', 'sangerSampleID')]
  colnames(metadata) <- c('donorID', 'externalID', 'tissue', 'sex', 'assay', 'sangerSampleID')
  metadata$channelID <- gsub('^SB_', '', metadata$sangerSampleID)
  
  # Merge metadata
  srat$sample_id <- gsub('^SB\\.', '', srat$orig.ident)
  srat@meta.data <- srat@meta.data %>%
    left_join(metadata, by = c("sample_id" = "channelID"))
  rownames(srat@meta.data) <- srat@meta.data$cellID
  
  srat_annotated <- annotate_dataset(
    seurat_obj = srat,
    dataset_name = "L061_SETBP1",
    output_dir = version_dir,
    marker_list = HEME_MARKERS,
    marker_list_short = HEME_MARKERS,
    run_celltypist = TRUE
  )
  
  # Log this run
  log_annotation_run(
    base_dir = output_base_dir,
    version_dir = version_dir,
    datasets = "L061",
    params = ANNOT_PARAMS()
  )
  
  return(srat_annotated)
}

#' Annotate L067 GATA1 leukemia
annotate_l067 <- function(
  output_base_dir = 'Results/04_public_scRNAseq/2601/L067/annotation',
  run_suffix = NULL,
  force_rerun = FALSE
) {
  message("\n==================== Annotating L067 ====================\n")
  
  # Create versioned output directory
  version_dir <- create_versioned_output_dir(
    base_dir = output_base_dir,
    run_suffix = run_suffix,
    create_current_symlink = TRUE
  )
  
  srat <- readRDS('Results/04_public_scRNAseq/2601/L067/L067__clean_noMTCells.RDS')
  
  # Add metadata
  projMan <- readxl::read_excel('~/lustre_mt22/projectManifest.xlsx', sheet = 'GOSH_others_included')
  projMan <- projMan[!is.na(projMan$assay) & !grepl('WGS', projMan$assay) & 
                       projMan$donorID == 'L067', ]
  
  metadata <- projMan[, c('donorID', 'Tissue', 'Sex', 'assay', 'sangerSampleID')]
  colnames(metadata) <- c('donorID', 'tissue', 'sex', 'assay', 'sangerSampleID')
  metadata$channelID <- gsub('^SB_', '', metadata$sangerSampleID)
  
  # Merge metadata
  srat$sample_id <- gsub('^SB\\.', '', srat$orig.ident)
  srat@meta.data <- srat@meta.data %>%
    left_join(metadata, by = c("sample_id" = "channelID"))
  rownames(srat@meta.data) <- srat@meta.data$cellID
  
  srat_annotated <- annotate_dataset(
    seurat_obj = srat,
    dataset_name = "L067_GATA1",
    output_dir = version_dir,
    marker_list = HEME_MARKERS,
    marker_list_short = HEME_MARKERS,
    run_celltypist = TRUE
  )
  
  # Log this run
  log_annotation_run(
    base_dir = output_base_dir,
    version_dir = version_dir,
    datasets = "L067",
    params = ANNOT_PARAMS()
  )
  
  return(srat_annotated)
}

#' Annotate Wang et al. 2021 thyroid dataset
annotate_wang2021 <- function(
  input_rds = 'Results/04_public_scRNAseq/2601/Wang21__rhoLimNone__clean_noMTCells.RDS',
  output_base_dir = 'Results/04_public_scRNAseq/2601/Wang21/annotation',
  run_suffix = NULL,
  force_rerun = FALSE
) {
  message("\n==================== Annotating Wang et al. 2021 ====================\n")
  
  # Create versioned output directory
  version_dir <- create_versioned_output_dir(
    base_dir = output_base_dir,
    run_suffix = run_suffix,
    create_current_symlink = TRUE
  )
  
  srat <- readRDS(input_rds)
  
  srat_annotated <- annotate_dataset(
    seurat_obj = srat,
    dataset_name = "Wang2021",
    output_dir = version_dir,
    marker_list = HEME_MARKERS,
    marker_list_short = HEME_MARKERS,
    run_celltypist = TRUE,
    celltypist_model = NULL
  )
  
  # Log this run
  log_annotation_run(
    base_dir = output_base_dir,
    version_dir = version_dir,
    datasets = "Wang2021 (combined)",
    params = ANNOT_PARAMS
  )
  
  return(srat_annotated)
}

##============================================================================##
##                             Main Execution                                 ##
##============================================================================##

#' Run annotation on all datasets
main_annotation <- function(
  datasets = c("setbp1_batches", "l061", "l067", "wang2021"),
  split_setbp1_first = TRUE) {
  
  results <- list()
  
  # Step 0: Split SETBP1 if needed
  if ("setbp1_batches" %in% datasets && split_setbp1_first) {
    message("\n=============== Step 0: Splitting SETBP1 ====================\n")
    results$setbp1_split <- tryCatch(
      split_setbp1_by_batch(),
      error = function(e) {
        message("Error splitting SETBP1: ", e$message)
        NULL
      }
    )
  }
  
  # Step 1: Annotate SETBP1 batches
  if ("setbp1_batches" %in% datasets) {
    results$setbp1_batches <- tryCatch(
      annotate_setbp1_batches(),
      error = function(e) {
        message("Error annotating SETBP1 batches: ", e$message)
        NULL
      }
    )
  }
  
  # Step 2: Annotate L061
  if ("l061" %in% datasets) {
    results$l061 <- tryCatch(
      annotate_l061(),
      error = function(e) {
        message("Error annotating L061: ", e$message)
        NULL
      }
    )
  }
  
  # Step 3: Annotate L067
  if ("l067" %in% datasets) {
    results$l067 <- tryCatch(
      annotate_l067(),
      error = function(e) {
        message("Error annotating L067: ", e$message)
        NULL
      }
    )
  }
  
  # Step 4: Annotate Wang2021
  if ("wang2021" %in% datasets) {
    results$wang2021 <- tryCatch(
      annotate_wang2021(),
      error = function(e) {
        message("Error annotating Wang2021: ", e$message)
        NULL
      }
    )
  }
  
  message("\n==================== All Annotations Complete ====================\n")
  return(results)
}


results <- main_annotation()


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
i=1
DimPlot(seurat_list[[i]],group.by = c('seurat_clusters'),
        label = T,label.box = T,label.size = 2) + 
  NoLegend()+ggtitle(names(seurat_list)[i])

table(seurat_list[[i]]$marker_celltype[seurat_list[[i]]$seurat_clusters == 8])

library(patchwork)
out_pdf <- "Results/05_combined_scRNAseq_annotation/sgs_pbmc_batch_subcluster_qc.pdf"
# Estimate height dynamically (helps with many marker sets)
dotplot_height <- max(4, length(pbmc_markers_short) * 0.25)
pdf(out_pdf, width = 13, height = 6 + dotplot_height)

for (i in seq_along(seurat_list)) {
  
  obj <- seurat_list[[i]]
  batch_name <- names(seurat_list)[i]
  
  ## ---- DimPlots ----
  p_clusters <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "seurat_clusters",
    repel = TRUE,
    label = TRUE,
    label.box = TRUE,
    label.size = 3
  ) +
    ggtitle("Seurat clusters") +
    theme_bw(base_size = 12)+
    NoLegend() 
  
  p_celltypes <- DimPlot(
    obj,
    reduction = "umap",
    group.by = "marker_celltype",
    repel = TRUE,
    label = TRUE,
    label.box = TRUE,
    label.size = 3
  ) +
    ggtitle("Marker cell types") +
    theme_bw(base_size = 12)+
    NoLegend() 
  
  top_row <- p_clusters + p_celltypes + plot_layout(ncol = 2)
  
  ## ---- DotPlot ----
  Idents(obj) <- "seurat_clusters"
  
  p_dot <- DotPlot(
    obj,
    features = unique(unlist(pbmc_markers_short, use.names = FALSE)),
    dot.scale = 6
  ) +
    RotatedAxis() +
    scale_color_viridis_c(option = "C") +
    ggtitle("Marker expression by cluster") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ## ---- Combine ----
  combined <- top_row / p_dot +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(
      title = paste("Batch:", batch_name),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
      )
    )
  
  print(combined)
}

dev.off()


for (b in names(seurat_list)) {
  seurat_obj$subcluster_batch[Cells(seurat_list[[b]])] <-
    paste0(b, "_", seurat_list[[b]]$seurat_clusters)
}