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
  #seurat_obj <- add_marker_module_scores(seurat_obj, marker_list)
  
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
      marker_result = marker_result,
      params=params
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
    
    output_dir <- file.path(version_dir, "SETBP1_inhouse", batch_name)
    
    srat_annotated <- tryCatch({
      annotate_dataset(
        seurat_obj = srat,
        dataset_name = paste0("SETBP1_", batch_name),
        output_dir = output_dir,
        marker_list = HEME_MARKERS,
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
    params = ANNOT_PARAMS()
  )
  
  message("\n==================== SETBP1 Batch Annotation Complete ====================\n")
  
  return(results)
}

#' Annotate L061 SETBP1 leukemia
annotate_l061 <- function(
  output_base_dir = 'Results/03_scRNAseq_annotation',
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
    output_dir = file.path(version_dir,"L061_SETBP1"),
    marker_list = HEME_MARKERS,
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
  output_base_dir = 'Results/03_scRNAseq_annotation',
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
    output_dir = file.path(version_dir,"L067_GATA1"),
    marker_list = HEME_MARKERS,
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
  output_base_dir = 'Results/03_scRNAseq_annotation',
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
    output_dir = file.path(version_dir,"Wang2021"),
    marker_list = HEME_MARKERS,
    run_celltypist = TRUE,
    celltypist_model = NULL
  )
  
  # Log this run
  log_annotation_run(
    base_dir = output_base_dir,
    version_dir = version_dir,
    datasets = "Wang2021 (combined)",
    params = ANNOT_PARAMS()
  )
  
  return(srat_annotated)
}

##============================================================================##
##                             Main Execution                               ####
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


results <- main_annotation(
  datasets = c("l067", "wang2021"),
  split_setbp1_first = FALSE
)


# Manually finalising annotations ----------------------------------------------

CELLTYPE_HIERARCHY <- tibble::tribble(
  ~celltype_lvl3, ~celltype_lvl2, ~lineage,

  # B lineage
  "naive_b_cells", "Naive B", "B",
  "b_cells", "Naive B", "B",
  "memory_b_cells", "Memory B", "B",
  "age_associated_b_cells", "Memory B", "B",
  "plasma_cells", "Plasma B", "B",
  "plasmablasts", "Plasma B", "B",

  # T lineage
  "cd8a_a", "Naive / Central Memory CD8 T", "T",
  "tcm_naive_cytotoxic_t_cells", "Naive / Central Memory CD8 T", "T",
  "cytotoxic_t_cells", "CD8 T","T",
  "helper_t_cells", "CD4 T","T",
  "tcm_naive_helper_t_cells", "Naive / Central Memory CD4 T", "T",

  "tem_temra_cytotoxic_t_cells", "Effector / Memory CD8 T", "T",
  "tem_trm_cytotoxic_t_cells", "Effector / Memory CD8 T", "T",

  "tem_effector_helper_t_cells", "Effector / Memory CD4 T", "T",
  "type_17_helper_t_cells", "Effector / Memory CD4 T", "T",

  "regulatory_t_cells", "Regulatory CD4 T", "T",

  "mait_cells", "MAIT cells", "T",
  "gamma-delta T cells", "gd T cells", "T",
  "gamma_delta_t_cells", "gd T cells", "T",
  "crtam_gamma_delta_t_cells", "gd T cells", "T",
  "nkt_cells","NKT","T",
  "T_CD8",'t_cd8','T',
  "T_CD4",'t_cd4','T',

  # NK / ILC
  "nk_cells", "NK", "NK/ILC",
  "cd16_nk_cells", "NK", "NK/ILC",
  "ilc3", "ILC", "NK/ILC",
  
  # Myeloid
  "classical_monocytes", "classical_monocytes", "Myeloid",
  "non_classical_monocytes", "non_classical_monocytes", "Myeloid",
  "macrophages","macrophages","Myeloid",

  # Dendritic cells
  "dc_precursor", "DC_precursor", "DC",
  "dc1", "DC1", "DC",
  "dc2", "DC2", "DC",
  "dc3", "DC3", "DC",
  "pdc", "pDC", "DC",

  # Other
  "mast_cells", "Mast", "Mast",
  "megakaryocytes_platelets", "Megakaryocyte", "Megakaryocyte",
  "platelets","Megakaryocyte", "Megakaryocyte",
  "unknown", "Unknown", "Unknown",
  "doublets", "doublets","doublets",
  "cycling_doublets", "cycling_doublets","doublets",
  "MDS",'MDS','MDS',
  "cycling_cells","cycling_cells","cycling_cells"
) %>%
  dplyr::mutate(
    celltype_lvl2_snakecase = celltype_lvl2 %>%
      str_to_lower() %>%
      str_replace_all("[^a-z0-9]+", "_") %>%
      str_replace_all("^_|_$", "")
  )

finalised_annot = list()


batch_1 = readRDS('Results/03_scRNAseq_annotation/current/SETBP1_inhouse/batch_1/SETBP1_batch_1_annotated.RDS')
batch_2 = readRDS('Results/03_scRNAseq_annotation/current/SETBP1_inhouse/batch_2/SETBP1_batch_2_annotated.RDS')
batch_3 = readRDS('Results/03_scRNAseq_annotation/current/SETBP1_inhouse/batch_3/SETBP1_batch_3_annotated.RDS')
batch_4 = readRDS('Results/03_scRNAseq_annotation/current/SETBP1_inhouse/batch_4/SETBP1_batch_4_annotated.RDS')

## SETBP1-inhouse - batch 1 ----------------------------------------------------
srat = batch_1
srat$annot = annot$annot_2601[match(srat$cellID,annot$cellID)]
DimPlot(srat, group.by = 'seurat_clusters',label = T,repel = T,label.box = T) + NoLegend()
DimPlot(srat,cells.highlight = srat$cellID[srat$annot == 'cycling_t_cells'])
FeaturePlot(srat, 'marker_ambiguity')

annot = srat@meta.data %>%
  dplyr::select(c(cellID,seurat_clusters,celltypist_label,celltypist_softmax_label,marker_celltype)) %>%
  dplyr::mutate(
    celltypist_label_og = celltypist_label,
    celltypist_softmax_label_og = celltypist_softmax_label,
    across(
      c(celltypist_label, celltypist_softmax_label, marker_celltype),
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )
annot = annot %>%
  dplyr::mutate(annot_2601 = dplyr::case_when(
    
    celltypist_label == celltypist_softmax_label ~ celltypist_label,
    celltypist_label_og %in% c("Double-positive thymocytes","Proliferative germinal center B cells") ~ "doublets",
    celltypist_label_og %in% c("Myelocytes") ~ "classical_monocytes",
    celltypist_label_og %in% c("Tem/Effector helper T cells PD1+") ~ "helper_t_cells",
    
    celltypist_label_og %in% c("B cells","Transitional B cells") ~ "memory_b_cells",
    celltypist_label_og %in% c("Age-associated B cells",
                               "CD16- NK cells","CD16+ NK cells",'CD8a/a',
                               "Classical monocytes","CRTAM+ gamma-delta T cells",
                               "DC2","ILC3","MAIT cells",
                               "Mast cells","Megakaryocytes/platelets",
                               "Memory B cells","NK cells","NKT cells","Non-classical monocytes",
                               "Plasma cells",
                               "Regulatory T cells","Tem/Effector helper T cells",
                               "Tem/Temra cytotoxic T cells","Tcm/Naive cytotoxic T cells","Tem/Trm cytotoxic T cells",
                               "Tcm/Naive helper T cells"
    ) ~ celltypist_label,
    .default = 'unknown'
  ))
table(annot$annot_2601)
table(annot$celltypist_label_og[annot$annot_2601=='unknown'],
      annot$seurat_clusters[annot$annot_2601=='unknown'])


df = as.data.frame(table(annot$celltypist_label[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$marker_celltype[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$celltypist_label_og[annot$celltypist_label == annot$celltypist_softmax_label]))
df = df[df$Freq >0,]

df = as.data.frame(table(annot$celltypist_label_og[annot$celltypist_label != annot$celltypist_softmax_label],
                         annot$celltypist_softmax_label_og[annot$celltypist_label != annot$celltypist_softmax_label]))

DimPlot(srat,cells.highlight =
          srat$cellID[
            srat$celltypist_label == "Proliferative germinal center B cells" & 
              srat$celltypist_softmax_label == 'unknown'])
FeaturePlot(srat,HEME_MARKERS$Plasma_B)






table(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3)
unique(annot$annot_2601[!annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3])
checkmate::assert_true(all(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3))
annot <- annot %>%
  mutate(celltype_lvl3 = annot_2601) %>%
  left_join(CELLTYPE_HIERARCHY, by = "celltype_lvl3")

# Add back to Seurat object
batch_1@meta.data = cbind(batch_1@meta.data,annot[match(batch_1$cellID,annot$cellID),c('annot_2601','celltype_lvl3',
                                                                                       'celltype_lvl2', 'lineage', 'celltype_lvl2_snakecase')])
# Add to final annotation
finalised_annot[['SETBP1_inhouse_batch_1']] = annot

write.csv(annot,'Results/03_scRNAseq_annotation/setbp1_inhouse_batch1_finalised_annotation.csv')

## SETBP1-inhouse - batch 2 ----------------------------------------------------
srat = batch_2
srat$annot = annot$annot_2601[match(srat$cellID,annot$cellID)]
DimPlot(srat, group.by = 'annot',label = T,repel = T,label.box = T) + NoLegend()
DimPlot(srat,cells.highlight = srat$cellID[srat$annot == 'cycling_t_cells'])
FeaturePlot(srat, 'marker_ambiguity')

annot = srat@meta.data %>%
  dplyr::select(c(cellID,seurat_clusters,celltypist_label,celltypist_softmax_label,marker_celltype)) %>%
  dplyr::mutate(
    celltypist_label_og = celltypist_label,
    celltypist_softmax_label_og = celltypist_softmax_label,
    across(
      c(celltypist_label, celltypist_softmax_label, marker_celltype),
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )
annot = annot %>%
  dplyr::mutate(annot_2601 = dplyr::case_when(
    seurat_clusters %in% c(2) ~ 'gamma-delta T cells',
    celltypist_label == celltypist_softmax_label ~ celltypist_label,
    celltypist_label_og %in% c("Age-associated B cells","B cells","Transitional B cells") ~ "b_cells",
    celltypist_label_og %in% c("CD16- NK cells",'CD8a/a',"Classical monocytes","CRTAM+ gamma-delta T cells",
                               "DC2","ILC3","Megakaryocytes/platelets",
                               "Memory B cells","Naive B cells","Plasma cells",
                               "Regulatory T cells","Tem/Effector helper T cells",
                               "Tem/Temra cytotoxic T cells","Tcm/Naive cytotoxic T cells",
                               "CD16+ NK cells"
                               ) ~ celltypist_label,
    celltypist_label_og %in% c("CD16+ NK cells","NK cells") & seurat_clusters == 15 ~ celltypist_label,
    celltypist_label_og %in% c("Erythrophagocytic macrophages","Monocytes","Myelocytes") ~ "classical_monocytes",
    celltypist_label_og %in% c("Non-classical monocytes") & celltypist_softmax_label_og == 'unknown' ~ "classical_monocytes",
    celltypist_label_og %in% c("Non-classical monocytes") & celltypist_softmax_label_og != 'unknown' ~ "non_classical_monocytes",
    celltypist_label_og %in% c("Proliferative germinal center B cells") & seurat_clusters == 17 ~ "memory_b_cells",
    celltypist_label_og %in% c("Tcm/Naive cytotoxic T cells","Tem/Trm cytotoxic T cells") & seurat_clusters == 12 ~ celltypist_label,
    celltypist_label_og %in% c("Tcm/Naive helper T cells","Type 17 helper T cells") & seurat_clusters != 12 ~ celltypist_label,
    seurat_clusters == 12 ~ 'cytotoxic_t_cells',
    celltypist_label_og %in% c("Tem/Trm cytotoxic T cells","MAIT cells") & seurat_clusters %in% c(4,9) ~ 'helper_t_cells',
    celltypist_label_og %in% c("MAIT cells") & seurat_clusters == 10 ~ 'mait_cells',
    celltypist_label_og %in% c("MAIT cells") & seurat_clusters == 15 ~ 'nk_cells',
    celltypist_label_og %in% c("Mast cells") & seurat_clusters %in% c(14) ~ 'classical_monocytes',
    seurat_clusters == 6 ~ 'megakaryocytes_platelets',
    celltypist_label_og %in% c("T(agonist)","Cycling T cells") & seurat_clusters %in% c(17,4) ~ 'helper_t_cells',
    .default = 'unknown'
  ))
table(annot$annot_2601)
table(annot$celltypist_label_og[annot$annot_2601=='unknown'],
      annot$seurat_clusters[annot$annot_2601=='unknown'])


df = as.data.frame(table(annot$celltypist_label[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$marker_celltype[annot$celltypist_label == annot$celltypist_softmax_label]))
df = df[df$Freq >0,]

df = as.data.frame(table(annot$celltypist_label_og[annot$celltypist_label != annot$celltypist_softmax_label],
                         annot$celltypist_softmax_label_og[annot$celltypist_label != annot$celltypist_softmax_label]))


table(srat$celltypist_label[srat$seurat_clusters == 3],
      srat$celltypist_softmax_label[srat$seurat_clusters == 3],
      srat$marker_celltype[srat$seurat_clusters == 3])
DimPlot(srat,cells.highlight =
          srat$cellID[
            srat$celltypist_label == 'MAIT cells' &
              srat$annot == 'unknown'])
FeaturePlot(srat,HEME_MARKERS$Plasma_B)






table(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3)
unique(annot$annot_2601[!annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3])
checkmate::assert_true(all(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3))
annot <- annot %>%
  mutate(celltype_lvl3 = annot_2601) %>%
  left_join(CELLTYPE_HIERARCHY, by = "celltype_lvl3")

# Add back to Seurat object
batch_2@meta.data = cbind(batch_2@meta.data,annot[match(batch_2$cellID,annot$cellID),c('annot_2601','celltype_lvl3',
                                                     'celltype_lvl2', 'lineage', 'celltype_lvl2_snakecase')])
# Add to final annotation
finalised_annot[['SETBP1_inhouse_batch_2']] = annot

write.csv(annot,'Results/03_scRNAseq_annotation/setbp1_inhouse_batch2_finalised_annotation.csv')

## SETBP1-inhouse - batch 3 ----------------------------------------------------
srat = batch_3
srat$annot = annot$annot_2601[match(srat$cellID,annot$cellID)]
DimPlot(srat, group.by = 'annot',label = T,repel = T,label.box = T) + NoLegend()
DimPlot(srat,cells.highlight = srat$cellID[srat$seurat_clusters == 14])
FeaturePlot(srat, 'marker_ambiguity')

annot = srat@meta.data %>%
  dplyr::select(c(cellID,seurat_clusters,celltypist_label,celltypist_softmax_label,marker_celltype)) %>%
  dplyr::mutate(
    celltypist_label_og = celltypist_label,
    celltypist_softmax_label_og = celltypist_softmax_label,
    across(
      c(celltypist_label, celltypist_softmax_label, marker_celltype),
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )
annot = annot %>%
  dplyr::mutate(annot_2601 = dplyr::case_when(
    seurat_clusters == 14 ~ 'memory_b_cells',
    celltypist_label == 'tcm_naive_helper_t_cells' & marker_celltype == 'classical_mono' & !seurat_clusters %in% c(17,19,10,4,21,13) ~ 'doublets',
    celltypist_label %in% c('tcm_naive_helper_t_cells','tcm_naive_cytotoxic_t_cells','naive_b_cells','memory_b_cells','cd16_nk_cells') & marker_celltype == 'platelets' ~ 'doublets',
    celltypist_label_og %in% c("Classical monocytes") & seurat_clusters %in% c(9,12,11) ~ "doublets",
    celltypist_label_og %in% c("Double-negative thymocytes","Double-positive thymocytes","DC",
                               "Epithelial cells","Erythrophagocytic macrophages","Follicular helper T cells",
                               "Plasma cells","Plasmablasts","Proliferative germinal center B cells") ~ "doublets",
    
    celltypist_label == celltypist_softmax_label ~ celltypist_label,
    celltypist_label_og %in% c("Age-associated B cells") ~ "memory_b_cells",
    celltypist_label_og %in% c("Small pre-B cells","Transitional B cells","B cells") ~ "naive_b_cells",
    celltypist_label_og %in% c("T(agonist)") ~ "cytotoxic_t_cells",
    celltypist_label_og %in% c("Alveolar macrophages","Cycling monocytes","Intermediate macrophages","Monocytes","Myelocytes") ~ "classical_monocytes",
    seurat_clusters == 10 ~ 'classical_monocytes',
    celltypist_label_og %in% c("Naive B cells") & seurat_clusters %in% c(5,9,10,12,15,20) ~ "doublets",
    celltypist_label_og %in% c("CD16- NK cells","CD16+ NK cells",
                               'CD8a/a',"Classical monocytes",
                               "CRTAM+ gamma-delta T cells","DC precursor","DC1",
                               "DC2","gamma-delta T cells",
                               "ILC3","Macrophages","MAIT cells","Megakaryocytes/platelets",
                               "Memory B cells","Naive B cells","NK cells","NKT cells",
                               "Non-classical monocytes","pDC",
                               "Regulatory T cells",
                               "Tcm/Naive cytotoxic T cells","Tem/Temra cytotoxic T cells","Tem/Trm cytotoxic T cells",
                               "Tcm/Naive helper T cells","Tem/Effector helper T cells"
    ) ~ celltypist_label,
    celltypist_label_og %in% c("Type 1 helper T cells","Type 17 helper T cells") ~ 'helper_t_cells',
    .default = 'unknown'
  ))
table(annot$annot_2601)
table(annot$celltypist_label_og[annot$annot_2601=='doublets'],
      annot$seurat_clusters[annot$annot_2601=='doublets'])


df = as.data.frame(table(annot$celltypist_label[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$marker_celltype[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$celltypist_label_og[annot$celltypist_label == annot$celltypist_softmax_label]))
df = df[df$Freq >0,]

df = as.data.frame(table(annot$celltypist_label_og[annot$celltypist_label != annot$celltypist_softmax_label],
                         annot$celltypist_softmax_label_og[annot$celltypist_label != annot$celltypist_softmax_label]))

DimPlot(srat,cells.highlight =
          srat$cellID[
            srat$celltypist_label == 'B cells' &
              srat$seurat_clusters == 16 ])

FeaturePlot(srat,HEME_MARKERS$Effector_CD8_T)


table(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3)
unique(annot$annot_2601[!annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3])
checkmate::assert_true(all(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3))
annot <- annot %>%
  mutate(celltype_lvl3 = annot_2601) %>%
  left_join(CELLTYPE_HIERARCHY, by = "celltype_lvl3")

# Add back to Seurat object
batch_3@meta.data = cbind(batch_3@meta.data,annot[match(batch_3$cellID,annot$cellID),c('annot_2601','celltype_lvl3',
                                                                                       'celltype_lvl2', 'lineage', 'celltype_lvl2_snakecase')])
# Add to final annotation
finalised_annot[['SETBP1_inhouse_batch_3']] = annot

write.csv(annot,'Results/03_scRNAseq_annotation/setbp1_inhouse_batch3_finalised_annotation.csv')

## SETBP1-inhouse - batch 4 ----------------------------------------------------
srat = batch_4
srat$annot = annot$annot_2601[match(srat$cellID,annot$cellID)]
DimPlot(srat, group.by = 'seurat_clusters',label = T,repel = T,label.box = T) + NoLegend()
DimPlot(srat,cells.highlight = srat$cellID[srat$seurat_clusters == 14])
FeaturePlot(srat, 'marker_ambiguity')

annot = srat@meta.data %>%
  dplyr::select(c(cellID,seurat_clusters,celltypist_label,celltypist_softmax_label,marker_celltype)) %>%
  dplyr::mutate(
    celltypist_label_og = celltypist_label,
    celltypist_softmax_label_og = celltypist_softmax_label,
    across(
      c(celltypist_label, celltypist_softmax_label, marker_celltype),
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )
annot = annot %>%
  dplyr::mutate(annot_2601 = dplyr::case_when(
    celltypist_label == 'classical_monocytes' & marker_celltype == 'platelets' ~ 'doublets',
    celltypist_label %in% c('memory_b_cells','naive_b_cells') & marker_celltype == 'memory_cd4_t' ~ 'doublets',
    celltypist_label %in% c('tcm_naive_cytotoxic_t_cells','tcm_naive_helper_t_cells') & marker_celltype == 'classical_mono' ~  'doublets',
    celltypist_label %in% c('tcm_naive_cytotoxic_t_cells','tcm_naive_helper_t_cells','tem_trm_cytotoxic_t_cells') & marker_celltype == 'platelets' ~ 'platelets',
    celltypist_label_og %in% c("Migratory DCs") ~ "dc2",
    celltypist_label_og %in% c("Trm cytotoxic T cells") ~ 'cytotoxic_t_cells',
    celltypist_label_og %in% c("B cells","Transitional B cells", "Pro-B cells") ~ "naive_b_cells",
    celltypist_label_og %in% c("Age-associated B cells") ~ "memory_b_cells",
    celltypist_label_og %in% c("Alveolar macrophages") ~ "classical_monocytes",
    seurat_clusters == 13 ~ 'gamma_delta_t_cells',
    seurat_clusters == 38 ~ 'plasma_cells',
    seurat_clusters %in% c(26,30) ~ 'doublets',
    
    celltypist_label == celltypist_softmax_label ~ celltypist_label,

    celltypist_label_og %in% c('CD16+ NK cells') & seurat_clusters %in% c(28,29,30,37,41) ~ 'doublets',
    celltypist_label_og %in% c('CD16- NK cells') & seurat_clusters %in% c(5,7,10,13,26,28,30) ~ 'doublets',
    celltypist_label_og %in% c('Cycling NK cells',"Double-negative thymocytes","Double-positive thymocytes",
                               "Early lymphoid/T lymphoid","Follicular helper T cells",
                               "Intermediate macrophages","Kupffer cells",
                               "Mid erythroid","Monocytes","Myelocytes",
                               "Plasmablasts","Pro-B cells","Proliferative germinal center B cells") ~ "doublets",
    celltypist_label_og %in% c('Memory B cells') & seurat_clusters %in% c(5,7,10,16,17,22,29,33,38) ~ 'doublets',
    celltypist_label_og %in% c('Naive B cells') & seurat_clusters %in% c(3,7,10,13,28,29) ~ 'doublets',
    celltypist_label_og %in% c('NK cells') & seurat_clusters %in% c(7,8,10,13,22,27,28,40) ~ 'doublets',
    celltypist_label_og %in% c('Plasma cells') & seurat_clusters %in% c(7,29) ~ 'doublets',
    celltypist_label_og %in% c('Treg(diff)') ~ 'regulatory_t_cells',
    celltypist_label_og %in% c('Transitional NK') ~ 'nk_cells',

    celltypist_label_og %in% c('CD8a/a',"Classical monocytes",
                               "CRTAM+ gamma-delta T cells","DC precursor","DC1",
                               "DC2","gamma-delta T cells",
                               "ILC3","Macrophages","MAIT cells","Megakaryocytes/platelets",
                               "Non-classical monocytes","pDC",
                               "Regulatory T cells",
                               "Tcm/Naive cytotoxic T cells","Tem/Temra cytotoxic T cells","Tem/Trm cytotoxic T cells",
                               "Tcm/Naive helper T cells","Tem/Effector helper T cells") ~ celltypist_label,
    
    celltypist_label_og %in% c("Type 1 helper T cells","Type 17 helper T cells") ~ 'helper_t_cells',
    .default = 'unknown'
  ))
table(annot$annot_2601)
table(annot$celltypist_label_og[annot$annot_2601=='doublets'],
      annot$seurat_clusters[annot$annot_2601=='doublets'])


df = as.data.frame(table(annot$celltypist_label[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$marker_celltype[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$celltypist_label_og[annot$celltypist_label == annot$celltypist_softmax_label]))
df = df[df$Freq >0,]

df = as.data.frame(table(annot$celltypist_label_og[annot$celltypist_label != annot$celltypist_softmax_label],
                         annot$celltypist_softmax_label_og[annot$celltypist_label != annot$celltypist_softmax_label]))

DimPlot(srat,cells.highlight =
          srat$cellID[
            srat$celltypist_label == 'Treg(diff)' &
              #srat$seurat_clusters == 40
              srat$celltypist_softmax_label == 'unknown'
              ])

table(annot$seurat_clusters[annot$celltypist_label_og=='Transitional B cells' & 
                              annot$celltypist_label!=annot$celltypist_softmax_label])
FeaturePlot(srat,HEME_MARKERS$GammaDelta_T)


table(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3)
unique(annot$annot_2601[!annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3])
checkmate::assert_true(all(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3))
annot <- annot %>%
  mutate(celltype_lvl3 = annot_2601) %>%
  left_join(CELLTYPE_HIERARCHY, by = "celltype_lvl3")

# Add back to Seurat object
batch_4@meta.data = cbind(batch_4@meta.data,annot[match(batch_4$cellID,annot$cellID),c('annot_2601','celltype_lvl3',
                                                                                       'celltype_lvl2', 'lineage', 'celltype_lvl2_snakecase')])
# Add to final annotation
finalised_annot[['SETBP1_inhouse_batch_4']] = annot

write.csv(annot,'Results/03_scRNAseq_annotation/setbp1_inhouse_batch4_finalised_annotation.csv')


## L061 ------------------------------------------------------------------------
l061 = readRDS('Results/03_scRNAseq_annotation/current/L061_SETBP1/L061_SETBP1_annotated.RDS')
srat = l061
srat$annot = annot$annot_2601[match(srat$cellID,annot$cellID)]
DimPlot(srat, group.by = 'annot',label = T,repel = T,label.box = T) + NoLegend()
DimPlot(srat,cells.highlight = srat$cellID[srat$annot == 'non_classical_monocytes'])
FeaturePlot(srat, 'marker_ambiguity')

annot = srat@meta.data %>%
  dplyr::select(c(cellID,seurat_clusters,celltypist_label,celltypist_softmax_label,marker_celltype)) %>%
  dplyr::mutate(
    celltypist_label_og = celltypist_label,
    celltypist_softmax_label_og = celltypist_softmax_label,
    across(
      c(celltypist_label, celltypist_softmax_label, marker_celltype),
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )
annot = annot %>%
  dplyr::mutate(annot_2601 = dplyr::case_when(
    seurat_clusters %in% c(7,1) ~ 'MDS',
    celltypist_label_og %in% c("Monocyte precursor","Myelocytes","Monocytes") ~ 'classical_monocytes',
    celltypist_label_og %in% c("Alveolar macrophages","Intermediate macrophages",
                               "Intestinal macrophages","Macrophages") ~ 'macrophages',
    celltypist_label_og %in% c("DC") ~ 'dc1',
    celltypist_label == celltypist_softmax_label ~ celltypist_label,
    
    celltypist_label_og %in% c("Age-associated B cells","B cells","CD16- NK cells","CD16+ NK cells") ~ 'doublets',
    celltypist_label_og %in% c("Memory B cells","Naive B cells") & seurat_clusters %in% c(2,3,4,6,8) & seurat_clusters!=10 ~ 'doublets',
    celltypist_label_og %in% c("NK cells","Regulatory T cells","Tcm/Naive cytotoxic T cells",
                               "Tcm/Naive helper T cells","Tem/Effector helper T cells",
                               "Tem/Temra cytotoxic T cells","Tem/Trm cytotoxic T cells",
                               "Treg(diff)","Trm cytotoxic T cells") & 
      seurat_clusters %in% c(10,4,6,9,12,13) ~ 'doublets',
    
    celltypist_label_og %in% c("Classical monocytes","CRTAM+ gamma-delta T cells",
                               "DC","DC Precursor","DC1","DC2","DC3","MAIT cells",
                               "Megakaryocytes/platelets","Memory B cells","Monocytes",
                               "Naive B cells","Non-classical monocytes","pDC",
                               "NK cells","Regulatory T cells","Tcm/Naive cytotoxic T cells",
                               "Tcm/Naive helper T cells","Tem/Effector helper T cells",
                               "Tem/Temra cytotoxic T cells","Tem/Trm cytotoxic T cells",
                               "Treg(diff)","DC precursor"
                               ) ~ celltypist_label,
    
    celltypist_label_og %in% c("Double-negative thymocytes","Double-positive thymocytes",
                               "Endothelial cells","Epithelial cells","Follicular helper T cells",
                               "ILC3","Large pre-B cells","Mast cells",
                               "Megakaryocyte-erythroid-mast cell progenitor","Plasma cells",
                               "T(agonist)","Fibroblasts","Mid erythroid"
                               ) ~ 'doublets',
    celltypist_label_og %in% c("Type 1 helper T cells","Type 17 helper T cells") ~ 'helper_t_cells',
    celltypist_label_og %in% c("Trm cytotoxic T cells") ~ 'cytotoxic_t_cells',
    .default = 'unknown'
  ))


table(annot$annot_2601)
table(annot$celltypist_label_og[annot$annot_2601=='unknown'],
      annot$seurat_clusters[annot$annot_2601=='unknown'])


df = as.data.frame(table(annot$celltypist_label[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$marker_celltype[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$celltypist_label_og[annot$celltypist_label == annot$celltypist_softmax_label]))
df = df[df$Freq >0,]

df = as.data.frame(table(annot$celltypist_label_og[annot$celltypist_label != annot$celltypist_softmax_label],
                         annot$celltypist_softmax_label_og[annot$celltypist_label != annot$celltypist_softmax_label]))

DimPlot(srat,cells.highlight =
          srat$cellID[
            srat$celltypist_label == 'NK cells' &
              #srat$seurat_clusters == 8
              srat$celltypist_softmax_label == 'unknown'
            ])

DotPlot(srat,features =  
          c(# Stemness / HOX
            "CD34","SOX4","MEIS1","HOXA9","HOXA10","HLF","LMO2",
            # Proliferation
            "MKI67","TOP2A","HMGB2","TYMS","CENPF",
            # Myeloid immaturity / inflammation
            "MPO","ELANE","AZU1","LGALS3",
            "S100A8","S100A9","IL1B","CXCL8",
            # Stress / transformation
            "HSP90AA1","DDIT3","ATF3")
)+RotatedAxis()


table(annot$seurat_clusters[annot$celltypist_label_og=='NK cells' & 
                              annot$celltypist_label!=annot$celltypist_softmax_label])
FeaturePlot(srat,HEME_MARKERS$GammaDelta_T)


table(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3)
unique(annot$annot_2601[!annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3])
checkmate::assert_true(all(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3))
annot <- annot %>%
  mutate(celltype_lvl3 = annot_2601) %>%
  left_join(CELLTYPE_HIERARCHY, by = "celltype_lvl3")

# Add back to Seurat object
l061@meta.data = cbind(l061@meta.data,annot[match(l061$cellID,annot$cellID),c('annot_2601','celltype_lvl3',
                                                                                       'celltype_lvl2', 'lineage', 'celltype_lvl2_snakecase')])
# Add to final annotation
finalised_annot[['L061_SETBP1']] = annot

write.csv(annot,'Results/03_scRNAseq_annotation/L061_SETBP1_finalised_annotation.csv')

## L067 ------------------------------------------------------------------------
l067 = readRDS('Results/03_scRNAseq_annotation/current/L067_GATA1/L067_GATA1_annotated.RDS')
srat = l067
srat$annot = annot$annot_2601[match(srat$cellID,annot$cellID)]
DimPlot(srat, group.by = 'seurat_clusters',label = T,repel = T,label.box = T) + NoLegend()
DimPlot(srat,cells.highlight = srat$cellID[srat$annot == 'unknown' & srat$celltypist_label == 'MEMP'])
FeaturePlot(srat, 'marker_ambiguity')

annot = srat@meta.data %>%
  dplyr::select(c(cellID,seurat_clusters,celltypist_label,celltypist_softmax_label,marker_celltype)) %>%
  dplyr::mutate(
    celltypist_label_og = celltypist_label,
    celltypist_softmax_label_og = celltypist_softmax_label,
    across(
      c(celltypist_label, celltypist_softmax_label, marker_celltype),
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )
annot = annot %>%
  dplyr::mutate(annot_2601 = dplyr::case_when(
    seurat_clusters %in% c(8) ~ 'MDS',
    celltypist_label_og %in% c("Early erythroid","MEMP") ~ 'MDS',
    celltypist_label_og %in% c("Pro-B cells") & marker_celltype == 'Memory_CD4_T' ~ 'doublets',
    celltypist_label_og %in% c("Age-associated B cells") ~ 'memory_b_cells',
    celltypist_label_og %in% c("Alveolar macrophages","Intermediate macrophages",
                               "Intestinal macrophages","Macrophages") ~ 'macrophages',
    celltypist_label_og %in% c("Monocytes") ~ 'classical_monocytes',
    celltypist_label_og %in% c("Pro-B cells","Small pre-B cells") ~ 'naive_b_cells',
    
    celltypist_label == celltypist_softmax_label ~ celltypist_label,
    celltypist_label_og %in% c("B cells","DC precursor","Double-negative thymocytes",
                               "Large pre-B cells","Macrophages","Plasma cells",
                               "Proliferative germinal center B cells") ~ 'doublets',
    
    celltypist_label_og %in% c("CD16- NK cells","CD16+ NK cells",
                               "CD8a/a","Classical monocytes",
                               "CRTAM+ gamma-delta T cells",
                               "DC","DC Precursor","DC1","DC2","DC3","ILC3",
                               "MAIT cells",
                               "Megakaryocytes/platelets","Memory B cells",
                               "Naive B cells","Non-classical monocytes","pDC",
                               "NK cells","Regulatory T cells","Tcm/Naive cytotoxic T cells",
                               "Tcm/Naive helper T cells","Tem/Effector helper T cells",
                               "Tem/Temra cytotoxic T cells","Tem/Trm cytotoxic T cells"
    ) ~ celltypist_label,
    celltypist_label_og %in% c("Type 1 helper T cells","Type 17 helper T cells") ~ 'helper_t_cells',
    celltypist_label_og %in% c("Trm cytotoxic T cells") ~ 'cytotoxic_t_cells',
    .default = 'unknown'
  ))


table(annot$annot_2601)
table(annot$celltypist_label_og[annot$annot_2601=='unknown'],
      annot$seurat_clusters[annot$annot_2601=='unknown'])


df = as.data.frame(table(annot$celltypist_label[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$marker_celltype[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$celltypist_label_og[annot$celltypist_label == annot$celltypist_softmax_label]))
df = df[df$Freq >0,]

df = as.data.frame(table(annot$celltypist_label_og[annot$celltypist_label != annot$celltypist_softmax_label],
                         annot$celltypist_softmax_label_og[annot$celltypist_label != annot$celltypist_softmax_label]))

DimPlot(srat,cells.highlight =
          srat$cellID[
            srat$celltypist_label == 'Tcm/Naive helper T cells' &
              #srat$seurat_clusters == 8
              srat$celltypist_softmax_label == 'unknown'
            ])

DotPlot(srat,features =  
          c(# Stemness / HOX
            "CD34","SOX4","MEIS1","HOXA9","HOXA10","HLF","LMO2",
            # Proliferation
            "MKI67","TOP2A","HMGB2","TYMS","CENPF",
            # Myeloid immaturity / inflammation
            "MPO","ELANE","AZU1","LGALS3",
            "S100A8","S100A9","IL1B","CXCL8",
            # Stress / transformation
            "HSP90AA1","DDIT3","ATF3")
)+RotatedAxis()


table(annot$seurat_clusters[annot$celltypist_label_og=='NK cells' & 
                              annot$celltypist_label!=annot$celltypist_softmax_label])
FeaturePlot(srat,HEME_MARKERS$Naive_CD8_T)


table(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3)
unique(annot$annot_2601[!annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3])
checkmate::assert_true(all(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3))
annot <- annot %>%
  mutate(celltype_lvl3 = annot_2601) %>%
  left_join(CELLTYPE_HIERARCHY, by = "celltype_lvl3")

# Add back to Seurat object
l067@meta.data = cbind(l067@meta.data,annot[match(l067$cellID,annot$cellID),c('annot_2601','celltype_lvl3',
                                                                              'celltype_lvl2', 'lineage', 'celltype_lvl2_snakecase')])
# Add to final annotation
finalised_annot[['L067_GATA1']] = annot

write.csv(annot,'Results/03_scRNAseq_annotation/L067_GATA1_finalised_annotation.csv')

## Wang21 ----------------------------------------------------------------------
wang21 = readRDS('Results/03_scRNAseq_annotation/current/Wang2021/Wang2021_annotated.RDS')

srat = wang21
srat$annot = annot$celltype_lvl2[match(srat$cellID,annot$cellID)]
DimPlot(srat, group.by = 'annot',label = T,repel = T,label.box = T) + NoLegend()
DimPlot(srat,cells.highlight = srat$cellID[srat$annot == 'unknown' & srat$celltypist_label == 'MEMP'])
FeaturePlot(srat, 'marker_ambiguity')

annot = srat@meta.data %>%
  dplyr::select(c(cellID,seurat_clusters,celltypist_label,celltypist_softmax_label,marker_celltype)) %>%
  dplyr::mutate(
    celltypist_label_og = celltypist_label,
    celltypist_softmax_label_og = celltypist_softmax_label,
    across(
      c(celltypist_label, celltypist_softmax_label, marker_celltype),
      ~ .x %>%
        str_to_lower() %>%
        str_replace_all("[^a-z0-9]+", "_") %>%
        str_replace_all("^_|_$", "")
    )
  )
annot = annot %>%
  dplyr::mutate(annot_2601 = dplyr::case_when(
    celltypist_label_og %in% c("Regulatory T cells","Tem/Trm cytotoxic T cells") & marker_celltype == 'Classical_Mono' ~ 'classical_monocytes',
    celltypist_label_og %in% c("Regulatory T cells") & marker_celltype == 'Platelets' ~ 'doublets',
    
    celltypist_label_og %in% c("Age-associated B cells", "Germinal center B cells") ~ 'memory_b_cells',
    celltypist_label_og %in% c("B cells") ~ 'naive_b_cells',
    
    celltypist_label_og %in% c("Alveolar macrophages","Intermediate macrophages",
                               "Intestinal macrophages","Macrophages") ~ 'macrophages',
    
    celltypist_label_og %in% c("Monocytes","Myelocytes") ~ 'classical_monocytes',
    celltypist_label_og %in% c("Pro-B cells","Small pre-B cells","Transitional B cells") ~ 'naive_b_cells',
    seurat_clusters == 10 ~ 'unknown',
    seurat_clusters == 18 ~ 'cycling_cells',
    celltypist_label == celltypist_softmax_label ~ celltypist_label,
    celltypist_label_og %in% c("B cells","DC precursor","Double-negative thymocytes",
                               "Large pre-B cells","Macrophages","Plasma cells",
                               "Proliferative germinal center B cells","CD8a/b(entry)",
                               "Double-positive thymocytes","Early lymphoid/T lymphoid",
                               "Erythrophagocytic macrophages","Mast cells",
                               "Proliferative germinal center B cells",
                               "gamma-delta T cells") ~ 'doublets',
  
    celltypist_label_og %in% c("CD16- NK cells","CD16+ NK cells") & 
      seurat_clusters %in% c(4,6,9,10,13) ~ 'doublets',
    celltypist_label_og %in% c("Memory B cells","Naive B cells") & seurat_clusters %in% c(0,1,5,10,24) ~ 'doublets',
    celltypist_label_og %in% c("Regulatory T cells") & seurat_clusters %in% c(5,28) ~ 'doublets',
    
    celltypist_label_og %in% c("CD16- NK cells","CD16+ NK cells",
                               "CD8a/a","Classical monocytes",
                               "CRTAM+ gamma-delta T cells",
                               "DC2","ILC3","MAIT cells","Megakaryocytes/platelets",
                               "Memory B cells","Naive B cells","NK cells","NKT cells",
                               "Non-classical monocytes","Plasma cells","Plasmablasts","pDC",
                               "Regulatory T cells","Tcm/Naive cytotoxic T cells",
                               "Tcm/Naive helper T cells","Tem/Effector helper T cells",
                               "Tem/Temra cytotoxic T cells","Tem/Trm cytotoxic T cells"
    ) ~ celltypist_label,
    celltypist_label_og %in% c("Follicular helper T cells","Treg(diff)",
                               "Type 1 helper T cells","Type 17 helper T cells") ~ 'helper_t_cells',
    celltypist_label_og %in% c("Trm cytotoxic T cells","Memory CD4+ cytotoxic T cells") ~ 'cytotoxic_t_cells',
    .default = 'unknown'
  ))


table(annot$annot_2601)
table(annot$seurat_clusters[annot$annot_2601=='unknown'])
table(annot$celltypist_label_og[annot$annot_2601=='unknown' & annot$seurat_clusters != 10],
      annot$seurat_clusters[annot$annot_2601=='unknown' & annot$seurat_clusters != 10])


df = as.data.frame(table(annot$celltypist_label[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$marker_celltype[annot$celltypist_label == annot$celltypist_softmax_label],
                         annot$celltypist_label_og[annot$celltypist_label == annot$celltypist_softmax_label]))
df = df[df$Freq >0,]

df = as.data.frame(table(annot$celltypist_label_og[annot$celltypist_label != annot$celltypist_softmax_label],
                         annot$celltypist_softmax_label_og[annot$celltypist_label != annot$celltypist_softmax_label]))

DimPlot(srat,cells.highlight =
          srat$cellID[
            srat$celltypist_label == "Treg(diff)" &
              #srat$seurat_clusters == 28
              #srat$marker_celltype == 'Memory_CD4_T'
              srat$celltypist_softmax_label == 'unknown'
            ])

table(annot$seurat_clusters[annot$celltypist_label_og=='Regulatory T cells' & 
                              annot$celltypist_label!=annot$celltypist_softmax_label])
FeaturePlot(srat,HEME_MARKERS$GammaDelta_T)


table(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3)
unique(annot$annot_2601[!annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3])
checkmate::assert_true(all(annot$annot_2601 %in% CELLTYPE_HIERARCHY$celltype_lvl3))
annot <- annot %>%
  mutate(celltype_lvl3 = annot_2601) %>%
  left_join(CELLTYPE_HIERARCHY, by = "celltype_lvl3")

# Add back to Seurat object
wang21@meta.data = cbind(wang21@meta.data,annot[match(wang21$cellID,annot$cellID),
                                                c('annot_2601','celltype_lvl3','celltype_lvl2', 'lineage', 'celltype_lvl2_snakecase')])
# Add to final annotation
finalised_annot[['Wang21']] = annot

write.csv(annot,'Results/03_scRNAseq_annotation/Wang21_finalised_annotation.csv')

final_annot_df = do.call(rbind,finalised_annot)

# Combine data -----------------------------------------------------------------
sgs_inhouse = readRDS('Results/03_scRNAseq_annotation/2505/SETBP1_SGS_clean_noMTCells_annot_2505.RDS')

## Add metadata to each dataset ------------------------------------------------
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
c(sgs_inhouse_annot, l061_annot, l067_annot, wang21_annot) %<-%
  lapply(c(sgs_inhouse, l061, l067, wang21),function(srat){
    sample_mdat_columns_to_keep = unique(c('sample_id',setdiff(colnames(sample_metadata),c(colnames(srat@meta.data),'number_of_cells_post_qc'))))
    
    srat$sample_id = gsub('^SB\\.','',srat$orig.ident)
    checkmate::assert_true(all(srat$cellID == rownames(srat@meta.data)))
    checkmate::assert_true(all(srat$sample_id %in% sample_metadata$sample_id))
    
    srat@meta.data = srat@meta.data %>%
      dplyr::left_join(sample_metadata[,sample_mdat_columns_to_keep], by = 'sample_id')
    rownames(srat@meta.data) = srat@meta.data$cellID
    srat
  })

# Add annotation for sgs_inhouse
sgs_inhouse_annot@meta.data = cbind(sgs_inhouse_annot@meta.data,
                                    final_annot_df[match(sgs_inhouse_annot$cellID,final_annot_df$cellID),
                                                   setdiff(colnames(final_annot_df),colnames(sgs_inhouse_annot@meta.data))])

## Combine data
sgs_pbmc = merge_seurat_objects(sgs_inhouse_annot,wang21_annot,keepAllGenes = F,genomeVersions = c('v38','v38'))
leuk_srat = merge_seurat_objects(l061_annot,l067_annot,keepAllGenes = F,genomeVersions = c('v38','v38'))
combined_srat = merge_seurat_objects(sgs_pbmc,leuk_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))

c(sgs_pbmc, leuk_srat, combined_srat) %<-%
  lapply(c(sgs_pbmc, leuk_srat, combined_srat),function(x){standard_clustering(x)})

c(sgs_pbmc, leuk_srat, combined_srat) %<-%
  lapply(c(sgs_pbmc, leuk_srat, combined_srat),function(srat){
    srat@meta.data = cbind(srat@meta.data,
                               final_annot_df[match(srat$cellID,final_annot_df$cellID),
                                              setdiff(colnames(final_annot_df),colnames(srat@meta.data))])
    srat
  })

## Finalise annotation
### SGS_PBMC ----
source('~/lustre_mt22/generalScripts/utils/misc.R')
table(sgs_pbmc$annot)

DimPlot(sgs_pbmc,group.by = 'annot',label = T,repel = T,
        label.box = T,cols = col25)+NoLegend()

DimPlot(sgs_pbmc,cells.highlight = 
          sgs_pbmc$cellID[sgs_pbmc$annot == 'unknown'])
FeaturePlot(sgs_pbmc,c('CD3D','CD8A','CD4'))
DimPlot(sgs_pbmc,cells.highlight = 
          sgs_pbmc$cellID[sgs_pbmc$seurat_clusters == 40])

table(sgs_pbmc$seurat_clusters[sgs_pbmc$annot=='unknown'])
table(sgs_pbmc$celltype_lvl2_snakecase[sgs_pbmc$seurat_clusters==32],
      sgs_pbmc$donor_id[sgs_pbmc$seurat_clusters==32])


sgs_pbmc$annot = as.character(sgs_pbmc$celltype_lvl2_snakecase)
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'dc_precursor'] = 'dc2'
sgs_pbmc$annot[sgs_pbmc$seurat_clusters == 47] = 'dc1'
sgs_pbmc$annot[sgs_pbmc$seurat_clusters == 47] = 'dc1'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase %in% c('macrophages','mast')] = 'doublets'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'megakaryocyte' & 
                 sgs_pbmc$seurat_clusters %in% c(0,1,2,8,17)] = 'doublets'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'classical_monocytes' & 
                 sgs_pbmc$seurat_clusters %in% c(0,4,8,11,16,18,20,21,22,23,26,38)] = 'doublets'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'gd_t_cells' & 
                 sgs_pbmc$seurat_clusters %in% c(16)] = 't_cd8'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'gd_t_cells' & 
                 sgs_pbmc$seurat_clusters %in% c(2,6,9,11,18,21,22,29,40,42,45)] = 'doublets'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'ilc' & 
                 sgs_pbmc$seurat_clusters %in% c(8,11,16,18,28)] = 'doublets'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'mait_cells' & 
                 sgs_pbmc$seurat_clusters %in% c(4,9,11,22,25,31,32,40)] = 'doublets'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'nkt' & 
                 sgs_pbmc$seurat_clusters %in% c(9,18)] = 'doublets'
sgs_pbmc$annot[sgs_pbmc$celltype_lvl2_snakecase == 'nkt' & 
                 sgs_pbmc$seurat_clusters %in% c(15)] = 'gd_t_cells'
sgs_pbmc$annot[sgs_pbmc$seurat_clusters %in% c(32) | sgs_pbmc$annot == 'cycling_cells'] = 'cycling_doublets'


sgs_pbmc$annot[grepl('cd4',sgs_pbmc$annot)] = 't_cd4'
sgs_pbmc$annot[grepl('cd8',sgs_pbmc$annot)] = 't_cd8'

checkmate::assert_true(all(sgs_pbmc$annot %in% CELLTYPE_HIERARCHY$celltype_lvl2_snakecase))
unique(sgs_pbmc$annot[!sgs_pbmc$annot %in% CELLTYPE_HIERARCHY$celltype_lvl2_snakecase])

sgs_pbmc$broad_annot = CELLTYPE_HIERARCHY$lineage[match(sgs_pbmc$annot,CELLTYPE_HIERARCHY$celltype_lvl2_snakecase)]
checkmate::assert_true(all(!is.na(sgs_pbmc$broad_annot)))

top_markers <- unique(unlist(HEME_MARKERS))
DotPlot(sgs_pbmc,group.by = 'annot',features = top_markers)+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1,size=6))

saveRDS(sgs_pbmc,'Results/03_scRNAseq_annotation/SGS_inhouse_Wang21_annot_2601.RDS')

### SGS_inhouse ----------------------------------------------------------------
cols_to_add = setdiff(colnames(sgs_pbmc@meta.data),colnames(sgs_inhouse_annot@meta.data))
checkmate::assert_true(all(sgs_inhouse_annot$cellID %in% sgs_pbmc$cellID))
sgs_inhouse_annot@meta.data = cbind(sgs_inhouse_annot@meta.data,
                                    sgs_pbmc@meta.data[match(sgs_inhouse_annot$cellID,sgs_pbmc$cellID),cols_to_add])

DimPlot(sgs_inhouse_annot,group.by = 'annot',label = T,repel = T,label.box = T,cols = col25)+ NoLegend()
DimPlot(sgs_inhouse_annot,cells.highlight = 
          sgs_inhouse_annot$cellID[sgs_inhouse_annot$annot %in% c('unknown','doublets') &
                                     sgs_inhouse_annot$seurat_clusters == 39])


t_cd4 = 0,1,
t_cd8 = 11,14
naive_b = 2,19,29,34
nk = 5
memory_b = 8
classical_monocytes = 15
mk = 16,39
nk_ilc = 17,21
18,9,27,34,36 - unsure
table(sgs_inhouse_annot$seurat_clusters[sgs_inhouse_annot$annot %in% c('unknown','doublets')])
sgs_inhouse_annot$annot_tmp = as.character(sgs_inhouse_annot$annot)
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(0,1)] = 'd_t_cd4'
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(11,14)] = 'd_t_cd8'
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(2,19,29,34)] = 'd_naive_b'
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(5)] = 'd_nk'
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(8)] = 'd_memory_b'
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(15)] = 'd_c_mono'
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(16,39)] = 'd_mk'
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(17)] = 'd_ilc'
sgs_inhouse_annot$annot_tmp[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                              sgs_inhouse_annot$seurat_clusters %in% c(18,9,27,34,36)] = paste0('d_',
                                                                                                sgs_inhouse_annot$seurat_clusters[sgs_inhouse_annot$annot %in% c('unknown','doublets') & 
                                                                                                                              sgs_inhouse_annot$seurat_clusters %in% c(18,9,27,34,36)])

DimPlot(sgs_inhouse_annot,cells.highlight = 
          sgs_pbmc$cellID[sgs_pbmc$annot == 'doublets'])

saveRDS(sgs_inhouse_annot,'Results/03_scRNAseq_annotation/SGS_inhouse_annot_2601.RDS')

# Remove unknown / doublets
sgs_inhouse_annot_clean = subset(sgs_inhouse_annot,
                                 subset = cellID %in% sgs_inhouse_annot$cellID[
                                   !sgs_inhouse_annot$annot %in% c('cycling_doublets','doublets','unknown')])
sgs_inhouse_annot_clean = standard_clustering(sgs_inhouse_annot_clean)
sgs_inhouse_annot_clean = RunUMAP(sgs_inhouse_annot_clean, dims=seq(70))
DimPlot(sgs_inhouse_annot_clean,group.by = 'annot',label = T,repel = T,
        label.box = T,cols = col25)+ NoLegend()

saveRDS(sgs_inhouse_annot_clean,'Results/03_scRNAseq_annotation/SGS_inhouse_annot_clean_2601.RDS')

### Leukaemia ------------------------------------------------------------------
table(leuk_srat$celltype_lvl2_snakecase)

DimPlot(leuk_srat,group.by = 'annot',label = T,repel = T,
        label.box = T,cols = col25)+NoLegend()

DimPlot(leuk_srat,cells.highlight = 
          leuk_srat$cellID[leuk_srat$annot == 'doublets' & 
                             leuk_srat$seurat_clusters %in% c(15)  ])
leuk_srat$annot_tmp = as.character(leuk_srat$annot)

DimPlot(leuk_srat,cells.highlight = 
          l061$cellID[l061$celltype_lvl2_snakecase == 'gd_t_cells'])

FeaturePlot(leuk_srat,c('CD3D','CD8A','CD4'))
DimPlot(leuk_srat,cells.highlight = 
          leuk_srat$cellID[leuk_srat$seurat_clusters == 40])

table(leuk_srat$seurat_clusters[leuk_srat$annot=='memory_b'])
table(leuk_srat$celltype_lvl2_snakecase[leuk_srat$seurat_clusters==32],
      leuk_srat$donor_id[leuk_srat$seurat_clusters==32])
table(leuk_srat$donor_id[leuk_srat$annot=='gd_t_cells'],
      leuk_srat$seurat_clusters[leuk_srat$annot=='gd_t_cells'])


leuk_srat$annot = as.character(leuk_srat$celltype_lvl2_snakecase)
leuk_srat$annot[leuk_srat$celltype_lvl2_snakecase == 'dc_precursor'] = 'dc2'
leuk_srat$annot[leuk_srat$celltype_lvl2_snakecase == 'dc3'] = 'doublets'
leuk_srat$annot[leuk_srat$celltype_lvl2_snakecase == 'classical_monocytes' & 
                  leuk_srat$seurat_clusters %in% c(1,2,3,12,16)] = 'doublets'
leuk_srat$annot[leuk_srat$celltype_lvl2_snakecase == 'gd_t_cells' & 
                  leuk_srat$seurat_clusters %in% c(4)] = 'doublets'
leuk_srat$annot[leuk_srat$celltype_lvl2_snakecase == 'pdc' & 
                  leuk_srat$seurat_clusters %in% c(4)] = 'doublets'
leuk_srat$annot[leuk_srat$annot == 'naive_b' & 
                   leuk_srat$seurat_clusters %in% c(1,6,9,14,16)] = 'doublets'
leuk_srat$annot[leuk_srat$annot == 'memory_b' & 
                  leuk_srat$seurat_clusters %in% c(4)] = 'doublets'
leuk_srat$annot[grepl('cd4',leuk_srat$annot)] = 't_cd4'
leuk_srat$annot[grepl('cd8',leuk_srat$annot)] = 't_cd8'
leuk_srat$annot[leuk_srat$annot == 'doublets' & leuk_srat$seurat_clusters %in% c( 0,2,3,5,15)] = 't_cd4'
leuk_srat$annot[leuk_srat$annot == 'doublets' & leuk_srat$seurat_clusters %in% c(4)] = 'classical_monocytes'
leuk_srat$annot[leuk_srat$annot == 'doublets' & leuk_srat$seurat_clusters %in% c(6)] = 'dc2'
leuk_srat$annot[leuk_srat$annot == 'doublets' & leuk_srat$seurat_clusters %in% c(16)] = 'pdc'
leuk_srat$annot[leuk_srat$annot == 'doublets' & leuk_srat$seurat_clusters %in% c(12)] = 'memory_b'
leuk_srat$annot[leuk_srat$annot == 'doublets' & leuk_srat$seurat_clusters %in% c(13,14)] = 'naive_b'
leuk_srat$annot[leuk_srat$annot == 'doublets' & leuk_srat$seurat_clusters %in% c(10)] = 'ilc'


checkmate::assert_true(all(leuk_srat$annot %in% CELLTYPE_HIERARCHY$celltype_lvl2_snakecase))
unique(leuk_srat$annot[!leuk_srat$annot %in% CELLTYPE_HIERARCHY$celltype_lvl2_snakecase])

leuk_srat$broad_annot = CELLTYPE_HIERARCHY$lineage[match(leuk_srat$annot,CELLTYPE_HIERARCHY$celltype_lvl2_snakecase)]
checkmate::assert_true(all(!is.na(leuk_srat$broad_annot)))

top_markers <- unique(unlist(HEME_MARKERS))
DotPlot(leuk_srat,group.by = 'annot_tmp',features = top_markers)+
  RotatedAxis()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1,size=6))

saveRDS(leuk_srat,'Results/03_scRNAseq_annotation/MDS_L061_L067_annot_2601.RDS')

# Remove unknown / doublets
leuk_srat_clean = subset(leuk_srat,
                         subset = cellID %in% leuk_srat$cellID[!leuk_srat$annot %in% c('doublets','unknown')])
leuk_srat_clean = standard_clustering(leuk_srat_clean)
leuk_srat_clean = RunUMAP(leuk_srat_clean, dims=seq(70),seed.use = 2397)
DimPlot(leuk_srat_clean,group.by = 'annot',label = T,repel = T,
        label.box = T,cols = col25)+ NoLegend()

saveRDS(leuk_srat_clean,'Results/03_scRNAseq_annotation/MDS_L061_L067_annot_clean_2601.RDS')

### All the data together ------------------------------------------------------
cols_to_add = setdiff(unique(colnames(leuk_srat@meta.data),colnames(sgs_pbmc@meta.data)),
                      colnames(combined_srat@meta.data))
df = rbind(leuk_srat@meta.data[,c('cellID',cols_to_add)],
           sgs_pbmc@meta.data[,c('cellID',cols_to_add)])
checkmate::assert_true(all(combined_srat$cellID %in% df$cellID))
combined_srat@meta.data = cbind(combined_srat@meta.data,df[match(combined_srat$cellID,df$cellID),cols_to_add])

table(combined_srat$annot)
DimPlot(combined_srat,group.by = 'annot',label = T,repel = T,label.box = T,cols = col25)+ NoLegend()
DimPlot(sgs_pbmc,cells.highlight = 
          combined_srat$cellID[combined_srat$annot == 'unknown'])

table(combined_srat$batch,combined_srat$donor_id)

inhouse_srat = subset(combined_srat,subset = batch !='external')

# OLD --------------------------------------------------------------------------
# outDir = 'Results/05_combined_scRNAseq_annotation'
# sc_datasets_fp = c(sgs_fp = 'Results/03_SETBP1_annotation/2505/SETBP1_SGS_clean_noMTCells_annot_2505.RDS',
#                    wang21_fp = 'Results/04_public_scRNAseq/2601/Wang21__rhoLimNone__clean_noMTCells.RDS',
#                    mds_l061 = 'Results/04_public_scRNAseq/2601/L061/L061__clean_noMTCells.RDS',
#                    mds_l067 = 'Results/04_public_scRNAseq/2601/L067/L067__clean_noMTCells.RDS'
# )
#
# checkmate::assert_true(all(file.exists(sc_datasets_fp)))
#
# if(!dir.exists(outDir)){
#   dir.create(outDir,recursive = T)
# }
#
# # Import datasets --------------------------------------------------------------
# c(sgs_srat, wang21_srat, l061_srat, l067_srat) %<-%
#   lapply(sc_datasets_fp,function(x){readRDS(x)})
#
# sgs_pbmc = merge_seurat_objects(sgs_srat,wang21_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# leuk_srat = merge_seurat_objects(l061_srat,l067_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
# combined_srat = merge_seurat_objects(sgs_pbmc,leuk_srat,keepAllGenes = F,genomeVersions = c('v38','v38'))
#
# c(sgs_pbmc, leuk_srat, combined_srat) %<-%
#   lapply(c(sgs_pbmc, leuk_srat, combined_srat),function(x){standard_clustering(x)})
#
# # Add metadata to each dataset -------------------------------------------------
# sample_metadata = readxl::read_excel('Results/Supplementary Table 1.xlsx',skip = 2) %>%
#   janitor::clean_names() %>%
#   dplyr::mutate(batch = dplyr::case_when(donor_id == 'BJ9' ~ 'batch_3',
#                                          donor_id == 'GOSH084' & sample_id %in% c('NB11528275','NB11528276','NB11528277','NB11528278') ~ 'batch_1',
#                                          donor_id == 'GOSH084' & sample_id %in% c('Leuk13234209','Leuk13234210','Leuk13234211') ~ 'batch_2',
#                                          donor_id %in% c('BJ111','BJ112','BJ113','BJ114') ~ 'batch_4',
#                                          donor_id == 'L061' ~ 'batch_5',
#                                          donor_id == 'L067' ~ 'batch_6',
#                                          .default = 'external'
#                                          ))
# table(sample_metadata$donor_id,sample_metadata$batch)
#
# # add sample-level metadata
# c(sgs_pbmc, leuk_srat, combined_srat) %<-%
#   lapply(c(sgs_pbmc, leuk_srat, combined_srat),function(srat){
#     sample_mdat_columns_to_keep = unique(c('sample_id',setdiff(colnames(sample_metadata),c(colnames(srat@meta.data),'number_of_cells_post_qc'))))
#
#     srat$sample_id = gsub('^SB\\.','',srat$orig.ident)
#     checkmate::assert_true(all(srat$cellID == rownames(srat@meta.data)))
#     checkmate::assert_true(all(srat$sample_id %in% sample_metadata$sample_id))
#
#     srat@meta.data = srat@meta.data %>%
#       dplyr::left_join(sample_metadata[,sample_mdat_columns_to_keep], by = 'sample_id')
#     rownames(srat@meta.data) = srat@meta.data$cellID
#     srat
#   })
#
# # Cell-type annotation ---------------------------------------------------------
# ## helper functions ------------------------------------------------------------
#
#
# res <- score_clusters_markers(
#   sgs_pbmc,
#   pbmc_markers_short
# )
#
# res$scores          # cluster × cell-type
# res$predicted       # best hit per cluster
# res$ambiguity       # confidence proxy
#
# seurat_obj$marker_celltype <- res$predicted[as.character(seurat_obj$seurat_clusters)]
# seurat_obj$marker_ambiguity <- res$ambiguity[as.character(seurat_obj$seurat_clusters)]
# seurat_obj <- AddModuleScore(
#   seurat_obj,
#   features = marker_list,
#   name = "ModuleScore"
# )
#
# score_cols <- grep("^ModuleScore", colnames(seurat_obj@meta.data), value = TRUE)
# names(score_cols) <- names(marker_list)
#
# DimPlot(seurat_obj,group.by = 'seurat_clusters',label = T,label.box = T)
# DimPlot(seurat_obj,group.by = 'marker_celltype',label = T,label.box = T) + NoLegend()
# FeaturePlot(seurat_obj,'marker_ambiguity')
#
# # Sub-cluster by batch to get better resolution
# table(seurat_obj$batch)
# seurat_list <- SplitObject(
#   seurat_obj,
#   split.by = "batch"
# )
# seurat_list <- lapply(seurat_list, standard_clustering)
# i=1
# DimPlot(seurat_list[[i]],group.by = c('seurat_clusters'),
#         label = T,label.box = T,label.size = 2) +
#   NoLegend()+ggtitle(names(seurat_list)[i])
#
# table(seurat_list[[i]]$marker_celltype[seurat_list[[i]]$seurat_clusters == 8])
#
# library(patchwork)
# out_pdf <- "Results/05_combined_scRNAseq_annotation/sgs_pbmc_batch_subcluster_qc.pdf"
# # Estimate height dynamically (helps with many marker sets)
# dotplot_height <- max(4, length(pbmc_markers_short) * 0.25)
# pdf(out_pdf, width = 13, height = 6 + dotplot_height)
#
# for (i in seq_along(seurat_list)) {
#
#   obj <- seurat_list[[i]]
#   batch_name <- names(seurat_list)[i]
#
#   ## ---- DimPlots ----
#   p_clusters <- DimPlot(
#     obj,
#     reduction = "umap",
#     group.by = "seurat_clusters",
#     repel = TRUE,
#     label = TRUE,
#     label.box = TRUE,
#     label.size = 3
#   ) +
#     ggtitle("Seurat clusters") +
#     theme_bw(base_size = 12)+
#     NoLegend()
#
#   p_celltypes <- DimPlot(
#     obj,
#     reduction = "umap",
#     group.by = "marker_celltype",
#     repel = TRUE,
#     label = TRUE,
#     label.box = TRUE,
#     label.size = 3
#   ) +
#     ggtitle("Marker cell types") +
#     theme_bw(base_size = 12)+
#     NoLegend()
#
#   top_row <- p_clusters + p_celltypes + plot_layout(ncol = 2)
#
#   ## ---- DotPlot ----
#   Idents(obj) <- "seurat_clusters"
#
#   p_dot <- DotPlot(
#     obj,
#     features = unique(unlist(pbmc_markers_short, use.names = FALSE)),
#     dot.scale = 6
#   ) +
#     RotatedAxis() +
#     scale_color_viridis_c(option = "C") +
#     ggtitle("Marker expression by cluster") +
#     theme_bw(base_size = 10) +
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1),
#       panel.grid.major = element_blank(),
#       panel.grid.minor = element_blank()
#     )
#
#   ## ---- Combine ----
#   combined <- top_row / p_dot +
#     plot_layout(heights = c(1, 1.2)) +
#     plot_annotation(
#       title = paste("Batch:", batch_name),
#       theme = theme(
#         plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
#       )
#     )
#
#   print(combined)
# }
#
# dev.off()
#
#
# for (b in names(seurat_list)) {
#   seurat_obj$subcluster_batch[Cells(seurat_list[[b]])] <-
#     paste0(b, "_", seurat_list[[b]]$seurat_clusters)
# }