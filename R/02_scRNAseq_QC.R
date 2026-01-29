##============================================================================##
##                    Single-cell RNA-seq Quality Control                     ##
##============================================================================##

# Libraries --------------------------------------------------------------------
library(tidyverse)
library(Seurat)
source("R/utils/sc_basicQC.R")

# Global parameters ------------------------------------------------------------
seed <- 2397
set.seed(seed)

# Helper functions -------------------------------------------------------------

#' Setup output directory
#' @param outDir Output directory path
setup_output_dir <- function(outDir) {
  if (!dir.exists(outDir)) {
    message(sprintf('Creating output directory: %s', outDir))
    dir.create(outDir, recursive = TRUE)
  }
}

#' Run QC pipeline wrapper
#' @param dataDirs Named vector of data directories
#' @param outPath Output path prefix
#' @param plotDir Plot directory prefix
#' @param qc_params List of QC parameters
#' @param metadata Optional metadata dataframe
#' @param matchBy Optional matching columns
#' @param cleanCountDir Optional clean count directory.
#' Specify NULL to save scrublet and soupX results in the same dataDir 
#' @return List containing cleaned Seurat object and QC summary
run_qc_pipeline <- function(dataDirs, outPath, plotDir, qc_params, 
                             metadata = NULL, matchBy = NULL, 
                             cleanCountDir = NULL) {
  
  # Check if output already exists
  cleanSrat_fp <- ifelse(
    qc_params$keepMTCells,
    paste0(outPath, '_clean_withMTCells.RDS'),
    paste0(outPath, '_clean_noMTCells.RDS')
  )
  
  if (file.exists(cleanSrat_fp) && qc_params$skipIfExists) {
    message(sprintf('Loading existing clean Seurat object from: %s', cleanSrat_fp))
    cleanSrat <- readRDS(cleanSrat_fp)
    return(list(seurat = cleanSrat, summary = NULL))
  }
  
  # Filter valid directories
  dataDirs <- dataDirs[file.exists(dataDirs)]
  dataDirs <- dataDirs[sapply(dataDirs, function(x) length(list.files(x)) > 0)]
  
  message(sprintf('Processing %d samples...', n_distinct(dataDirs)))
  
  # Run basicQC
  message('\nPerforming scRNAseq QC...')
  QC.output <- basicQC(
    dataDirs = dataDirs,
    maxMT = qc_params$maxMT,
    minGenes = qc_params$minGenes,
    minUMIs = qc_params$minUMIs,
    maxBadFrac = qc_params$maxBadFrac,
    numPCs = qc_params$numPCs,
    clusteringRes = qc_params$clusteringRes,
    cleanCountDir = cleanCountDir,
    skipScrub = qc_params$skipScrub,
    skipSoup = qc_params$skipSoup,
    scrubScoreMax = qc_params$scrubScoreMax,
    scrubPath = qc_params$scrubPath,
    metadata = metadata,
    matchBy = matchBy,
    scPath = qc_params$scPath,
    rho_max_limit = qc_params$rho_max_limit,
    outPath = outPath,
    skipIfExists = qc_params$skipIfExists,
    doPlot = qc_params$doPlot,
    plotDir = plotDir,
    verbose = qc_params$verbose,
    is10X = qc_params$is10X
  )
  
  cleanSrat <- QC.output[[1]]
  df.out <- QC.output[[2]]
  
  # Save QC summary
  write.csv(df.out, paste0(outPath, '_qc_summary.csv'), row.names = FALSE)
  
  return(list(seurat = cleanSrat, summary = df.out))
}

##============================================================================##
##                         DATASET 1: SETBP1 In-house                         ##
##============================================================================##


# QC SETBP1 single-cell RNASeq dataset
process_setbp1_inhouse <- function() {
  message("\n==================== Processing SETBP1 In-house Dataset ====================\n")
  
  outDir <- 'Results/02_SETBP1_scQC/2505'
  setup_output_dir(outDir)
  
  params <- QC_PARAMS(minGenes = 200,
                      minUMIs = 300,
                      skipSoup = FALSE,
                      rho_max_limit = NULL)
  
  # Setup data directories
  dataDirs = file.path(list.files('Data/inhouse_SETBP1_10X',full.names = T,pattern = 'GRCh38-2020-A'),'/filtered_feature_bc_matrix')
  # Remove MDS-SETBP1 sample
  dataDirs = dataDirs[!grepl('SB_Leuk13645525',dataDirs)]
  # Remove BJ9-solidTumour snRNA-seq samples
  dataDirs = dataDirs[!grepl('MY_200531_13839752|MY_200531_13839753|MY_200531_13839754',dataDirs)]    
  
  names(dataDirs) = gsub('^.*_CG_|_GRCh38-2020-A.*$','',dataDirs)    
  names(dataDirs) = gsub('^.*_MY_','MY_',names(dataDirs))    
  names(dataDirs) = gsub('^.*_SB_','SB_',names(dataDirs))
  names(dataDirs) = gsub('_','.',names(dataDirs))    
  
  dataDirs=dataDirs[file.exists(dataDirs)]
  dataDirs=dataDirs[sapply(dataDirs, function(x){length(list.files(x))>0})]
  print(n_distinct(dataDirs))
  checkmate::assert_true(length(dataDirs) > 0)
  # Run QC pipeline
  outPath <- file.path(outDir, 'SETBP1_SGS')
  plotDir <- file.path(outDir, 'SETBP1_SGS_')
  
  result <- run_qc_pipeline(
    dataDirs = dataDirs,
    outPath = outPath,
    plotDir = plotDir,
    qc_params = params,
    cleanCountDir = dataDirs
  )
  
  message("SETBP1 In-house dataset processing complete.\n")
  return(result$seurat)
}

##============================================================================##
##                         DATASET 2: L061 SETBP1 Leukemia                    ##
##============================================================================##

process_l061_setbp1 <- function() {
  message("\n============== Processing L061 SETBP1 Dataset =================\n")
  
  outDir <- 'Results/04_public_scRNAseq/2601/L061'
  setup_output_dir(outDir)
  
  # Dataset-specific parameters
  params <- QC_PARAMS(rho_max_limit = 0.02)
  
  # Setup data directory
  dataDirs <- 'Data/L061_pleukPBMC/cellranger612_count_47089_SB_Leuk13645525_GRCh38-2020-A/filtered_feature_bc_matrix'
  names(dataDirs) <- gsub('.*SB_|_GR.*$', '', basename(dirname(dataDirs)))
  
  # Load metadata
  projMan <- readxl::read_excel('~/lustre_mt22/projectManifest.xlsx', sheet = 'SETBP1')
  projMan <- projMan[!is.na(projMan$assay) & grepl('scRNA', projMan$assay) & 
                       projMan$DonorID == 'L061', ]
  
  metadata <- projMan[, c('DonorID', 'externalID', 'Tissue', 'Sex', 'assay', 'sangerSampleID')]
  colnames(metadata) <- c('donorID', 'externalID', 'tissue', 'sex', 'assay', 'sangerSampleID')
  metadata$channelID <- gsub('^SB_', '', metadata$sangerSampleID)
  
  # Run QC pipeline
  outPath <- file.path(outDir, 'L061')
  plotDir <- file.path(outDir, 'L061_')
  
  result <- run_qc_pipeline(
    dataDirs = dataDirs,
    outPath = outPath,
    plotDir = plotDir,
    qc_params = params,
    metadata = metadata,
    matchBy = c('orig.ident', 'channelID'),
    cleanCountDir = NULL
  )
  
  message("L061 SETBP1 dataset processing complete.\n")
  return(result$seurat)
}

##============================================================================##
##                         DATASET 3: L067 GATA1                              ##
##============================================================================##

process_l067_gata1 <- function() {
  message("\n==================== Processing L067 GATA1 Dataset ====================\n")
  
  outDir <- 'Results/04_public_scRNAseq/2601/L067'
  setup_output_dir(outDir)
  
  # Dataset-specific parameters (same as L061)
  params <- QC_PARAMS(rho_max_limit = 0.02)
  
  # Setup data directory
  dataDirs <- 'Data/L067_pleukPBMC/cellranger700_count_47089_SB_Leuk13645530_GRCh38-2020-A/filtered_feature_bc_matrix'
  names(dataDirs) <- gsub('.*SB_|_GR.*$', '', basename(dirname(dataDirs)))
  
  # Load metadata
  projMan <- readxl::read_excel('~/lustre_mt22/projectManifest.xlsx', sheet = 'GOSH_others_included')
  projMan <- projMan[!is.na(projMan$assay) & !grepl('WGS', projMan$assay) & 
                       projMan$donorID == 'L067', ]
  
  metadata <- projMan[, c('donorID', 'Tissue', 'Sex', 'assay', 'sangerSampleID')]
  colnames(metadata) <- c('donorID', 'tissue', 'sex', 'assay', 'sangerSampleID')
  metadata$channelID <- gsub('^SB_', '', metadata$sangerSampleID)
  
  # Run QC pipeline
  outPath <- file.path(outDir, 'L067')
  plotDir <- file.path(outDir, 'L067_')
  
  result <- run_qc_pipeline(
    dataDirs = dataDirs,
    outPath = outPath,
    plotDir = plotDir,
    qc_params = params,
    metadata = metadata,
    matchBy = c('orig.ident', 'channelID'),
    cleanCountDir = NULL
  )
  
  message("L067 GATA1 dataset processing complete.\n")
  return(result$seurat)
}

##============================================================================##
##                         DATASET 4: Wang et al. 2021                        ##
##============================================================================##

process_wang2021 <- function(rho_max_limit = NULL) {
  message("\n============== Processing Wang et al. 2021 Dataset ============\n")
  
  outDir <- 'Results/04_public_scRNAseq/2601'
  setup_output_dir(outDir)
  
  # Dataset-specific parameters
  params <- QC_PARAMS(skipSoup = TRUE,
                      rho_max_limit = rho_max_limit)
  
  # Setup paths based on rho_max_limit
  if (!is.null(rho_max_limit)) {
    outPath <- paste0(file.path(outDir, 'Wang21_rhoLim'), rho_max_limit, '_')
    plotDir <- paste0(file.path(outDir, 'Wang21_rhoLim'), rho_max_limit, '_')
  } else {
    outPath <- file.path(outDir, 'Wang21_rhoLimNone_')
    plotDir <- file.path(outDir, 'Wang21_rhoLimNone_')
  }
  
  # Setup data directories
  dataDirs <- list.dirs('Data/published_scRNAseq/Wang21')
  dataDirs <- dataDirs[grepl('H\\d$', dataDirs)]
  names(dataDirs) <- basename(dataDirs)
  
  # Run QC pipeline
  result <- run_qc_pipeline(
    dataDirs = dataDirs,
    outPath = outPath,
    plotDir = plotDir,
    qc_params = params,
    cleanCountDir = NULL
  )
  
  message("Wang et al. 2021 dataset processing complete.\n")
  return(result$seurat)
}

##============================================================================##
##                         MAIN EXECUTION                                     ##
##============================================================================##

main <- function(datasets = c("setbp1", "l061", "l067", "wang2021")) {
  
  results <- list()
  
  if ("setbp1" %in% datasets) {
    results$setbp1 <- tryCatch(
      process_setbp1_inhouse(),
      error = function(e) {
        message("Error processing SETBP1 dataset: ", e$message)
        NULL
      }
    )
  }
  
  if ("l061" %in% datasets) {
    results$l061 <- tryCatch(
      process_l061_setbp1(),
      error = function(e) {
        message("Error processing L061 dataset: ", e$message)
        NULL
      }
    )
  }
  
  if ("l067" %in% datasets) {
    results$l067 <- tryCatch(
      process_l067_gata1(),
      error = function(e) {
        message("Error processing L067 dataset: ", e$message)
        NULL
      }
    )
  }
  
  if ("wang2021" %in% datasets) {
    results$wang2021 <- tryCatch(
      process_wang2021(rho_max_limit = NULL),
      error = function(e) {
        message("Error processing Wang2021 dataset: ", e$message)
        NULL
      }
    )
  }
  
  message("\n=================== All Datasets Processed ====================\n")
  return(results)
}

# OLD --------------------------------------------------------------------------

##-------------------------------##
##   SETBP1 Global parameters  ####
##-------------------------------##

seed = 2397


maxMT = 30
minGenes = 200
minUMIs = 300
maxBadFrac = 0.5
numPCs = 75
clusteringRes = 10
skipScrub = F
skipSoup=F
scrubScoreMax = 0.5
scrubPath='./scrubletScores.tsv'
scPath="./strainedCounts"
doPlot=T
verbose = T
skipIfExists=F
keepMTCells=T

