library(Seurat)
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes

# Seurat object handling -------------------------------------------------------
standard_clustering = function(srat,nPCs=75,clusteringRes=1,skipCCS=FALSE,s.genes=NULL,g2m.genes=NULL,runHarmony=F,harmonyVar = NULL,doPlot=T,...){
  require(Seurat)
  require(cowplot)
  srat = NormalizeData(srat)
  srat = FindVariableFeatures(srat)
  srat = ScaleData(srat,...)
  srat = RunPCA(srat, npcs = nPCs)
  
  if(runHarmony){
    # check the harmony variables provided
    if(is.null(harmonyVar) | !all(harmonyVar %in% colnames(srat@meta.data))){
      stop('No or not all variables for Harmony provided exist...')
    }else{
      require(harmony)
      
      srat = RunHarmony(srat,harmonyVar, plot_convergence = TRUE)
      # To directly access the new Harmony embeddings, use the Embeddings command.
      harmony_embeddings <- Embeddings(srat, 'harmony')
    
      if(doPlot){
        options(repr.plot.height = 5, repr.plot.width = 12)
        p1 <- DimPlot(object = srat, reduction = "harmony", pt.size = .1, group.by = harmonyVar[1])
        p2 <- VlnPlot(object = srat, features = "harmony_1", group.by = harmonyVar[1], pt.size = .1)
        plot_grid(p1,p2)  
      }
    }
  }
  #ElbowPlot(srat, ndims = nPCs)

  # Always cluster without Harmony first
  srat = FindNeighbors(srat, dims=seq(nPCs))
  srat = FindClusters(srat, resolution = clusteringRes)
  srat = RunUMAP(srat, dims=seq(nPCs))
  
  if(!runHarmony | is.null(harmonyVar)){doPlot=F}
  if(doPlot){
    p = list()
    for(var in harmonyVar){
      p[[var]] = DimPlot(object = srat, reduction = "umap", pt.size = .1, group.by = var) + ggtitle('Pre-Harmony')
    }
    plot_grid(plotlist = p)
  }
  
  if(runHarmony){
    srat = FindNeighbors(srat, dims=seq(nPCs), reduction = ifelse(runHarmony,"harmony",'pca'))
    srat = FindClusters(srat, resolution = clusteringRes, reduction = ifelse(runHarmony,"harmony",'pca'))
    srat = RunUMAP(srat, dims=seq(nPCs), reduction = ifelse(runHarmony,"harmony",'pca'))
    
    if(doPlot){
      pHarm = list()
      for(var in harmonyVar){
        pHarm[[var]] = DimPlot(object = srat, reduction = "umap", pt.size = .1, group.by = var) + ggtitle('Post-Harmony')
      }
      plot_grid(plotlist = pHarm)
    }
  }
  


  # CellCycle Scoring
  if(!skipCCS){
    if(is.null(s.genes)){
      s.genes <- cc.genes$s.genes  
    }
    if(is.null(g2m.genes)){
      g2m.genes <- cc.genes$g2m.genes  
    }
    
    srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)  
  }
  
  return(srat)
}


## Most of the time, Seurat::merge() is fine. But because it keeps the union of genes between the 2 objects
# this might be problematic if the 2 objects were generated and pre-processed differently (eg. published external data or different genome versions etc.)
# also, should always match on geneID instead of gene names

merge_seurat_objects = function(srat1, srat2, keepAllGenes=F, genomeVersions = NULL){
  require(tidyverse)
  ## Check genome version used to map data from each of the 2 objects
  if(is.null(genomeVersions) | length(genomeVersions) != 2){
    stop('Please provide genome versions used for the 2 seurat objects')
  }
  
  ## Map gene names via ensID
  if(genomeVersions[1] == genomeVersions[2]){
    message(sprintf('Same genome version (%s) was used for both objects - No further geneID processing is needed',unique(genomeVersions)))
  }else if(genomeVersions[1] != genomeVersions[2]){
    message(sprintf('Mapping geneID for both objects to v38'))
  }
  
  ## use Seurat function to generate merged metadata
  srat1@meta.data$cellID = rownames(srat1@meta.data)
  srat2@meta.data$cellID = rownames(srat2@meta.data)
  merged.mdat = plyr::rbind.fill(srat1@meta.data,srat2@meta.data)
  
  if(sum(duplicated(merged.mdat$cellID)) > 0){
    rownames(merged.mdat) =  make.names(merged.mdat$cellID, unique=TRUE)  
  }else{
    rownames(merged.mdat) = merged.mdat$cellID
  }
  
  #rownames(merged.mdat) =  gsub('^X','',rownames(merged.mdat))
  #sum(duplicated(colnames(srat2)))
  
  
  # If we wish to check that this is equivalent to what Seurat::merge() does to merge metadata, here's the code for it
  # I've checked though and its good!
  #require(Seurat)
  #merged.srat = merge(srat1,srat2)
  #merged.mdat2 = merged.srat@meta.data
  #merged.mdat2$cellID = rownames(merged.mdat2)
  #indx <- sapply(merged.mdat, is.factor)
  #merged.mdat[indx] <- lapply(merged.mdat[indx], function(x) as.character(x))
  #dplyr::all_equal(merged.mdat,merged.mdat2)
  require(Seurat)
  
  if(keepAllGenes){
    # Use Seurat::merge() method
    merged.srat = merge(srat1,srat2)
    return(merged.srat)
  }
  
  ## If NOT keepAllGenes - ie. only keep genes that are common between the 2 srat objects. 
  
  ## get common genes
  # removing genes with 0 counts everywhere
  gene1 = rownames(srat1@assays$RNA@counts)
  gene2 = rownames(srat2@assays$RNA@counts)
  common_genes = intersect(gene1,gene2)
  
  ## Calculate the fraction of total counts attributed to those genes that got removed
  srat1_unique_genes = rownames(srat1@assays$RNA@counts)[!rownames(srat1@assays$RNA@counts) %in% common_genes]
  srat1_frac_count_lost = sum(srat1@assays$RNA@counts[rownames(srat1@assays$RNA@counts) %in% srat1_unique_genes,])/sum(srat1@assays$RNA@counts)  
  
  srat2_unique_genes = rownames(srat2@assays$RNA@counts)[!rownames(srat2@assays$RNA@counts) %in% common_genes]
  srat2_frac_count_lost = sum(srat2@assays$RNA@counts[rownames(srat2@assays$RNA@counts) %in% srat2_unique_genes,])/sum(srat2@assays$RNA@counts)  
  
  message(sprintf('# genes in srat1: %d \n# genes in srat2: %d \n# genes in common: %d\n\nFrac count lost in srat1: %f \nFrac count lost in srat2: %f \n',
                  n_distinct(rownames(srat1@assays$RNA@counts)),
                  n_distinct(rownames(srat2@assays$RNA@counts)),
                  n_distinct(common_genes),
                  srat1_frac_count_lost,
                  srat2_frac_count_lost))
  
  
  
  cnt.mtx1 = srat1@assays$RNA@counts[rownames(srat1@assays$RNA@counts) %in% common_genes,]
  cnt.mtx2 = srat2@assays$RNA@counts[rownames(cnt.mtx1),]
  if(nrow(cnt.mtx1) != nrow(cnt.mtx2)){
    stop('There is issue with the merge!')
  }
  mtx = cbind(cnt.mtx1,cnt.mtx2)
  srat = CreateSeuratObject(mtx,meta.data = merged.mdat)
  
  return(srat)
}

# Annotation helpers -----------------------------------------------------------
ANNOT_PARAMS <- function(
  # CellTypist parameters
  celltypist_model = 'Immune_All_Low',
  celltypist_conf_high = 0.95,
  celltypist_conf_low = 0.85,
  celltypist_diff_threshold = 0.1,
  
  # Marker scoring parameters
  min_pct = 0.1,
  scale_across_clusters = TRUE,
  weight_pct = 1,
  weight_expr = 1,
  
  # Integration / ambiguity handling
  ambiguity_threshold = 0.2,
  
  # Plotting
  dotplot_dot_scale = 6,
  umap_label_size = 3
) {
  
  ## CellTypist parameters
  checkmate::qassert(celltypist_model, "S1")
  checkmate::qassert(celltypist_conf_high, "N(0,1]")
  checkmate::qassert(celltypist_conf_low,  "N(0,1]")
  checkmate::qassert(celltypist_diff_threshold, "N[0,1]")
  
  if (celltypist_conf_low >= celltypist_conf_high) {
    stop("celltypist_conf_low must be < celltypist_conf_high")
  }
  
  ## Marker scoring parameters
  checkmate::qassert(min_pct, "N[0,1]")
  checkmate::qassert(scale_across_clusters, "B1")
  checkmate::qassert(weight_pct,  "N[0,]")
  checkmate::qassert(weight_expr, "N[0,]")
  
  if (weight_pct == 0 && weight_expr == 0) {
    stop("At least one of weight_pct or weight_expr must be > 0")
  }
  
  ## Integration / ambiguity
  checkmate::qassert(ambiguity_threshold, "N[0,1]")
  
  ## Plotting
  checkmate::qassert(dotplot_dot_scale, "N(0,]")
  checkmate::qassert(umap_label_size, "N(0,]")
  
  ## return
  list(
    celltypist_model = celltypist_model,
    celltypist_conf_high = celltypist_conf_high,
    celltypist_conf_low = celltypist_conf_low,
    celltypist_diff_threshold = celltypist_diff_threshold,
    
    min_pct = min_pct,
    scale_across_clusters = scale_across_clusters,
    weight_pct = weight_pct,
    weight_expr = weight_expr,
    
    ambiguity_threshold = ambiguity_threshold,
    
    dotplot_dot_scale = dotplot_dot_scale,
    umap_label_size = umap_label_size
  )
}


#' Score clusters by marker genes
#' @param seurat_obj Seurat object
#' @param marker_list Named list of marker genes
#' @param params Scoring parameters
#' @return List with scores, predictions, and ambiguity
score_clusters_markers <- function(
  seurat_obj,
  marker_list,
  params = ANNOT_PARAMS(),
  assay = "RNA",
  slot = "data"
) {
  message("\n==== Scoring clusters by markers ====")
  
  Idents(seurat_obj) <- "seurat_clusters"
  
  expr <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  clusters <- levels(Idents(seurat_obj))
  
  # Filter to existing genes
  genes <- intersect(unique(unlist(marker_list)), rownames(expr))
  if (length(genes) == 0)
    stop("None of the marker genes found in expression matrix")
  
  message(sprintf("Found %d/%d marker genes in dataset",
                  length(genes), length(unique(unlist(marker_list)))))
  # initialize
  avg_exp <- matrix(0, nrow = length(genes), ncol = length(clusters),
                    dimnames = list(genes, clusters))
  pct_exp <- matrix(0, nrow = length(genes), ncol = length(clusters),
                    dimnames = list(genes, clusters))
  
  # compute cluster-level stats
  for (cl in clusters) {
    cells <- WhichCells(seurat_obj, idents = cl)
    mat <- expr[genes, cells, drop = FALSE]
    
    avg_exp[, cl] <- Matrix::rowMeans(mat)
    pct_exp[, cl] <- Matrix::rowMeans(mat > 0)
  }
  
  # remove uninformative genes
  keep <- matrixStats::rowMaxs(pct_exp) >= params$min_pct
  avg_exp <- avg_exp[keep, , drop = FALSE]
  pct_exp <- pct_exp[keep, , drop = FALSE]
  
  message(sprintf("Retained %d/%d genes after min_pct filter",
                  sum(keep), length(keep)))
  
  # specificity normalization
  if (params$scale_across_clusters) {
    avg_exp <- t(scale(t(avg_exp)))
    avg_exp[is.na(avg_exp)] <- 0
  }
  
  # combined gene-level score
  gene_score <- (params$weight_expr * avg_exp) * (params$weight_pct * pct_exp)
  
  # aggregate per cell type
  celltype_scores <- sapply(marker_list, function(g) {
    g <- intersect(g, rownames(gene_score))
    if (length(g) == 0) return(rep(0, ncol(gene_score)))
    colMeans(gene_score[g, , drop = FALSE])
  })
  
  rownames(celltype_scores) <- clusters
  
  # ambiguity metric
  ambiguity <- apply(celltype_scores, 1, function(x) {
    sx <- sort(x, decreasing = TRUE)
    if (length(sx) < 2) return(NA)
    sx[1] - sx[2]
  })
  
  # Predicted cell type per cluster
  predicted <- apply(celltype_scores, 1, function(x) {
    if (all(x == 0)) return("unknown")
    names(which.max(x))
  })
  
  list(
    scores = as.data.frame(celltype_scores),
    predicted = predicted,
    ambiguity = ambiguity,
    gene_score = gene_score,
    avg_exp = avg_exp,
    pct_exp = pct_exp
  )
}

## CellTypist helpers ----------------------------------------------------------
#' Improved softmax with numerical stability
#' @param logits Matrix of logit scores (cells × cell types)
#' @return Matrix of probabilities
stable_softmax <- function(logits) {
  if (is.null(dim(logits))) {
    # Single vector
    logits <- logits - max(logits)
    exp_logits <- exp(logits)
    return(exp_logits / sum(exp_logits))
  } else {
    # Matrix
    t(apply(logits, 1, function(x) {
      x <- x - max(x)
      exp_x <- exp(x)
      exp_x / sum(exp_x)
    }))
  }
}

#' Assign cell labels from probability matrix
#' @param prob_mat Matrix of probabilities (cells × cell types)
#' @param conf_high High confidence threshold
#' @param conf_low Low confidence threshold
#' @param diff_threshold Minimum difference between top 2 predictions
#' @return Named vector of cell type labels
assign_labels_from_probs <- function(
  prob_mat,
  conf_high = 0.95,
  conf_low = 0.85,
  diff_threshold = 0.1
) {
  
  apply(prob_mat, 1, function(x) {
    sorted_x <- sort(x, decreasing = TRUE)
    max_score <- sorted_x[1]
    second_max_score <- if(length(sorted_x) > 1) sorted_x[2] else 0
    diff <- max_score - second_max_score
    
    top_label <- names(x)[which.max(x)]
    second_label <- if(length(sorted_x) > 1) {
      names(x)[order(x, decreasing = TRUE)[2]]
    } else ""
    
    # Decision tree
    if (sum(x == max_score) > 1) {
      # Tie between multiple cell types
      tied <- names(x[x == max_score])
      return(paste0("ambiguous:", paste(tied, collapse = ":")))
    }
    
    if (max_score >= conf_high && diff >= diff_threshold) {
      return(top_label)
    }
    
    if (max_score >= conf_high && diff < diff_threshold) {
      return(paste0("ambiguous:", top_label, ":", second_label))
    }
    
    if (max_score >= conf_low && diff >= diff_threshold) {
      return(paste0("lowConf_", top_label))
    }
    
    if (max_score >= conf_low && diff < diff_threshold) {
      return(paste0("lowConf_ambiguous:", top_label, ":", second_label))
    }
    
    return("unknown")
  })
}

#' Run CellTypist annotation
#' @param seurat_obj Seurat object
#' @param model CellTypist model name or path
#' @param params List of annotation parameters
#' @return List with logits, softmax, and labels
run_celltypist_annotation <- function(
  seurat_obj,
  model = NULL, # for Immune_All_Low.pkl
  params = ANNOT_PARAMS()
) {
  
  message("\n==== Running CellTypist annotation ====")
  
  mtx <- seurat_obj@assays$RNA@counts
  
  message("Predicting cell labels...")
  lr_output <- runCelltypist(cnts = mtx)
  
  message("Computing stable softmax...")
  softmax_p <- stable_softmax(lr_output$logitMat)
  rownames(softmax_p) <- rownames(lr_output$logitMat)
  
  message("Assigning labels...")
  softmax_label <- assign_labels_from_probs(
    softmax_p,
    conf_high = params$celltypist_conf_high,
    conf_low = params$celltypist_conf_low,
    diff_threshold = params$celltypist_diff_threshold
  )
  
  # Calculate confidence metrics
  max_probs <- apply(softmax_p, 1, max)
  entropy <- apply(softmax_p, 1, function(p) {
    p <- p[p > 0]  # avoid log(0)
    -sum(p * log(p))
  })
  
  list(
    logits = lr_output$logitMat,
    softmax = softmax_p,
    predicted_labels = lr_output$labMat$predicted_labels,
    softmax_labels = softmax_label,
    max_prob = max_probs,
    entropy = entropy
  )
}

## Marker gene lists -----------------------------------------------------------
pbmc_markers_comprehensive <- list(
  # Hematopoietic stem/progenitor
  HSC_MPP = c("CD34","CD38","THY1","CD90","MLLT3","CD133","CD45RA"),
  Progenitors = c("CD105","CD44","CD73","PTPRC","CD29","STRO1","ALDH1A1","KDR","CD117"),
  
  # B cell lineage
  PreProB = c("IGLL1","CD79B","VPREB1","EBF1","CD24","CD38","CD99","TCL1A","ARPP21"),
  ProB = c("DNTT","CD79B","RAG1","EBF1","VPREB3","IGLL5","PAX5"),
  PreB = c("MME","CD79A","TCL1A","IGLL5","CD24","CD38"),
  Naive_B = c("CD19","CD27","IgD","MS4A1","CD38","CD24","PAX5","JCHAIN"),
  Memory_B = c("CD19","CD27","IgD","MS4A1","CD38","CD24","JCHAIN"),
  Plasma_B = c("CD27","SDC1","XBP1","JCHAIN","MS4A1","CD38"),
  
  # T cell lineage
  Naive_CD4_T = c("CD3D","CD3E","CD3G","CD4","CD45RA","CCR7","CD27","CD127"),
  Naive_CD8_T = c("CD3D","CD3E","CD3G","CD8A","CD8B","CD45RA","CCR7","CD27","CD127"),
  Memory_CD4_T = c("CD3D","CD3E","CD3G","CD4","CD45RO","CD27","CCR7"),
  Memory_CD8_T = c("CD3D","CD3E","CD3G","CD8A","CD8B","CD45RO","CD27","CCR7"),
  Effector_CD4_T = c("CD3D","CD3E","CD3G","CD4","PRF1","GZMA","GZMB","NKG7"),
  Effector_CD8_T = c("CD3D","CD3E","CD3G","CD8A","PRF1","GZMA","GZMB","NKG7"),
  Regulatory_T = c("FOXP3","IL2RA","CD127"),
  GammaDelta_T = c("TRDV2","TRGV9","TRGC1"),
  MAIT_T = c("SLC4A10","TRAV1-2","KLRB1"),
  
  # NK and NKT
  NK = c("NKG7","KLRD1","CD56","GNLY","FCGR3A"),
  NKT_CD8 = c("CD8","CD56","NKG7","GZMB","GNLY"),
  NKT_CD4 = c("CD4","CD56","NKG7","GZMB","GNLY"),
  
  # Monocytes and macrophages
  Classical_Mono = c("CD14","ITGAM","ITGB2","HLA-DR","LYZ","S100A8","S100A9"),
  NonClassical_Mono = c("CD16","FCGR3A","CD14","S100A8","S100A9","C1QC","CST3"),
  Intermediate_Mono = c("CD14","CD16","IL1B","S100A8","S100A9","C1QC"),
  Macrophages = c("CD68","CD163","MSR1","APOE","CD206","CD80","CD86"),
  
  # Dendritic cells
  cDC1 = c("CLEC9A","ITGAX","CD1C"),
  cDC2 = c("CD1C","ITGAX","IL3RA"),
  pDC = c("CLEC4C","IL3RA","ITGAX"),
  Myeloid_DC = c("ITGAX","CD1C","CD83","CD86","HLA-DRA","HLA-DRB1","IL3RA"),
  
  # Granulocytes
  Neutrophils = c("CXCR2","CSF3R","FCGR3B","MPO","ELANE","PRTN3","S100A8","S100A9"),
  Eosinophils = c("SIGLEC8","PRG2","EPX","IL5RA","CD193","CLC"),
  Basophils = c("CD123","CD203c","CD63","HDC","MS4A2","FCER1A"),
  Mast_cells = c("KIT","TPSAB1","CD203c","FCER1A"),
  
  # Megakaryocytes / Platelets
  Megakaryocytes = c("PF4","ITGA2B","CD41","CD42b","PPBP","CXCR4"),
  Platelets = c("PPBP","PF4","ITGA2B","CD41","CD42b"),
  
  # Erythroid
  Early_Erythroid = c("GATA1","KLF1","ALAS2"),
  Mid_Erythroid = c("HBA1","HBA2","BPGM","HBB"),
  Erythroid_precursor = c("GYPA","HBD","TFRC","CD36")
)


pbmc_markers_short <- list(
  HSC_MPP = c("CD34","CD38","THY1"),
  ProB = c("DNTT","CD79B","EBF1"),
  PreB = c("MME","CD79A","TCL1A"),
  Naive_B = c("CD19","MS4A1","CD27"),
  Memory_B = c("CD19","CD27"),
  Plasma_B = c("CD27","SDC1","XBP1"),
  
  Naive_CD4_T = c("CD3D","CD4","CCR7","CD45RA"),
  Naive_CD8_T = c("CD3D","CD8A","CCR7","CD45RA"),
  Memory_CD4_T = c("CD3D","CD4","CD45RO","CD27"),
  Memory_CD8_T = c("CD3D","CD8A","CD45RO","CD27"),
  Effector_T = c("PRF1","GZMA","GZMB","NKG7"),
  Regulatory_T = c("FOXP3","IL2RA"),
  GammaDelta_T = c("TRDV2","TRGV9"),
  MAIT_T = c("SLC4A10","TRAV1-2"),
  
  NK = c("NKG7","KLRD1","GNLY"),
  
  Classical_Mono = c("CD14","LYZ","S100A8"),
  NonClassical_Mono = c("FCGR3A","CD16","C1QC"),
  Macrophages = c("CD68","CD163","MSR1"),
  
  cDC1 = c("CLEC9A","ITGAX"),
  cDC2 = c("CD1C","ITGAX"),
  pDC = c("CLEC4C","IL3RA"),
  
  Neutrophils = c("CXCR2","FCGR3B","MPO"),
  Eosinophils = c("SIGLEC8","PRG2"),
  Basophils = c("CD123","TPSAB1"),
  Mast_cells = c("KIT","TPSAB1"),
  
  Megakaryocytes = c("PF4","ITGA2B","PPBP"),
  Platelets = c("PPBP","PF4"),
  
  Early_Erythroid = c("GATA1","KLF1","ALAS2"),
  Mid_Erythroid = c("HBA1","HBB")
)

# Comprehensive PBMC/BM markers
HEME_MARKERS <- list(
  # ============ Stem/Progenitor Cells ============
  HSC = c("CD34", "THY1", "KIT", "CD90"),
  MPP = c("CD34", "CD38", "FLT3"),
  CMP = c("CD34", "CD38", "IL3RA"),
  GMP = c("CD34", "CD38", "IL3RA", "CSF3R"),
  MEP = c("CD34", "CD38", "IL3RA", "CD36"),
  
  # ============ B Cell Lineage ============
  ProB = c("DNTT", "CD79B", "EBF1", "IGLL1"),
  PreB = c("MME", "CD79A", "TCL1A", "VPREB1"),
  Naive_B = c("CD19", "MS4A1", "TCL1A", "IGHD"),
  Memory_B = c("CD19", "MS4A1", "CD27", "IGHG1"),
  Plasma_B = c("CD27", "SDC1", "XBP1", "MZB1", "JCHAIN"),
  
  # ============ T Cell Lineage ============
  Naive_CD4_T = c("CD3D", "CD4", "CCR7", "SELL", "LEF1"),
  Naive_CD8_T = c("CD3D", "CD8A", "CCR7", "SELL", "LEF1"),
  Memory_CD4_T = c("CD3D", "CD4", "IL7R", "CD27"),
  Memory_CD8_T = c("CD3D", "CD8A", "GZMK", "CD27"),
  Effector_CD8_T = c("CD8A", "GZMB", "PRF1", "FGFBP2"),
  Regulatory_T = c("FOXP3", "IL2RA", "IKZF2", "CTLA4"),
  GammaDelta_T = c("TRDV2", "TRGV9", "TRDC"),
  MAIT_T = c("SLC4A10", "TRAV1-2", "KLRB1"),
  
  # ============ NK Cells ============
  NK = c("NKG7", "KLRD1", "GNLY", "NCAM1"),
  NK_CD56bright = c("NCAM1", "FCGR3A", "XCL1"),
  NK_CD56dim = c("FCGR3A", "PRF1", "GZMB"),
  
  # ============ Myeloid Lineage ============
  Classical_Mono = c("CD14", "LYZ", "S100A8", "S100A9", "FCN1"),
  Intermediate_Mono = c("CD14", "FCGR3A", "HLA-DRA"),
  NonClassical_Mono = c("FCGR3A", "MS4A7", "CDKN1C"),
  Macrophages = c("CD68", "CD163", "MSR1", "MRC1"),
  
  # ============ Dendritic Cells ============
  cDC1 = c("CLEC9A", "XCR1", "CADM1"),
  cDC2 = c("CD1C", "FCER1A", "CLEC10A"),
  pDC = c("CLEC4C", "IL3RA", "GZMB", "TCF4"),
  
  # ============ Granulocytes ============
  Neutrophils = c("CXCR2", "FCGR3B", "MPO", "CSF3R"),
  Eosinophils = c("SIGLEC8", "PRG2", "CLC", "EPX"),
  Basophils = c("IL3RA", "TPSAB1", "HDC", "MS4A2"),
  Mast_cells = c("KIT", "TPSAB1", "CPA3", "HDC"),
  
  # ============ Megakaryocyte/Platelet ============
  Megakaryocytes = c("PF4", "ITGA2B", "PPBP", "GP9"),
  Platelets = c("PPBP", "PF4", "CAVIN2"),
  
  # ============ Erythroid ============
  Early_Erythroid = c("GATA1", "KLF1", "ALAS2", "TFRC"),
  Late_Erythroid = c("HBA1", "HBA2", "HBB", "GYPA"),
  
  # ============ Leukemia/Blast Markers ============
  Myeloblast = c("CD34", "CD117", "MPO", "CD33"),
  Lymphoblast = c("CD34", "CD10", "TDT", "CD19"),
  Immature = c("CD34", "CD38", "CD123"),
  
  # ============ Proliferating ============
  Cycling = c("MKI67", "TOP2A", "PCNA", "STMN1")
)

#' Add module scores to Seurat object
#' @param seurat_obj Seurat object
#' @param marker_list Named list of marker genes
#' @return Seurat object with module scores added
add_marker_module_scores <- function(seurat_obj, marker_list) {
  
  message("\n==== Adding module scores ====")
  
  for (celltype in names(marker_list)) {
    safe_name <- gsub("[^[:alnum:]_]", "_", celltype)
    seurat_obj <- AddModuleScore(
      seurat_obj,
      features = list(marker_list[[celltype]]),
      name = paste0(safe_name, "_score"),
      assay = "RNA"
    )
  }
  
  return(seurat_obj)
}


#' Integrate CellTypist and marker-based annotations
#' @param seurat_obj Seurat object with both annotations
#' @param celltypist_col Column name for CellTypist predictions
#' @param marker_col Column name for marker predictions
#' @param ambiguity_col Column name for ambiguity scores
#' @param ambiguity_threshold Threshold for calling ambiguous
#' @return Vector of consensus labels
integrate_annotations <- function(
  seurat_obj,
  celltypist_col = "celltypist_label",
  marker_col = "marker_celltype",
  ambiguity_col = "marker_ambiguity",
  ambiguity_threshold = 0.2
) {
  
  message("\n==== Integrating annotations ====")
  
  md <- seurat_obj@meta.data
  
  consensus <- sapply(1:nrow(md), function(i) {
    ct_label <- md[[celltypist_col]][i]
    mk_label <- md[[marker_col]][i]
    ambig <- md[[ambiguity_col]][i]
    
    # Parse CellTypist label
    ct_is_ambig <- grepl("ambiguous|lowConf|unknown", ct_label)
    ct_clean <- gsub("^lowConf_|^ambiguous:", "", ct_label)
    ct_clean <- strsplit(ct_clean, ":")[[1]][1]
    
    # Decision logic
    if (!ct_is_ambig && !is.na(ambig) && ambig > ambiguity_threshold) {
      # Both agree and confident
      if (tolower(ct_clean) == tolower(mk_label)) {
        return(paste0("consensus:", ct_clean))
      } else {
        return(paste0("partial:", ct_clean, "|", mk_label))
      }
    } else if (!ct_is_ambig) {
      return(paste0("celltypist:", ct_clean))
    } else if (!is.na(ambig) && ambig > ambiguity_threshold) {
      return(paste0("marker:", mk_label))
    } else {
      return("ambiguous")
    }
  })
  
  return(consensus)
}

#' Log annotation run information
#' @param base_dir Base output directory
#' @param version_dir Version directory
#' @param datasets Datasets being processed
#' @param params Parameters used
#' @param notes Optional notes about this run
log_annotation_run <- function(
  base_dir,
  version_dir,
  datasets,
  params = ANNOT_PARAMS,
  notes = NULL
) {
  
  log_file <- file.path(base_dir, "run_log.txt")
  
  log_entry <- c(
    sprintf("\n========================================"),
    sprintf("Annotation Run: %s", Sys.time()),
    sprintf("Version Directory: %s", basename(version_dir)),
    sprintf("========================================"),
    sprintf("Datasets: %s", paste(datasets, collapse = ", ")),
    sprintf("\nParameters:"),
    sprintf("  CellTypist model: %s", params$celltypist_model),
    sprintf("  High confidence threshold: %.2f", params$celltypist_conf_high),
    sprintf("  Low confidence threshold: %.2f", params$celltypist_conf_low),
    sprintf("  Min pct for markers: %.2f", params$min_pct),
    sprintf("  Ambiguity threshold: %.2f", params$ambiguity_threshold)
  )
  
  if (!is.null(notes)) {
    log_entry <- c(log_entry, sprintf("\nNotes: %s", notes))
  }
  
  log_entry <- c(log_entry, sprintf("========================================\n"))
  
  # Append to log file
  cat(paste(log_entry, collapse = "\n"), file = log_file, append = TRUE)
  
  message(sprintf("Run logged to: %s", log_file))
}

##----------------------------##
##   Visualization Functions  ####
##----------------------------##

#' Plot annotation QC
#' @param seurat_obj Seurat object
#' @param output_prefix Output file prefix
#' @param celltypist_result CellTypist results list
#' @param marker_result Marker scoring results list
#' @import patchwork
plot_annotation_qc <- function(
  seurat_obj,
  output_prefix,
  celltypist_result,
  marker_result,
  params
) {
  message("\n==== Generating QC plots ====")
  
  pdf(paste0(output_prefix, "_annotation_qc.pdf"), width = 25, height = 25)
  
  # ---- Page 1: CellTypist confidence ----
  p1 <- FeaturePlot(seurat_obj, features = "celltypist_max_prob",
                    cols = c("lightgrey", "blue")) +
    ggtitle("CellTypist Max Probability")
  
  p2 <- FeaturePlot(seurat_obj, features = "celltypist_entropy",
                    cols = c("blue", "lightgrey")) +
    ggtitle("CellTypist Entropy (uncertainty)")
  
  p3 <- DimPlot(seurat_obj, group.by = "celltypist_label",
                label = TRUE, repel = TRUE, label.size = 3) +
    NoLegend() + ggtitle("CellTypist Labels")
  
  p4 <- DimPlot(seurat_obj, group.by = "celltypist_softmax_label",
                label = TRUE, repel = TRUE, label.size = 3) +
    NoLegend() + ggtitle("CellTypist Softmax Labels")
  
  print((p1 | p2) / (p3 | p4))
  
  # ---- Page 2: Marker scoring ----
  p5 <- FeaturePlot(seurat_obj, features = "marker_ambiguity",
                    cols = c("lightgrey", "darkgreen")) +
    ggtitle("Marker Score Ambiguity")
  
  p6 <- DimPlot(seurat_obj, group.by = "marker_celltype",
                label = TRUE, repel = TRUE, label.size = 3) +
    NoLegend() + ggtitle("Marker-based Cell Types")
  
  p7 <- DimPlot(seurat_obj, group.by = "seurat_clusters",
                label = TRUE, label.size = 3) +
    NoLegend() + ggtitle("Seurat Clusters")
  
  print((p5 | p6) / p7)
  
  # ---- Page 3: Consensus ----
  p8 <- DimPlot(seurat_obj, group.by = "consensus_celltype",
                label = TRUE, repel = TRUE, label.size = 2.5) +
    ggtitle("Consensus Annotation") + NoLegend()
  
  print(p8)
  
  # ---- Page 4: CellTypist logits heatmap ----
  Idents(seurat_obj) <- "seurat_clusters"
  
  col_fun <- circlize::colorRamp2(
    c(-2, 0, 2),
    c("blue", "white", "red")
  )
  
  # Sample cells if too many (for visualization)
  max_cells_plot <- 5000
  if (ncol(seurat_obj) > max_cells_plot) {
    set.seed(123)
    cells_plot <- sample(colnames(seurat_obj), max_cells_plot)
    logit_subset <- celltypist_result$logits[cells_plot, ]
    cluster_subset <- seurat_obj$seurat_clusters[cells_plot]
  } else {
    logit_subset <- celltypist_result$logits
    cluster_subset <- seurat_obj$seurat_clusters[
      match(rownames(celltypist_result$logits), colnames(seurat_obj))
      ]
  }
  
  ht <- Heatmap(
    t(logit_subset),
    name = "Logit",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_title = "Cells",
    row_title = "Cell Types (CellTypist)",
    heatmap_legend_param = list(direction = "horizontal"),
    column_split = cluster_subset,
    column_title_gp = gpar(fontsize = 8)
  )
  
  draw(ht, heatmap_legend_side = "bottom")
  
  # ---- Page 5: Marker score heatmap ----
  ht2 <- Heatmap(
    t(marker_result$scores),
    name = "Score",
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    column_title = "Clusters",
    row_title = "Cell Types (Markers)",
    heatmap_legend_param = list(direction = "horizontal")
  )
  
  draw(ht2, heatmap_legend_side = "bottom")
  
  # ---- Page 6: Dotplot of key markers ----
  Idents(seurat_obj) <- "consensus_celltype"
  
  top_markers <- unique(unlist(HEME_MARKERS))
  top_markers <- intersect(top_markers, rownames(seurat_obj))
  
  if (length(top_markers) > 0) {
    p9 <- DotPlot(
      seurat_obj,
      features = top_markers,
      dot.scale = params$dotplot_dot_scale
    ) +
      RotatedAxis() +
      scale_color_viridis_c(option = "C") +
      ggtitle("Key Marker Expression by Consensus Cell Type") +
      theme(axis.text.x = element_text(size = 8))
    
    print(p9)
  }
  
  # ---- Page 7: Dotplot by clusters ----
  Idents(seurat_obj) <- "seurat_clusters"
  
  if (length(top_markers) > 0) {
    p10 <- DotPlot(
      seurat_obj,
      features = top_markers,
      dot.scale = params$dotplot_dot_scale
    ) +
      RotatedAxis() +
      scale_color_viridis_c(option = "C") +
      ggtitle("Key Marker Expression by Cluster") +
      theme(axis.text.x = element_text(size = 8))
    
    print(p10)
  }
  
  dev.off()
  
  message(sprintf("QC plots saved to: %s_annotation_qc.pdf", output_prefix))
}