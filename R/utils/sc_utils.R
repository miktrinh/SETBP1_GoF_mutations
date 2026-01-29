library(Seurat)
s.genes = cc.genes.updated.2019$s.genes
g2m.genes = cc.genes.updated.2019$g2m.genes



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

# annotation helpers -----------------------------------------------------------
score_clusters_markers <- function(
  seurat_obj,
  marker_list,
  assay = "RNA",
  slot = "data",
  min_pct = 0.1,
  scale_across_clusters = TRUE,
  weight_pct = 1,
  weight_expr = 1
) {
  Idents(seurat_obj) <- "seurat_clusters"
  
  expr <- GetAssayData(seurat_obj, assay = assay, slot = slot)
  clusters <- levels(Idents(seurat_obj))
  
  # genes actually present
  genes <- intersect(unique(unlist(marker_list)), rownames(expr))
  if (length(genes) == 0)
    stop("None of the marker genes found in expression matrix")
  
  # initialize
  avg_exp <- matrix(0, nrow = length(genes), ncol = length(clusters))
  pct_exp <- matrix(0, nrow = length(genes), ncol = length(clusters))
  rownames(avg_exp) <- rownames(pct_exp) <- genes
  colnames(avg_exp) <- colnames(pct_exp) <- clusters
  
  # compute cluster-level stats
  for (cl in clusters) {
    cells <- WhichCells(seurat_obj, idents = cl)
    mat <- expr[genes, cells, drop = FALSE]
    
    avg_exp[, cl] <- Matrix::rowMeans(mat)
    pct_exp[, cl] <- Matrix::rowMeans(mat > 0)
  }
  
  # remove uninformative genes
  keep <- matrixStats::rowMaxs(pct_exp) >= min_pct
  avg_exp <- avg_exp[keep, , drop = FALSE]
  pct_exp <- pct_exp[keep, , drop = FALSE]
  
  # specificity normalization
  if (scale_across_clusters) {
    avg_exp <- t(scale(t(avg_exp)))
    avg_exp[is.na(avg_exp)] <- 0
  }
  
  # combined gene-level score
  gene_score <- (weight_expr * avg_exp) * (weight_pct * pct_exp)
  
  # aggregate per cell type
  celltype_scores <- sapply(marker_list, function(g) {
    g <- intersect(g, rownames(gene_score))
    if (length(g) == 0) return(rep(0, ncol(gene_score)))
    colMeans(gene_score[g, , drop = FALSE])
  })
  
  rownames(celltype_scores) <- clusters
  
  # ambiguity metric
  score_df <- as.data.frame(celltype_scores)
  ambiguity <- apply(score_df, 1, function(x) {
    sx <- sort(x, decreasing = TRUE)
    if (length(sx) < 2) return(NA)
    sx[1] - sx[2]
  })
  
  list(
    scores = score_df,
    predicted = apply(score_df, 1, function(x) names(which.max(x))),
    ambiguity = ambiguity
  )
}

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
