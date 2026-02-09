# Revision 2602 Figures 

# Libraries --------------------------------------------------------------------
library(ggplot2)
source('R/utils/misc.R')

# Global variables -------------------------------------------------------------
plotDir = "Results/Revision_2602"
if(!dir.exists(plotDir)){
  dir.create(plotDir,recursive = T)
}

ANNOT_COLS <- c(
  classical_monocytes       = "#4E79A7",
  non_classical_monocytes   = "#A0CBE8",
  macrophages               = "#59A14F",
  dc1                       = "#8CD17D",
  dc2                       = "#B6992D",
  pdc                       = "#F1CE63",
  
  t_cd4                     = "#E15759",
  t_cd8                     = "#ECAFAD",
  gd_t_cells                = "#D37295",
  mait_cells                = "#741b47",
  nkt_cells                 = "#7A6F9B",
  
  nk                        = "#79706E",
  ilc                       = "#BAB0AC",
  
  naive_b                   = "#9FCCC9",
  memory_b                  = "#4BA3A2",
  plasma_b                  = "#306918",
  
  megakaryocyte             = "#B07AA1",
  mds                       = "#9C755F"
)

ANNOT_COLS_BROAD <- c(
  "B cells"        = "#9FCCC9",
  
  "CD16+ Mono"     = "#4E79A7",
  "Mono/Mac"       = "#A0CBE8",
  "Myeloid"        = "#A0CBE8",
  
  "DC1"            = "#8CD17D",
  "DC2"            = "#B6992D",
  "pDC"            = "#F1CE63",
  
  "T cells"        = "#e6a19e",
  "NK cells"       = "#79706E",
  
  "NK_T"           = "#741b47",
  
  "Plasma cells"   = "#306918",
  
  "MegK"           = "#B07AA1",
  
  "MDS"            = "#bf6a02"
  
)

LINEAGE_COLS <- c(
  'T'            = "#e6a19e",
  'NK/ILC'       = "#79706E",
  B              = "#9FCCC9",
  DC             = "#B6992D",
  Myeloid        = "#4E79A7",
  MDS            = "#9C755F",
  Megakaryocyte  = "#B07AA1"
)

donorID_col=c('PD54849' = pal34H[9],
              'PD66174' = pal34H[27],
              # SGS
              'PD60412' = brewer.pal(n=10,'Paired')[7],
              'PD66162' = pal34H[18],
              'PD66163' = pal34H[17],
              # Normal
              'PD66164' = grey(0.4),
              'GSM5160432' = grey(0.85),
              'GSM5160434' = grey(0.7),
              'GSM5160435' = grey(0.55),
              # MDS
              'PD61858' = pal34H[3],
              'PD61857' = pal37H[14]
)
SETBP1_mutation_col = c('Wild type SETBP1' = colAlpha(grey(0.6),0.3),
                        'Germline hotspot SETBP1 mutation' = colAlpha('#cc5d0f',0.55),
                        'Germline hotspot-adjacent SETBP1 mutation' = colAlpha('#1b4891',0.55),
                        "Somatic hotspot SETBP1 mutation" = colAlpha('orange',0.55))
scales::show_col(SETBP1_mutation_col)
combined_srat_clean$setbp1_mutation_status[combined_srat_clean$donor_id %in% c('L061') & 
                                             combined_srat_clean$annot != 'mds'] = 'Wild type SETBP1'

umap_theme <- theme_void() +
  theme(
    legend.position = 'none',
    plot.title   = element_text(color = "black", hjust = 0.5)
  )


# Combined Integrated UMAP -----------------------------------------------------
integrated_umap_fp = 'Results/03_scRNAseq_annotation/SGS_inhouse_MDS_L061_L067_annot_clean_HARM_2602_mdat.csv'
sgs_inhouse_umap_fp = 'Results/03_scRNAseq_annotation/SGS_inhouse_annot_clean_2601_mdat.csv'
integrated_sgs_inhouse_umap_fp = 'Results/03_scRNAseq_annotation/SGS_inhouse_annot_clean_HARM_2602_mdat.csv'
mds_umap_fp = 'Results/03_scRNAseq_annotation/MDS_L061_L067_annot_clean_2602_mdat.csv'

plot_umap = function(plotDir,
                     dataset = c('sgs_inhouse','mds','combined'), 
                     colour_by_variables = c('annot','broad_annot','lineage','donorID','setbp1_status'), 
                     is_integrated=FALSE,
                     show_legend = FALSE,
                     alpha=1,
                     down_sample_fraction = 1,
                     width = 4.1,
                     height = 3.7){
  checkmate::qassert(plotDir,'S1')
  checkmate::qassert(colour_by_variables,'S+')
  checkmate::qassert(is_integrated,'B1')
  checkmate::qassert(show_legend,'B1')
  checkmate::qassert(alpha,'N1[0,1]')
  checkmate::qassert(down_sample_fraction,'N1[0,1]')
  checkmate::qassert(width,'N1')
  checkmate::qassert(height,'N1')
  dataset <- match.arg(dataset)
  
  mdat_fp = switch(dataset,
                   'sgs_inhouse' = ifelse(is_integrated,integrated_sgs_inhouse_umap_fp,sgs_inhouse_umap_fp),
                   'mds' = mds_umap_fp,
                   'combined' = integrated_umap_fp)
  selected_columns = c('cellID','annot','final_annot_broad','broad_annot','wgs_id','setbp1_mutation_status','UMAP_1','UMAP_2')
  
  mdat = read.csv(mdat_fp,row.names = 1) %>%
    dplyr::select(selected_columns) %>% 
    dplyr::mutate(donorID = substr(wgs_id,1,nchar(wgs_id)-1))
  if(dataset == 'combined'){
    mdat$final_annot_broad[mdat$final_annot_broad %in% c('CD16+ Mono','Mono/Mac',
                                                         'DC1','DC2','pDC')] = "Myeloid"
    mdat$final_annot_broad[mdat$final_annot_broad %in% c('NK_T')] = "T cells"
  }
  if(down_sample_fraction < 1 ){
    dd = mdat[sample(1:nrow(mdat),nrow(mdat)*down_sample_fraction),]
  }else{
    dd = mdat[sample(1:nrow(mdat),nrow(mdat)),]
  }
  
  
  
  for(colour_by in colour_by_variables){
    if(colour_by == 'annot' & dataset == 'combined'){
      dd2 = rbind(dd[!dd$annot %in% c('nkt_cells','ilc'),],
                  dd[dd$annot %in% c('nkt_cells','ilc'),]
                 )
      checkmate::assert_true(nrow(dd2) == nrow(dd))
      dd = dd2
    }
    cols = switch(colour_by,
                  'annot' = ANNOT_COLS[as.character(dd$annot)],
                  'broad_annot' = ANNOT_COLS_BROAD[as.character(dd$final_annot_broad)],
                  'lineage' = LINEAGE_COLS[as.character(dd$broad_annot)],
                  'donorID' = donorID_col[as.character(dd$donorID)],
                  'setbp1_status' = SETBP1_mutation_col[as.character(dd$setbp1_mutation_status)])
    
    cols = colAlpha(cols,rep(alpha,length(cols)))
    plot_prefix = ifelse(is_integrated, paste0(dataset,'_integrated_UMAP_',colour_by),
                         paste0(dataset,'_UMAP_',colour_by))
    
    plotFun_annot = function(noFrame=FALSE,noPlot=FALSE){
      par(mar=c(0,0,0.1,0))
      
      plot(dd$UMAP_1,dd$UMAP_2,
           las=1,
           type='n',
           cex.main = 0.85,xaxt='n',yaxt='n',
           xlab='',ylab='',
           frame.plot=F)
      
      if(!noPlot){
        points(dd$UMAP_1,dd$UMAP_2,
               col = cols,
               pch = 19,
               cex=0.07)
      }
    }
    
    saveFig(file.path(plotDir,plot_prefix),plotFun_annot,rawData=mdat,
            width = width,height = height,res = 500)
  }
}
  
  




plot_umap(plotDir=plotDir,
          dataset = 'sgs_inhouse', 
          colour_by_variables = c('annot','broad_annot','lineage','donorID','setbp1_status'), 
          #colour_by_variables = c('broad_annot'), 
          is_integrated=TRUE,
          show_legend = FALSE,
          alpha=0.7,
          down_sample_fraction = 1,
          width = 4.1,
          height = 3.7)

plot_umap(plotDir=plotDir,
          dataset = c('mds'), 
          colour_by_variables = c('annot','broad_annot','donorID','setbp1_status'), 
          is_integrated=FALSE,
          show_legend = FALSE,
          alpha=1,
          down_sample_fraction = 1,
          width = 4.1,
          height = 3.7)

plot_umap(plotDir=plotDir,
          dataset = 'combined', 
          colour_by_variables = c('annot','broad_annot','lineage','donorID','setbp1_status'), 
          #colour_by_variables = c('broad_annot'), 
          is_integrated=TRUE,
          show_legend = FALSE,
          alpha=0.5,
          down_sample_fraction = 1,
          width = 4.1,
          height = 3.7)





plotFun_annot_seurat = function(noFrame=FALSE,noPlot=FALSE){
    p <- Seurat::DimPlot(combined_srat_clean,group.by = 'annot',cols = ANNOT_COLS,
                    label = T,repel = T,label.box = T,label.size = 2)+
      theme_void()
    
    if(!show_legend){
      p <- p + NoLegend()
    }
    
    print(p)
  }
  
  saveFig(file.path(plotDir,'integrated_UMAP_annot_seuratstyle'),plotFun_annot_seurat,
          rawData=mdat,width = 4.5,height = 4.5,res = 500)
  
  plotFun_sgs_seurat = function(noFrame=FALSE,noPlot=FALSE){
    p <- Seurat::DimPlot(sgs_inhouse_annot_clean,group.by = 'annot',cols = ANNOT_COLS,
                         label = T,repel = T,label.box = T,label.size = 2) + 
      NoLegend()+
      umap_theme
    print(p)
  }
  
  saveFig(file.path(plotDir,'sgs_inhouse_UMAP_annot_seuratstyle'),plotFun_sgs_seurat,
          rawData=mdat,width = 4.5,height = 4.5,res = 500)
  saveFig(file.path(plotDir,'mds_UMAP_annot_seuratstyle'),plotFun_sgs_seurat,
          rawData=mdat,width = 4.5,height = 4.5,res = 500)
  
  
  
  
  plotFun_donorID_seurat = function(noFrame=FALSE,noPlot=FALSE){
    p <- Seurat::DimPlot(combined_srat_clean,group.by = 'donorID',cols = donorID_col,
                         label = F,repel = T,label.box = T,label.size = 2) + 
      theme_void() 
    print(p)
  }
  
  saveFig(file.path(plotDir,'integrated_UMAP_donorID_seuratstyle'),plotFun_donorID_seurat,
          rawData=mdat,width = 5,height = 4.5,res = 500)
  
  
  DimPlot(combined_srat_clean,cells.highlight = combined_srat_clean$cellID[combined_srat_clean$annot == 'mds' &
                                                                             combined_srat_clean$donor_id == 'L067'])
  plotFun_disease_seurat = function(noFrame=FALSE,noPlot=FALSE){
    p <- Seurat::DimPlot(combined_srat_clean,group.by = 'setbp1_mutation_status',cols = SETBP1_mutation_col,
                         label = F,repel = T,label.box = T,label.size = 2) + 
      theme_void() 
    print(p)
  }
  
  saveFig(file.path(plotDir,'integrated_UMAP_disease_seuratstyle'),plotFun_disease_seurat,
          rawData=mdat,width = 7,height = 4.5,res = 500)
  
  
  
## DotPlot ---------------------------------------------------------------------
markers = (c('PTPRC', # Immune cells
             # MDS
             'SETBP1','CD34','HMGA2','MEG3','RBPMS', 
             'MLLT3','PRSS57', 'GATA2','GATA1',
             # Megakaryocyte
             'ITGA2B', 'PPBP','PF4',
             # Monocytes
             'CD14',"LYZ","S100A8","FCN1",'CD68','HMOX1','ITGAX',#Macrophages
             'FCGR3A',# Monocytes
             
             
             'IRF8',	#DC.precursor
             'FLT3', #DCs
             'CLEC9A',#DC1
             'CLEC10A','CD1C', # DC2 
             'IL3RA', # DC2
             'CLEC4C', #pDC
             
             # B-cells
             'MS4A1',
             'CD79A','CD19',"IGHD","IL4R",# naive_b
             'TCL1A',# naive_b_cell = CD19+ CD27-
             # memory_b_cell = CD19+ and CD27 + 
             "IGHA1",# pre-b-cells
             
             "JCHAIN",'CD27', # Plasma cells
             
             'CD3D','CD8A', #Early.lymphoid_T.lymphocyte
             #'TRDV2','TRGV9', # gamma delta T-cell
             'PRF1','GZMA', # effector T cells
             "SLC4A10", "TRAV1-2", #MAIT
             "TRDV2", "TRGV9",  # Gamma-delta
             'ZBTB16','NKG7','KLRD1', #NK
             'IL7R',"LTB","KIT" # ILC
             
))

combined_srat_fp = 'Results/03_scRNAseq_annotation/SGS_inhouse_MDS_L061_L067_annot_clean_HARM_2601.RDS'
srat = combined_srat_clean
supFig1_SETBP1.MDS_dotPlot = function(){
  library(viridis)
  
  broad_annot_level = c('MDS','Ery','MegK',
                  'Mono/Mac','CD16+ Mono','DC1','DC2','pDC',
                  'B cells','Plasma cells','T cells','NK cells')
  detailed_annot_level = c('mds_PD61858','mds_PD61857','megakaryocyte',
                           'classical_monocytes','non_classical_monocytes',
                           'dc1','dc2','pdc','naive_b','memory_b','plasma_b',
                           't_cd4','t_cd8','mait_cells','gd_t_cells',
                           'nkt_cells','nk','ilc')
  # Corresponding paper-ready names
  paper_names <- c(
    "MDS_PD61858","MDS_PD61857", "Megakaryocyte",
    "Classical Mono", "Non-Classical Mono",
    "DC1", "DC2", "pDC", "Naive B", "Memory B", "Plasma B",
    "CD4 T", "CD8 T", "MAIT", "γδ T",
    "NKT", "NK", "ILC"
  )
  
  # Create a mapping table
  cell_type_mapping <- data.frame(
    code_style = detailed_annot_level,
    paper_style = paper_names,
    stringsAsFactors = FALSE
  )
  
  # Import the seurat object
  srat = readRDS(combined_srat_fp)
  srat$annot[srat$annot == 'mds'] = paste0('mds_',srat$donorID[srat$annot == 'mds'])
  
  checkmate::assert_true(all(srat$final_annot_broad %in% broad_annot_level))
  checkmate::assert_true(all(srat$annot %in% detailed_annot_level))
  srat$final_annot_broad = factor(srat$final_annot_broad,broad_annot_level)
  srat$annot = factor(srat$annot,detailed_annot_level)
  srat$annot2 = cell_type_mapping$paper_style[match(srat$annot,cell_type_mapping$code_style)]
  srat$annot2 = factor(srat$annot2,paper_names)
  
  direction = 'vertical'
  
  if(direction == 'vertical'){
    genes = rev(markers)
  }else{
    genes = markers
  }
  
  plotFun = function(noFrame=FALSE,noPlot=FALSE){
    Idents(srat) = srat$annot2
    p = DotPlot(srat,#idents = unique(srat$finalAnn_broad),
                #cols = c("#EBFFE5", "#244b05"),
                #cols = c(colAlpha('#F1F5FA',1),'#425580'),
                cols = c(colAlpha(grey(0.95),0.8),'black'),
                #cols = c(grey(0.99), grey(0.2)),
                #group.by = 'seurat_clusters',
                #idents = unique(srat$finalAnn[srat$finalAnn != 'others']),
                features = genes)+RotatedAxis() 
    
    
    if(direction == 'vertical'){
      p = p + coord_flip() + 
        scale_y_discrete(position = "right")+
        theme(axis.text.x = element_text(size=9,angle = 90,vjust = 1,hjust = 0),
              axis.text.y = element_text(size=9,face = "italic"),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top') + xlab('') + ylab('')
    }else if(direction == 'horizontal'){
      p = p +
        theme(axis.text.y = element_text(size=15),
              axis.text.x = element_text(size=8,angle = 90,vjust = 0.5,hjust = 1),
              legend.title = element_text(size=8),
              legend.text = element_text(size=8),
              legend.position = 'top') + xlab('') + ylab('')
    }
    print(p)
  }
  
  if(direction == 'vertical'){
    saveFig(file.path(plotDir,'SupFig1_SETBP1.MDS_Celltype_DotPlot_vertical'),plotFun,width = 4.7,height =11,res = 500)  
  }else{
    saveFig(file.path(plotDir,'SupFig1_SETBP1.MDS_Celltype_DotPlot_horizontal'),plotFun,width = 9.5,height = 4.6,res = 500)  
  }
  
}


# R1: NK-kB pathways in lymphoid lineage ---------------------------------------
library(enrichR)
library(zeallot)

listEnrichrSites()
setEnrichrSite("Enrichr") # Human genes
websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021",
         'HDSigDB_Human_2021','KEGG_2021_Human','MSigDB_Hallmark_2020','MSigDB_Oncogenic_Signatures',
         'TF_Perturbations_Followed_by_Expression','WikiPathways_2019_Human',
         'Cancer_Cell_Line_Encyclopedia','CCLE_Proteomics_2020','CellMarker_Augmented_2021',
         'Disease_Perturbations_from_GEO_down','Disease_Perturbations_from_GEO_up',
         'Drug_Perturbations_from_GEO_2014','Drug_Perturbations_from_GEO_down','Drug_Perturbations_from_GEO_up','DrugMatrix','IDG_Drug_Targets_2022',
         'Elsevier_Pathway_Collection','Enrichr_Submissions_TF-Gene_Coocurrence','Gene_Perturbations_from_GEO_down','Gene_Perturbations_from_GEO_up')

germline_sig = read.csv('Results/TableS2_germline.csv')
somatic_sig = read.csv('Results/TableS3_somatic.csv')
gene_module_list = list(
  g_lymph_up = germline_sig$Gene.symbol[germline_sig$Lineage == 'Lymphoid' & 
                                                 germline_sig$DE.Direction == "Up-regulated in germline SETBP1-mutated cells"],
  g_lymph_down = germline_sig$Gene.symbol[germline_sig$Lineage == 'Lymphoid' & 
                                          germline_sig$DE.Direction == "Down-regulated in germline SETBP1-mutated cells"],
  g_myeloid_up = germline_sig$Gene.symbol[germline_sig$Lineage == 'Myeloid' & 
                                          germline_sig$DE.Direction == "Up-regulated in germline SETBP1-mutated cells"],
  g_myeloid_down = germline_sig$Gene.symbol[germline_sig$Lineage == 'Myeloid' & 
                                            germline_sig$DE.Direction == "Down-regulated in germline SETBP1-mutated cells"],
  s_up = somatic_sig$Gene.symbol[somatic_sig$DE.Direction == 'Upregulated in hotspot-SETBP1-mutated aCML'],
  s_down = somatic_sig$Gene.symbol[somatic_sig$DE.Direction == 'Downregulated in hotspot-SETBP1-mutated aCML'],
  s_specific_up = somatic_sig$Gene.symbol[somatic_sig$DE.Direction == 'Upregulated in hotspot-SETBP1-mutated aCML' &
                                            somatic_sig$is_in_Shared_SETBP1_module == FALSE],
  s_myeloid_shared_up = somatic_sig$Gene.symbol[somatic_sig$DE.Direction == 'Upregulated in hotspot-SETBP1-mutated aCML' &
                                            somatic_sig$is_in_Shared_SETBP1_module == TRUE],
  g_myeloid_specific_up = germline_sig$Gene.symbol[germline_sig$Lineage == 'Myeloid' & 
                                                     germline_sig$DE.Direction == "Up-regulated in germline SETBP1-mutated cells" &
                                                     germline_sig$is_in_Shared_SETBP1_module == FALSE]
)
gene_module_list[['g_lymph_specific_up']] = germline_sig$Gene.symbol[germline_sig$Lineage == 'Lymphoid' & 
                                                 germline_sig$DE.Direction == "Up-regulated in germline SETBP1-mutated cells" &
                                                 !germline_sig$Gene.symbol %in% gene_module_list[['s_myeloid_shared_up']]]

# Run enrichR on germline-lymphoid gene module
c( g_lymph_up_enrichR,
   g_lymph_down_enrichR,
   g_myeloid_up_enrichR,
   g_myeloid_down_enrichR,
   s_up_enrichR,
   s_down_enrichR,
   s_specific_up_enrichR,
   s_myeloid_shared_up_enrichR,
   g_myeloid_specific_up_enrichR,
   g_lymph_specific_up_enrichR) %<-%
  lapply(gene_module_list,function(x) enrichr(x,dbs))

enrichr_res = g_lymph_specific_up_enrichR
plotEnrich(enrichr_res[[6]][enrichr_res[[6]]$Adjusted.P.value < 0.05,], showTerms = 50, numChar = 40, y = "Count", orderBy = "Adjusted.P.value")


allTerms_up = extract_enrichR.res(enrichR_result = enrichr_res,
                                  module_direction = 'up',
                                  pVal_cutoff = 0.1,
                                  min_overlapFrac = 0.04,
                                  db_toExtract = names(enrichr_res)[c(6)],
                                  nTerm = 15)

## Plot up-MSigDB_Hallmark_2020
module_type = 'SETBP1_germline_lymphoid_upregulated'
plotFun_enrichR_res = function(noFrame=FALSE,noPlot=FALSE){
  p3=ggplot(allTerms_up,aes(Combined.Score,Term))+
    geom_segment(aes(x=0,xend=(Combined.Score),y=Term,yend=Term),col=grey(0.7))+
    geom_point(aes(size = nGene,col=overlap2))+
    facet_grid(db~.,scales = 'free_y',space = 'free_y')+
    xlab('Combined Score')+ylab('') +
    ggtitle(module_type)+
    scale_color_gradient(low=grey(0.8),high = 'black')+
    theme_classic(base_size = 11)+
    scale_x_log10()+
    theme(panel.border = element_rect(fill=F,colour = 'black'),axis.line = element_blank(),
          strip.background=element_blank(),
          strip.text.y = element_text(size=10,colour = 'black'),
          axis.ticks = element_line(colour = 'black',linewidth = 0.2),
          axis.text.x = element_text(size = 10,angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
          axis.text.y = element_text(size=10),
          axis.text = element_text(colour = 'black'),
          legend.text = element_text(colour = 'black',size = 10),
          legend.title = element_text(colour = 'black',size = 12))+
    coord_flip()
  
  print(p3)
}
saveFig(file.path(plotDir,paste0(module_type,'_UP_enrichR')),plotFun_enrichR_res,
        rawData=allTerms_up,width = 5,height = 5.5,res = 500)



# Connectivity Map -------------------------------------------------------------


# BEAT AML ---------------------------------------------------------------------
res <- curate_beataml_setbp1_samples(
  donor_metadata_csv   = "Data/BEAT_AML/All_sample_data.csv",
  sample_metadata_xlsx = "Data/BEAT_AML/beataml_waves1to4_sample_mapping.xlsx",
  wes_mutations_txt    = "Data/BEAT_AML/beataml_wes_wv1to4_mutations_dbgap.txt",
  clinical_metadata_xlsx = "Data/BEAT_AML/beataml_wv1to4_clinical.xlsx",
  output_csv = "Data/BEAT_AML/BEATAML.v2_curated_selected_samples_mdat.csv"
)

curate_beataml_setbp1_samples <- function(
  donor_metadata_csv,
  sample_metadata_xlsx,
  wes_mutations_txt,
  clinical_metadata_xlsx,
  output_csv
) {
  ## ---------------------------
  ## Dependencies
  ## ---------------------------
  stopifnot(
    requireNamespace("data.table", quietly = TRUE),
    requireNamespace("readxl", quietly = TRUE),
    requireNamespace("stringr", quietly = TRUE)
  )
  
  library(data.table)
  
  ## ---------------------------
  ## Input checks
  ## ---------------------------
  files <- c(
    donor_metadata_csv,
    sample_metadata_xlsx,
    wes_mutations_txt,
    clinical_metadata_xlsx
  )
  
  missing_files <- files[!file.exists(files)]
  if (length(missing_files) > 0) {
    stop("Missing input files:\n", paste(missing_files, collapse = "\n"))
  }
  
  ## ---------------------------
  ## Donor-level metadata
  ## ---------------------------
  beat_donor <- fread(donor_metadata_csv, skip = 3)
  
  required_donor_cols <- c("PatientID", "corresponding_rnaseq_sample")
  if (!all(required_donor_cols %in% names(beat_donor))) {
    stop("Donor metadata missing required columns: ",
         paste(setdiff(required_donor_cols, names(beat_donor)), collapse = ", "))
  }
  
  beat_donor <- beat_donor[corresponding_rnaseq_sample != ""]
  
  ## ---------------------------
  ## Sample-level metadata
  ## ---------------------------
  beat_sample <- as.data.table(readxl::read_excel(sample_metadata_xlsx))
  
  required_sample_cols <- c("dbgap_subject_id", "dbgap_rnaseq_sample")
  if (!all(required_sample_cols %in% names(beat_sample))) {
    stop("Sample metadata missing required columns: ",
         paste(setdiff(required_sample_cols, names(beat_sample)), collapse = ", "))
  }
  
  beat_sample <- beat_sample[!is.na(dbgap_rnaseq_sample)]
  
  ## ---------------------------
  ## Merge donor + sample metadata
  ## ---------------------------
  setkey(beat_donor, PatientID, corresponding_rnaseq_sample)
  setkey(beat_sample, dbgap_subject_id, dbgap_rnaseq_sample)
  
  beat_mdat <- beat_sample[beat_donor, nomatch = 0]
  
  if (nrow(beat_mdat) == 0) {
    stop("No overlapping donor/sample records after merge.")
  }
  
  ## ---------------------------
  ## WES mutation data
  ## ---------------------------
  beat_mut <- fread(wes_mutations_txt)
  
  required_mut_cols <- c(
    "dbgap_sample_id", "symbol", "hgvsp_short",
    "hgvsc", "protein_position"
  )
  
  if (!all(required_mut_cols %in% names(beat_mut))) {
    stop("WES mutation file missing required columns: ",
         paste(setdiff(required_mut_cols, names(beat_mut)), collapse = ", "))
  }
  
  beat_mut[, `:=`(
    setbp1_mut = fifelse(symbol == "SETBP1", hgvsp_short, "none"),
    setbp1_mut_in_hotspot = fcase(
      symbol == "SETBP1" & protein_position == "1070/1596", "outside_hotspot",
      symbol == "SETBP1" & protein_position != "1070/1596", "inside_hotspot",
      default = "none"
    )
  )]
  
  beat_mut_sum <- beat_mut[
    , .(
      n_mut = uniqueN(hgvsc),
      setbp1_mut = paste(unique(setbp1_mut[setbp1_mut != "none"]), collapse = ","),
      setbp1_mut_in_hotspot = paste(
        unique(setbp1_mut_in_hotspot[setbp1_mut_in_hotspot != "none"]),
        collapse = ","
      )
    ),
    by = dbgap_sample_id
    ]
  
  beat_mut_sum[setbp1_mut == "", setbp1_mut := "none"]
  beat_mut_sum[setbp1_mut_in_hotspot == "", setbp1_mut_in_hotspot := "none"]
  
  ## Keep only samples with matched RNA/DNA
  beat_mut_sum <- beat_mut_sum[
    dbgap_sample_id %in% beat_mdat$corresponding_dnaseq_sample
    ]
  
  ## ---------------------------
  ## Clinical mutation data
  ## ---------------------------
  clinical_mut <- as.data.table(readxl::read_excel(clinical_metadata_xlsx))
  
  required_clin_cols <- c(
    "dbgap_subject_id",
    "dbgap_dnaseq_sample",
    "variantSummary"
  )
  
  if (!all(required_clin_cols %in% names(clinical_mut))) {
    stop("Clinical metadata missing required columns: ",
         paste(setdiff(required_clin_cols, names(clinical_mut)), collapse = ", "))
  }
  
  clinical_mut <- clinical_mut[!is.na(dbgap_dnaseq_sample)]
  
  clinical_mut[, setbp1_mut := fifelse(
    stringr::str_detect(variantSummary, "SETBP1"),
    sub("\\|.*$", "",
        stringr::str_extract(variantSummary, "SETBP1.*")),
    "none"
  )]
  
  clinical_mut[, setbp1_mut_in_hotspot := fcase(
    setbp1_mut != "none" & stringr::str_detect(setbp1_mut, "D868N"),
    "inside_hotspot",
    setbp1_mut != "none",
    "outside_hotspot",
    default = "none"
  )]
  
  clinical_sum <- clinical_mut[
    , .(
      setbp1_mut_long = paste(unique(setbp1_mut), collapse = ","),
      setbp1_mut_in_hotspot = paste(
        unique(setbp1_mut_in_hotspot), collapse = ","
      )
    ),
    by = dbgap_dnaseq_sample
    ]
  
  clinical_sum[, setbp1_mut := fifelse(
    setbp1_mut_long == "none",
    "none",
    gsub(" .*$|;.*$|,.*$", "",
         stringr::str_extract(setbp1_mut_long, "p\\..*"))
  )]
  
  clinical_sum[
    setbp1_mut_long == "SETBP1 (H1100R; MAF 51%)",
    setbp1_mut := "p.H1100R"
    ]
  
  ## ---------------------------
  ## Union clinical + WES mutation calls
  ## ---------------------------
  setkey(beat_mut_sum, dbgap_sample_id)
  setkey(clinical_sum, dbgap_dnaseq_sample)
  
  beat_mut_sum <- clinical_sum[
    beat_mut_sum,
    on = c(dbgap_dnaseq_sample = "dbgap_sample_id")
    ]
  
  beat_mut_sum[
    i.setbp1_mut != "none" & setbp1_mut == "none",
    `:=`(
      setbp1_mut = i.setbp1_mut,
      setbp1_mut_in_hotspot = i.setbp1_mut_in_hotspot
    )
    ]
  
  beat_mut_sum[, c(
    "i.setbp1_mut",
    "i.setbp1_mut_in_hotspot",
    "setbp1_mut_long"
  ) := NULL]
  
  ## ---------------------------
  ## Attach mutation info to main metadata
  ## ---------------------------
  setkey(beat_mdat, corresponding_dnaseq_sample)
  setkey(beat_mut_sum, dbgap_dnaseq_sample)
  
  beat_mdat <- beat_mut_sum[
    beat_mdat,
    on = c(dbgap_dnaseq_sample = "corresponding_dnaseq_sample")
    ]
  
  setnames(
    beat_mdat,
    c("n_mut", "setbp1_mut", "setbp1_mut_in_hotspot"),
    c("wes_n_mut", "wes_clin_setbp1_mut",
      "wes_clin_setbp1_mut_in_hotspot")
  )
  
  ## ---------------------------
  ## Merge full clinical metadata
  ## ---------------------------
  setnames(
    clinical_mut,
    old = names(clinical_mut),
    new = paste0("clin:", names(clinical_mut))
  )
  
  setkey(
    clinical_mut,
    "clin:dbgap_subject_id",
    "clin:dbgap_dnaseq_sample"
  )
  setkey(beat_mdat, PatientID, corresponding_dnaseq_sample)
  
  beat_mdat_clin <- clinical_mut[
    beat_mdat,
    on = c(
      `clin:dbgap_subject_id` = "PatientID",
      `clin:dbgap_dnaseq_sample` = "corresponding_dnaseq_sample"
    )
    ]
  
  ## ---------------------------
  ## Sample selection
  ## ---------------------------
  healthy <- beat_mdat_clin[rna_control != "No"]
  healthy[, sampleCat := "Healthy_control"]
  
  disease_setbp1 <- beat_mdat_clin[
    !is.na(wes_clin_setbp1_mut) & wes_clin_setbp1_mut != "none"
    ][, sampleCat := "Disease_SETBP1"]
  
  disease_young <- beat_mdat_clin[
    !is.na(`clin:ageAtDiagnosis`) &
      AgeCategory == "young" &
      `clin:ageAtDiagnosis` <= 10
    ][, sampleCat := "Disease_noneSETBP1"]
  
  beat_mdat_final <- rbindlist(
    list(disease_setbp1, disease_young, healthy),
    use.names = TRUE,
    fill = TRUE
  )
  
  ## ---------------------------
  ## Output
  ## ---------------------------
  fwrite(beat_mdat_final, output_csv)
  
  invisible(beat_mdat_final)
}


library(readxl)
### Generate metadata
# Donor-level metadata
beat_donor_mdat = read.csv('Data/BEAT_AML/All_sample_data.csv',skip = 3)
# Keep donor/samples with RNAseq data only 
beat_donor_mdat.sub = beat_donor_mdat[beat_donor_mdat$corresponding_rnaseq_sample != '',]
# --> results in 648 unique donorID, of which only 134 are in the beat_sample_mdat.sub object below
n_distinct(beat_donor_mdat.sub$PatientID)

# Sample-level metadata
beat_sample_mdat = read_excel('Data/BEAT_AML/beataml_waves1to4_sample_mapping.xlsx')
## check for number of overlapping patient ID
table(unique(beat_sample_mdat$dbgap_subject_id) %in% unique(beat_donor_mdat$PatientID))
unique(beat_sample_mdat$dbgap_subject_id)[!unique(beat_sample_mdat$dbgap_subject_id) %in% unique(beat_donor_mdat$PatientID)]
#3 Patient dbGapID are not present in All_sample_data.csv file. Ignore these cases as they arent SETBP1 mutation samples / or no DNA info anyway
# Keep donor/samples with RNAseq data only 
beat_sample_mdat.sub = beat_sample_mdat[!is.na(beat_sample_mdat$dbgap_rnaseq_sample),]
# --> results in 653 unique donorID (patientId)
n_distinct(beat_sample_mdat.sub$dbgap_subject_id)

# Merge donor and sample metadata
beat_mdat = merge(beat_donor_mdat,beat_sample_mdat,by.x=c('PatientID','corresponding_rnaseq_sample'),by.y=c('dbgap_subject_id','dbgap_rnaseq_sample'))
n_distinct(beat_mdat$PatientID) # 648 Patient ID

### BEAT_AML v1 includes samples from wave 1+2
### BEAT-AML v2 includes samples from wave 1+2, and additional samples from wave 3+4
### All relevant BEAT-AML v2 data (including their supplementary table S1) can be downloaded from github page: https://biodev.github.io/BeatAML2/
### Dataset includes Healthy control and Diseased samples. 
### Some patients have samples collected at multiple timepoints (eg. diagnosis, relapse, remission, residual etc.). However, not all timepoints have both RNA and DNA data
table(beat_donor_mdat$DiseaseStageAtSpecimenCollection)

### In terms of variants analysis:
### there are two files where variants/mutations were reported:
# 1. clinical summary: variants reported here are derived solely from clinical genotyping results (FAQ: https://biodev.github.io/BeatAML2/)
# 2. WES/targeted sequencing output: variants reported after processing and filtering steps as detailed in table S6-Variant Filtering from the Tyner et al. 2018 paper
### This means that the variants for the same DNA sampleID may not match. Clinically, SETBP1 may have not been included in the gene panel. WES data analysis may have filtered out the SETBP1 mutations detected clinically
### Additionally, not all DNA samples have matching RNA samples

### Strategy for selecting samples with SETBP1 mutation:
# 1. samples must have matching RNA and DNA data 
# 2. the DNA data was reported as having a SETBP1 mutations in either the "clinical summary" file, or WES_mutation file, as well as evidenced from paper 1 Tyner et al. 2018 (if patients were from wave 1+2)
# --> ie. taking the union of all the available information 
# This strategy results in 12 RNA samples, with supporting DNA data with SETBP1 mutation, from 11 different patients.
#   5 samples only have support from clinical data, not from WES_mutation. All of these 5 samples have SETBP1 mutations outside of the hotspot
#   1 sample have support from BOTH clinical and WES_mutations
#   6 samples only have support from WES_mutation data only





#### WES Mutations data
beat_mutation_mdat = read.delim('Data/BEAT_AML/beataml_wes_wv1to4_mutations_dbgap.txt')
# Check for SETBP1 mutation info
beat_mutation_mdat$setbp1_mut = ifelse(beat_mutation_mdat$symbol == 'SETBP1',beat_mutation_mdat$hgvsp_short,'none')
beat_mutation_mdat$setbp1_mut_in_hotspot = ifelse((beat_mutation_mdat$symbol == 'SETBP1'),ifelse(beat_mutation_mdat$protein_position == '1070/1596','outside_hotspot','inside_hotspot'),'none')
# Summarise mutation info by DNA sampleID
beat_mutation_mdat.summary = beat_mutation_mdat %>% group_by(dbgap_sample_id) %>% summarise(n_mut = n_distinct(hgvsc),setbp1_mut_in_hotspot=paste0(unique(setbp1_mut_in_hotspot),collapse = ','),
                                                                                            setbp1_mut=paste0(unique(setbp1_mut),collapse = ','))
beat_mutation_mdat.summary$setbp1_mut_in_hotspot = ifelse(grepl(',',beat_mutation_mdat.summary$setbp1_mut_in_hotspot),gsub('none,','',beat_mutation_mdat.summary$setbp1_mut_in_hotspot),as.character(beat_mutation_mdat.summary$setbp1_mut_in_hotspot))
beat_mutation_mdat.summary$setbp1_mut = ifelse(grepl(',',beat_mutation_mdat.summary$setbp1_mut),gsub('none,','',beat_mutation_mdat.summary$setbp1_mut),as.character(beat_mutation_mdat.summary$setbp1_mut))

# Keep only relevant DNA sampleID (these are samples with matching RNAsamples)
beat_mutation_mdat.sub = beat_mutation_mdat[beat_mutation_mdat$dbgap_sample_id %in% beat_mdat$corresponding_dnaseq_sample,]


#### Clinical Summary data
clinical_mut = read_excel('Data/BEAT_AML/beataml_wv1to4_clinical.xlsx')
n_distinct(clinical_mut$dbgap_subject_id) # 805 patients
# two patientIDs are not in donor_mdat, these are the same 2 patientID which are present in sample_mdat but NOT donor_mdat
table(clinical_mut$dbgap_subject_id[!clinical_mut$dbgap_subject_id %in% beat_donor_mdat$PatientID])

# Check for SETBP1 mutation info
clinical_mut$setbp1_mut = ifelse(grepl('SETBP1',clinical_mut$variantSummary),gsub("\\|.*$",'',str_extract(clinical_mut$variantSummary,pattern = 'SETBP1.*')),'none')
clinical_mut$setbp1_mut_in_hotspot = ifelse((clinical_mut$setbp1_mut != 'none'),ifelse(grepl('D868N',clinical_mut$setbp1_mut),'inside_hotspot','outside_hotspot'),'none')
clinical_mut.summary = clinical_mut %>% filter(!is.na(dbgap_dnaseq_sample)) %>% 
  group_by(dbgap_dnaseq_sample) %>% 
  summarise(setbp1_mut_in_hotspot=paste0(unique(setbp1_mut_in_hotspot),collapse = ','),
            setbp1_mut_long=paste0(unique(setbp1_mut),collapse = ','),
            setbp1_mut=ifelse(setbp1_mut_long == 'none','none',gsub(' .*$|;.*$|,.*$','',str_extract(setbp1_mut_long,'p\\..*'))))
clinical_mut.summary$setbp1_mut[clinical_mut.summary$setbp1_mut_long == 'SETBP1 (H1100R; MAF 51%)'] = 'p.H1100R'

# check that all 6 samples with SETBP1mut from clinical summary data are present in beat_mutation_mdat.summary
table(clinical_mut.summary$dbgap_dnaseq_sample[clinical_mut.summary$setbp1_mut != 'none'] %in% beat_mutation_mdat.summary$dbgap_sample_id)

# Update beat_mutation_mdat.summary
clinical_mut.summary$wes_setbp1_mut = beat_mutation_mdat.summary$setbp1_mut[match(clinical_mut.summary$dbgap_dnaseq_sample,beat_mutation_mdat.summary$dbgap_sample_id)]
w = clinical_mut.summary$dbgap_dnaseq_sample[(clinical_mut.summary$wes_setbp1_mut == 'none') & (clinical_mut.summary$setbp1_mut !='none')]

beat_mutation_mdat.summary$setbp1_mut[match(w,beat_mutation_mdat.summary$dbgap_sample_id)] = clinical_mut.summary$setbp1_mut[match(w,clinical_mut.summary$dbgap_dnaseq_sample)]
beat_mutation_mdat.summary$setbp1_mut_in_hotspot[match(w,beat_mutation_mdat.summary$dbgap_sample_id)] = clinical_mut.summary$setbp1_mut_in_hotspot[match(w,clinical_mut.summary$dbgap_dnaseq_sample)]


# Add SETBP1 mutation info to main metadata
beat_mdat$wes_n_mut = beat_mutation_mdat.summary$n_mut[match(beat_mdat$corresponding_dnaseq_sample,beat_mutation_mdat.summary$dbgap_sample_id)]
beat_mdat$wes_clin_setbp1_mut = beat_mutation_mdat.summary$setbp1_mut[match(beat_mdat$corresponding_dnaseq_sample,beat_mutation_mdat.summary$dbgap_sample_id)]
beat_mdat$wes_clin_setbp1_mut_in_hotspot = beat_mutation_mdat.summary$setbp1_mut_in_hotspot[match(beat_mdat$corresponding_dnaseq_sample,beat_mutation_mdat.summary$dbgap_sample_id)]

## Add additional information from clinical summary
colnames(clinical_mut) = paste0('clin:',colnames(clinical_mut))
beat_mdat_clin = merge(beat_mdat, clinical_mut,by.x=c('PatientID','corresponding_dnaseq_sample'),by.y=c('clin:dbgap_subject_id','clin:dbgap_dnaseq_sample'),all.x=T)



### Healthy control samples
# no age information for healthy control samples...
healthy_samples = beat_mdat_clin[beat_mdat_clin$rna_control != 'No',]



# Select samples to keep: 
# Keep samples with SETBP1 mutations 
# AND samples without SETBP1 mutation but in the 'young' category (10yrs or younger)
# I have checked that all young samples do not have SETBP1 mutations

# Overall, there are 62 samples from BEAT-AML v2 selected
set.seed(143)
df1 = beat_mdat_clin[!is.na(beat_mdat_clin$wes_clin_setbp1_mut) & beat_mdat_clin$wes_clin_setbp1_mut != 'none',]
df1$sampleCat = 'Disease_SETBP1'
df2 = beat_mdat_clin %>% filter(!is.na(`clin:ageAtDiagnosis`) & beat_mdat_clin$AgeCategory == 'young' & beat_mdat_clin$`clin:ageAtDiagnosis` <=18)
df2$sampleCat = 'Disease_noneSETBP1'

healthy_samples$sampleCat = 'Healthy_control'

beat_mdat.sub = rbind(df1,df2, 
                      healthy_samples)


write.csv(beat_mdat.sub,'Data/BEAT_AML/BEATAML.v2_cureated_selected_samples_mdat.csv')
write.csv(beat_mdat_clin,'Data/BEAT_AML/BEATAML.v2_cureated_ALL_samples_mdat.csv')

table(beat_mdat.sub$sampleCat,beat_mdat.sub$AgeCategory)

beat_sce_path = 'Data/BEAT_AML/BEAT_AML_sce.RDS'
if(!file.exists(beat_sce_path)){}
### Generate raw count
beat_raw_count = read.delim('Data/BEAT_AML/beataml_waves1to4_counts_dbgap.txt')
# Subset raw count matrix to include only samples with SETBP1 mutations
beat_raw_count = beat_raw_count[,colnames(beat_raw_count) %in% c(beat_mdat.sub$corresponding_rnaseq_sample,colnames(beat_raw_count)[1:4])]
rownames(beat_raw_count) = beat_raw_count$stable_id


### Ensembl build 75 gene models on GRCh37
geneAnn_BEAT_fp = 'Data/BEAT_AML/geneAnn_ens75_hg19_feb2014.csv'
if(file.exists(geneAnn_BEAT_fp)){
  geneAnn = read.table(geneAnn_BEAT_fp)
}else{
  library(biomaRt)
  library(EDASeq)
  mart=useMart("ensembl",host="https://feb2014.archive.ensembl.org", dataset = "hsapiens_gene_ensembl")
  geneAnn = getBM(c('ensembl_gene_id','external_gene_id','chromosome_name'),mart = mart)
  geneAnn_bp_length = getGeneLengthAndGCContent(geneAnn$ensembl_gene_id, "hsa")
  geneAnn = merge(geneAnn,geneAnn_bp_length,by.x = 'ensembl_gene_id',by.y=0)
  geneAnn = geneAnn[!is.na(geneAnn$length),]
  write.table(geneAnn,'Data/BEAT_AML/geneAnn_ens75_hg19_feb2014.csv',col.names = T)  
}

# Write individual bulk count file for deconvolution
beat_raw_count = merge(geneAnn,beat_raw_count,by.x = 'ensembl_gene_id',by.y=0)
rowData = beat_raw_count[,c(1:8)]

# save this table of raw count for BEAT-AML v2

sample_list = c()
rawCnt = tibble()
for(sample in colnames(beat_raw_count)[-c(1:8)]){
  sampleID = paste0('BEAT_AML_',sample,'_counts.txt')
  sample_list = c(sample_list,sampleID)
  tmp = beat_raw_count[,c('ensembl_gene_id','length',sample)]
  colnames(tmp)[1:2]=c('ensID',	'geneLength')
  
  if(!file.exists(paste0('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/','BEAT_AML_',sample,'_counts.txt'))){
    write.table(tmp,paste0('~/lustre_mt22/SETBP1/Results/June/CellSignalAnalysis/bulkData/','BEAT_AML_',sample,'_counts.txt'),sep = '\t',quote = F,row.names = F)  
  }
  
  if(ncol(rawCnt)==0){
    rawCnt = tmp
  }else{
    rawCnt = merge(rawCnt,tmp,by = c('ensID',	'geneLength'),all=T)
  }
  
}

# Prepare raw count and TPM count
rawCnt = column_to_rownames(rawCnt,var = 'ensID')
rawCnt = rawCnt[,-c(1)]

#if(!file.exists('~/lustre_mt22/SETBP1/Data/BEAT_AML/BEAT_AML_sub_tpm.csv')){
if(!file.exists('BEAT_AML_sub_tpm.csv')){
  geneLen = rowData$length
  
  rpk = apply(rawCnt, 2, function(x) x/(as.vector(geneLen)/1000))
  #normalize by the sample size using rpk values
  tpm = apply(rpk, 2, function(x) x / (sum(as.numeric(x)) / 1e6))  %>% as.data.frame()
  colNames = colnames(tpm)
  tpm$ensID = tmp$ensID
  tpm$geneLength = tmp$geneLength
  tpm = cbind(tpm[,c(ncol(tpm)-1,ncol(tpm))],tpm[,colNames])
  #write.csv(tpm,'~/lustre_mt22/SETBP1/Data/BEAT_AML/BEAT_AML_sub_tpm.csv')
  write.csv(tpm,'BEAT_AML_sub_tpm.csv')
}else{
  #tpm = read.delim('~/lustre_mt22/SETBP1/Data/BEAT_AML/BEAT_AML_sub_tpm.csv',sep = ',')
  tpm = read.delim('BEAT_AML_sub_tpm.csv',sep = ',')
  if('X' %in% colnames(tpm)){
    tpm=tpm[,-c(colnames(tpm) == 'X')]
  }
}

# Extract SETBP1 expression level in TPM
setbp1_beatAML = tpm[tpm$ensID == 'ENSG00000152217',-c(1,2)] %>% t()

### Generate metadata
m = match(colnames(rawCnt),beat_mdat.sub$corresponding_rnaseq_sample)
beat_mdat.sub = beat_mdat.sub[m,]

beat_mdat = data.frame(source = rep('BEAT_AML',nrow(beat_mdat.sub)),
                       sample_file = sample_list,
                       sampleID = beat_mdat.sub$corresponding_rnaseq_sample,
                       caseID = beat_mdat.sub$PatientID,
                       age = beat_mdat.sub$AgeCategory,
                       clin_age = beat_mdat.sub$`clin:ageAtDiagnosis`,
                       sex = beat_mdat.sub$Gender,
                       cancerType = beat_mdat.sub$SpecificDxAtAcquisition,
                       cancerSubClass = beat_mdat.sub$FAB_BlastMorphology,
                       setbp1_mut = beat_mdat.sub$wes_clin_setbp1_mut,
                       setbp1_mut_inHotspot = beat_mdat.sub$wes_clin_setbp1_mut_in_hotspot,
                       setbp1_tpm = setbp1_beatAML[,1])

# Generate a singleCellExperiment instance
beat_sce <- SingleCellExperiment(list(counts_raw=rawCnt, counts_tpm=tpm[,-c(1,2)]),
                                 colData=beat_mdat,
                                 rowData=rowData,
                                 metadata=list(genome_version = 'ensemble_release_75_feb2014',
                                               beat_donor_mdat = beat_donor_mdat,
                                               beat_sample_mdat = beat_sample_mdat,
                                               beat_mutation_mdat = beat_mutation_mdat))
saveRDS(beat_sce,beat_sce_path)    
}else{
  beat_sce = readRDS(beat_sce_path)
  beat_mdat = as.data.frame(colData(beat_sce))
}
beat_mdat2 = beat_mdat[,colnames(beat_mdat) !='age']
colnames(beat_mdat2)[colnames(beat_mdat2) == 'clin_age'] = 'age'
bulk_samples = rbind(bulk_samples,beat_mdat2)
}
## Score the gene signature in BEATAML samples - with and without SETBP1 -------
beat_mdat_clin$sampleCat = ifelse(beat_mdat_clin$SampleID %in% beat_mdat.sub$SampleID,
                                  beat_mdat.sub$sampleCat[match(beat_mdat_clin$SampleID, beat_mdat.sub$SampleID)],
                                  paste0("Disease_noneSETBP1::",beat_mdat_clin$AgeCategory)
                                  #'-'
                                  )
table(beat_mdat_clin$sampleCat)
data = beat_mdat_clin
beat_mdat = data.frame(source = rep('BEAT_AML',nrow(data)),
                       sampleID = data$corresponding_rnaseq_sample,
                       caseID = data$PatientID,
                       age = data$AgeCategory,
                       clin_age = data$`clin:ageAtDiagnosis`,
                       sex = data$Gender,
                       cancerType = data$SpecificDxAtAcquisition,
                       cancerSubClass = data$FAB_BlastMorphology,
                       setbp1_mut = data$wes_clin_setbp1_mut,
                       setbp1_mut_inHotspot = data$wes_clin_setbp1_mut_in_hotspot,
                       sample_cat = data$sampleCat)
# harmonised counts
beat_harm_count = fread('Data/BEAT_AML/beataml_waves1to4_norm_exp_dbgap.txt')
beat_harm_count = beat_harm_count[,colnames(beat_harm_count) %in% c(beat_mdat.sub$corresponding_rnaseq_sample,colnames(beat_harm_count)[1:4]),with=FALSE]
rownames(beat_harm_count) = beat_harm_count$stable_id

mtx <- beat_harm_count[,
                       colnames(beat_harm_count) %in% beat_mdat_clin$corresponding_rnaseq_sample,
                       with=FALSE] %>% as.matrix()
rownames(mtx) = beat_harm_count$display_label

## ----- 4. Rank genes -----
ranked <- singscore::rankGenes(mtx)

## ----- 5. Check module validity -----
check_module_validity(module_list = gene_module_list, module_name = names(gene_module_list), min_genes=50)

## ----- 6. Get gene sets for scoring -----
moduleList <- list('g_lymphoid' = list('up' = gene_module_list[['g_lymph_up']],
                                      'down' = gene_module_list[['g_lymph_down']]),
                  'g_myeloid' = list('up' = gene_module_list[['g_myeloid_up']],
                                     'down' = gene_module_list[['g_myeloid_down']]),
                  's_myeloid' = list('up' = gene_module_list[['s_up']],
                                     'down' = gene_module_list[['s_down']]),
                  's_specific_up' = list('up' = gene_module_list[['s_up']]),
                  's_myeloid_shared_up' = list('up' = gene_module_list[['s_myeloid_shared_up']]),
                  'g_myeloid_specific_up' = list('up' = gene_module_list[['g_myeloid_specific_up']])
                  )

## ----- 7. Calculate module scores -----
allScore <- data.frame()

for (i in seq_along(moduleList)) {
  module_name <- names(moduleList)[i]
  
  if(all(c('up','down') %in% names(moduleList[[i]]))){
    dims <- singscore::simpleScore(ranked, upSet = moduleList[[i]]$up, downSet = moduleList[[i]]$down)
    dims <- dims[, c('TotalScore', 'TotalDispersion'), drop = FALSE]
  } else {
    # For other modules, simpleScore with only upSet (no downSet)
    # Validate module genes count here optionally
    if (length(moduleList[[i]]$up) == 0) {
      warning(sprintf("Module '%s' is empty. Skipping.", module_name))
      next
    }
    if (length(moduleList[[i]]$up) < min_genes) {
      warning(sprintf("Module '%s' has fewer than %d genes. Skipping.", module_name, min_genes))
      next
    }
    dims <- singscore::simpleScore(ranked, upSet = moduleList[[i]]$up)
  }
  
  # Merge scores with metadata
  scoredf <- merge(beat_mdat, dims, by.x = 'sampleID', by.y = 0)
  scoredf$module_name <-module_name
  
  # Append to allScore
  allScore <- rbind(allScore, scoredf)
}

# Box plot of the score
allScore$setbp1_mut_inHotspot[is.na(allScore$setbp1_mut_inHotspot)] = 'WT'
allScore$cancerSubClass
allScore$dbgap_sample_id = data$corresponding_dnaseq_sample[match(allScore$sampleID,data$corresponding_rnaseq_sample)]
allScore$n_mut = beat_mutation_mdat.summary$n_mut[match(allScore$dbgap_sample_id,beat_mutation_mdat.summary$dbgap_sample_id)]
ggplot(allScore[allScore$module_name == 's_myeloid_shared_up',],aes(TotalScore,n_mut))+
  geom_point()+
  facet_wrap(vars(sample_cat))

df = allScore[allScore$cancerType %in% beat_mdat$cancerType[beat_mdat$sample_cat == 'Disease_SETBP1'] &
                allScore$module_name == 's_myeloid_shared_up'  ,]

module_type = 's_myeloid_shared_up'
df = allScore[allScore$module_name == 's_myeloid_shared_up',]
df$is_healthy = factor(ifelse(grepl('ealthy',df$cancerType),'Healthy_control','Disease'),
                       c('Healthy_control','Disease'))
df$cancerType[df$cancerType == 'unknown'] = 'Unknown'

p <- ggplot(df,
       aes(reorder(cancerType,TotalScore,median),TotalScore))+
  geom_boxplot(outlier.colour = 'white',position = 'dodge', alpha = 1,
               width=0.5,linewidth=0.25,col='black') +
  geom_quasirandom(size=0.3,width = 0.18,alpha=0.7,col=grey(0.7,0.5))+
  geom_quasirandom(data = df[df$setbp1_mut_inHotspot %in% c('inside_hotspot','outside_hotspot'),],
                   aes(col=setbp1_mut_inHotspot),size=1.4,width = 0.18,alpha=1)+
  scale_color_manual(values = c('WT'=grey(0.7,0.5),
                                'none'=grey(0.7,0.5),
                                'outside_hotspot'='#2D4372',
                                'inside_hotspot'=col25[5])) +
  facet_grid(module_name~is_healthy,scales = 'free',space='free')+
  theme_classic(base_size = 8)+
  ggtitle('')+xlab('')+ylab('Module score')+
  theme(#legend.position = 'none',
    panel.border = element_rect(fill=F,colour = 'black',linewidth = 0.2),
    axis.line = element_blank(),
    strip.background=element_rect(linewidth=0),
    strip.text = element_text(size = 10,colour = 'black'),
    axis.ticks = element_line(colour = 'black',linewidth = 0.2),
    axis.text.x = element_text(angle = 90, vjust = 0.5,hjust = 1,colour = 'black'),
    axis.text = element_text(colour = 'black'),
    axis.title = element_text(size=12))

plotFun_beat_aml = function(noFrame=FALSE,noPlot=FALSE){
  print(p)
}
saveFig(file.path(plotDir,paste0(module_type,'_singScore_BeatAML2')),
        plotFun_beat_aml,rawData=df,width = 6.4,height = 5,res = 500)


# BEAT-AML Drug sensitivity ----------------------------------------------------
drug_syn = readxl::read_excel('Data/BEAT_AML/beataml_drug_families.xlsx')
drug_family = readxl::read_excel('Data/BEAT_AML/beataml_drug_families.xlsx',sheet=2)
auc_dt = fread('Data/BEAT_AML/beataml_probit_curve_fits_v4_dbgap.txt')

auc_dt = cbind(auc_dt,beat_mdat[match(auc$dbgap_rnaseq_sample,beat_mdat$sampleID),])
auc_dt$setbp1_status = ifelse(is.na(auc_dt$sample_cat),'unknown',
                              ifelse(auc_dt$sample_cat == 'Disease_SETBP1','yes','no'))
dt = auc_dt %>% 
  group_by(setbp1_status,inhibitor) %>% 
  summarise(med = median(auc)) 
dt2 = dt %>% 
  #filter(setbp1_status != 'unknown') %>% 
  pivot_wider(id_cols = 'inhibitor',names_from = 'setbp1_status',values_from = 'med') %>% 
  mutate(delta = no-yes)
dt2 = dt2[order(dt2$delta,decreasing = T),]

auc.top10 = auc_dt[auc_dt$inhibitor %in% dt2$inhibitor[1:10],]

ggplot(auc_dt[auc_dt$inhibitor == 'Panobinostat',],
       aes(setbp1_status,auc,fill=sample_cat))+
  geom_boxplot(outlier.shape = NA)+
  geom_quasirandom(size=0.5)+
  facet_wrap(vars(inhibitor))+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))

auc_dt$s_myeloid_shared_score = score$TotalScore[match(auc_dt$dbgap_rnaseq_sample,score$sampleID)]
ggplot(auc_dt[auc_dt$inhibitor == 'Panobinostat' & auc_dt$sample_cat != '-',],
       aes(s_myeloid_shared_score,auc,col=sample_cat))+
  geom_point(aes(shape = setbp1_mut_inHotspot),size=3)+
  facet_wrap(vars(inhibitor))+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))

# Find which drug correlates with the gene signature score
score = allScore[allScore$module_name == 's_myeloid_shared_up',]
cor_score = data.frame()
table(score$sampleID %in% auc_dt$dbgap_rnaseq_sample)
table(auc_dt$dbgap_rnaseq_sample %in% score$sampleID)

samples_to_keep = intersect(score$sampleID,auc_dt$dbgap_rnaseq_sample)
samples_to_keep = samples_to_keep[samples_to_keep %in% score$sampleID[score$sample_cat == 'Disease_SETBP1']]
for(s in unique(auc_dt$inhibitor)){
  d = auc_dt[auc_dt$inhibitor == s & auc_dt$dbgap_rnaseq_sample %in% samples_to_keep,]
  sig_score = score[match(d$dbgap_rnaseq_sample,score$sampleID),]
  if(nrow(cor_score) == 0){
    cor_score = data.frame(drug = s, 
                           cor = cor(d$auc,sig_score$TotalScore),
                           n_sample = n_distinct(d$dbgap_rnaseq_sample))
  }
  cor_score = rbind(cor_score,data.frame(drug = s, 
                                         cor = cor(d$auc,sig_score$TotalScore),
                                         n_sample = n_distinct(d$dbgap_rnaseq_sample)))
}

cor_score_sub = cor_score[cor_score$n_sample >= 20,]
cor_score_sub = cbind(cor_score_sub,drug_family[match(cor_score_sub$drug,drug_family$inhibitor),])

cor_score_sub$drug[!is.na(cor_score_sub$cor) & cor_score_sub$cor <= -0.4]
cor_score_sub$drug[!is.na(cor_score_sub$cor) & cor_score_sub$cor >= 0.4]
