# Import top signatures --------------------------------------------------------
table(germline_sig$DE.Direction,germline_sig$Lineage)
table(germline_sig$Lineage,germline_sig$is_in_Shared_SETBP1_module)
table(somatic_sig$DE.Direction,somatic_sig$is_in_Shared_SETBP1_module)

cl_2 = read.delim('Results/06_cmap/cellinfo_beta.txt')
compound = read.delim('Results/06_cmap/compoundinfo_beta.txt')
compound_dict = read.delim(file.path(outDir,'compoundinfo_beta.txt'))
table(cp$pert_id %in% compound_dict$pert_id)
table(compound_dict$pert_id %in% cp$pert_id)
dplyr::n_distinct(cp$pert_id)
n_distinct(cp$moa)
View(table(cp$moa))
moa = unique(cp$moa)
table(grepl('NFKB|nfkb|NF-KB',moa))
cp = cp[grepl('NFKB|nfkb|NF-KB',cp$moa),]

# Plot distribution of ranking by MOA
hist(cp_toKeep$norm_cs_ranking[grepl('CDK inhibitor',cp_toKeep$moa)])

dt = data.table::setDT(cp.sub)
cp_agg <- dt[
  ,.(
    pert_iname = first(pert_iname),
    moa        = first(moa),
    norm_cs_mean = mean(norm_cs, na.rm = TRUE),
    norm_cs_median = median(norm_cs, na.rm = TRUE),
    n_cell_lines =.N
  ),
  by = pert_id
  ]

# Choose which aggregated score you want to rank by:
cp_agg[, `:=`(
  score = -norm_cs_mean,
  is_nfkb = grepl("nfkb|NFkB|NF-?κB|NFKB", moa, ignore.case = TRUE))] 

library(data.table)
# Split MOA on '|' and unnest
cp_long <- cp_agg %>%
  mutate(moa = ifelse(is.na(moa) | moa == "", NA_character_, moa)) %>%
  tidyr::separate_rows(moa, sep = "\\|") %>%
  filter(!is.na(moa) & moa != "") %>%
  as.data.table()

med_score = cp_long %>% group_by(moa) %>% summarise(med = median(score))
med_score = med_score[order(med_score$med),]
med_score$idx = seq(1:nrow(med_score))

ggplot(cp_agg[cp_agg$moa %in% c(head(med_score$moa,10),
                                med_score$moa[420:440]
                                tail(med_score$moa,40),
                                med_score$moa[grepl('NFKB',med_score$moa)]),],aes(score,reorder(moa,score, FUN = median)))+
  geom_boxplot(aes(fill = is_nfkb))+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust = 0.5,hjust = 1))


# For GSEA-like analysis we need a ranked vector: names = compound IDs, values = scores
# Filter out any NA scores
cp_long_clean <- cp_long[!is.na(score)]

# Choose one row per compound for the ranked list (we just need one score per compound)
# If compounds have multiple MOAs, we still only want one score per compound in the ranking
cp_ranked <- cp_long_clean[,.(score = first(score)),
                           by = pert_id
                           ]

# Order by descending score (high score = more strongly opposing/mimicking your module, depending on definition)
cp_ranked <- cp_ranked[order(-score)]

# Create the named vector
ranks <- cp_ranked$score
names(ranks) <- cp_ranked$pert_id

# Build a list: MOA -> vector of pert_id
moa_sets <- cp_long_clean[,.(pert_ids = list(unique(pert_id))), by = moa]
moa_sets_list <- setNames(moa_sets$pert_ids, moa_sets$moa)

# Example: look at size of each MOA set
moa_sizes <- sapply(moa_sets_list, length)
moa_sizes[order(-moa_sizes)][1:20]

min_size <- 3
moa_sets_list_filt <- moa_sets_list[moa_sizes >= min_size]

library(fgsea)

set.seed(123)

fgsea_res <- fgsea(
  pathways = moa_sets_list_filt,
  stats    = ranks,
  minSize  = min_size,
  maxSize  = 500,  # tune if needed
  nperm    = 10000 # more permutations = smoother p-values, but slower
)

# Order by NES (normalized enrichment score) or p-value
fgsea_res <- fgsea_res[order(padj, NES)]
head(fgsea_res)

# Look for NFkB-related MOAs
nfkb_idx <- grepl("nfkb|NFkB|NF-?κB|NFKB", fgsea_res$pathway, ignore.case = TRUE)
fgsea_res[nfkb_idx][order(padj)]
# Suppose your NFkB-related MOA set is named exactly "NFkB inhibitor" (adjust pattern as needed)
nfkb_moa_name <- fgsea_res$pathway[which.min(fgsea_res$padj[nfkb_idx])]  # or choose by hand
nfkb_moa_name <- "NFkB inhibitor"

# Get the set of pert_ids for that MOA
nfkb_set <- moa_sets_list_filt[[nfkb_moa_name]]

# Plot enrichment
plotEnrichment(nfkb_set, ranks) +
  ggplot2::ggtitle(paste0("GSEA-like enrichment for MOA: ", nfkb_moa_name))

library(ggplot2)

# Add rank and quantile groups to cp_ranked
cp_ranked[, rank :=.I]
n <- nrow(cp_ranked)
cp_ranked[, bin := cut(rank,
                       breaks = c(0, 0.05*n, 0.2*n, n),
                       labels = c("top_5pct", "5-20pct", "rest"))]

# Merge MOAs back in
cp_ranked_long <- merge(cp_ranked, cp_long_clean[,.(pert_id, moa)], by = "pert_id", all.x = TRUE)
cp_ranked_long <- cp_ranked_long[!is.na(moa)]

# Example: Over-representation of MOA in top 5%
tab_top <- cp_ranked_long[bin == "top_5pct",.N, by = moa]
tab_bg  <- cp_ranked_long[,.N, by = moa]

# Build a simple enrichment statistic: proportion in top vs background
enrich_dt <- merge(tab_top, tab_bg, by = "moa", suffixes = c("_top", "_all"))
enrich_dt[, prop_top := N_top / sum(N_top)]
enrich_dt[, prop_all := N_all / sum(N_all)]
enrich_dt[, enrichment := prop_top / prop_all]

# Optional: filter to MOAs with enough compounds
enrich_dt <- enrich_dt[N_all >= 3]

# Highlight NFkB-related MOAs
enrich_dt[, is_nfkb := grepl("nfkb|NFkB|NF-?κB|NFKB", moa, ignore.case = TRUE)]

ggplot(enrich_dt, aes(x = enrichment, y = reorder(moa, enrichment), color = is_nfkb)) +
  geom_point() +
  scale_color_manual(values = c(`TRUE` = "red", `FALSE` = "grey50")) +
  theme_bw() +
  labs(x = "Enrichment in top 5% (prop_top / prop_all)",
       y = "MOA",
       color = "NFkB-related",
       title = "MOA enrichment among top-ranked compounds")
# Look at dist ----------------------------------------------------------------