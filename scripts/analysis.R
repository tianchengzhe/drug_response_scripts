# rm(list=ls()), restart R, then execute rmarkdown::render('analysis.R') to generate the html report
# set up working space, load library
setwd('/Users/chti4479/scRNAseq_chen/');
library(dplyr); library(Seurat); library(cellrangerRkit);
library(GeneOverlap);
set.seed(1);

## Step 1. Load Raw Data and QC
gene_bc_matrix_T <- load_cellranger_matrix('/spencerlab/spencergrp/chya7855/RNA-seq/scRNA-seq_1st_library_resequence/T');
gene_bc_matrix_UT <- load_cellranger_matrix('/spencerlab/spencergrp/chya7855/RNA-seq/scRNA-seq_1st_library_resequence/UT'); 
gene_bc_matrix_aggr <- load_cellranger_matrix('/spencerlab/spencergrp/chya7855/RNA-seq/scRNA-seq_1st_library_resequence/aggr'); 

# genes should be detected in both datasets, and no two probes maps to the same gene
use_genes <- setdiff(intersect(get_nonzero_genes(gene_bc_matrix_T), get_nonzero_genes(gene_bc_matrix_UT)), which(duplicated(fData(gene_bc_matrix_T)[,'symbol'])))
expr_matrix_aggr <- exprs(gene_bc_matrix_aggr);
rm(gene_bc_matrix_T, gene_bc_matrix_UT);

# criteria 1. high library size
if_outlier_counts <- scater::isOutlier(Matrix::colSums(expr_matrix_aggr), nmads=5, type="lower", log=TRUE);
# criteria 2. high feature size
if_outlier_features <- scater::isOutlier(Matrix::colSums(expr_matrix_aggr>0), nmad=5, type="lower", log=TRUE);
# criteria 3. low MT expression
if_outlier_MT <- scater::isOutlier(Matrix::colSums(expr_matrix_aggr[grepl("^MT-", fData(gene_bc_matrix_aggr)[,2]),])/Matrix::colSums(expr_matrix_aggr), nmads=5, type="higher");
use_cells <- !(if_outlier_counts | if_outlier_features | if_outlier_MT); 
gene_bc_matrix_aggr <- gene_bc_matrix_aggr[use_genes, use_cells];

# extract data
aggr_data <- exprs(gene_bc_matrix_aggr);
rownames(aggr_data) <- fData(gene_bc_matrix_aggr)[, 'symbol'];
T_cell_id <- which(grepl("-1", colnames(aggr_data)));
UT_cell_id <- which(grepl("-2", colnames(aggr_data)));

# load pre-computed tSNE coordinates, run by cellranger on filtered data (verified to be deterministic)
analysis_results_aggr_filtered <- load_cellranger_analysis_results('/spencerlab/spencergrp/chya7855/RNA-seq/scRNA-seq_1st_library_resequence/aggr_filtered');
cellranger_tSNE <- analysis_results_aggr_filtered$tsne;

## Step 2. Compute Proliferation Probabilities and Classification
gene.list <- c('CCNA2','HJURP','UBE2C','KIF23','TOP2A','CDK1','PLK1','TPX2','NUF2','HMMR',
               'CDC20','NUSAP1','KIF11','BUB1B','TACC3','BIRC5','KIF20B','CDKN3','ARHGEF39','CCNB1',
               'ESCO2','RRM2','GTSE1','MKI67','MCM5','FAM83D','ATAD2','CENPA','AURKB','PRC1',
               'MAD2L1','GMNN','CDC45','BRCA1','CCNF','CENPE','BUB1','CCNB2','CCNE1','CCNE2',
               'CDC6','MCM2','CDC25A','CDK2','E2F1','MYBL2','MCM4','MCM6','RFC3','RFC4','KIF4A');

# function to compute probabilities
compute_prob_P <- function(expr_matrix, gene.list, IF_percent_P, epsilon) {
  # compute the probability of detecting each gene
  expr_matrix <- as.matrix(expr_matrix[gene.list, ]);
  prob_proliferative_genes <- Matrix::rowMeans(expr_matrix > 0); ## WAS >
  prob_detecting_genes <- (prob_proliferative_genes-epsilon*(1-IF_percent_P))/IF_percent_P;
  prob_detecting_genes[prob_detecting_genes > 1-epsilon] <- 1-epsilon;
  prob_detecting_genes[prob_detecting_genes < epsilon] <- epsilon;
  
  # compute gene profile for P and Q cases
  prob_P_given_geneprofile <- apply(expr_matrix, 2, function(curr_expr) {
    prob_geneprofile_P <- (curr_expr>0)*prob_detecting_genes + (curr_expr==0)*(1-prob_detecting_genes);
    prob_geneprofile_Q <- (curr_expr>0)*epsilon + (curr_expr==0)*(1-epsilon);
    return(prod(prob_geneprofile_P)*IF_percent_P / (prod(prob_geneprofile_P)*IF_percent_P + prod(prob_geneprofile_Q)*(1-IF_percent_P)));
  });
}

# compute probabilities of being proliferative
IF_percent_P_UT <- 0.95; # from IF data, untreated
IF_percent_P_T <- 0.18; # from IF data, treated
epsilon <- 0.01; # probability of detecing a gene in a quiescent cell (proliferative gene markers)
prob_P_given_geneprofile_aggr_T <- compute_prob_P(aggr_data[,T_cell_id], gene.list, IF_percent_P_T, epsilon);
prob_P_given_geneprofile_aggr_UT <- compute_prob_P(aggr_data[,UT_cell_id], gene.list, IF_percent_P_UT, epsilon);

# plot subpopulations by category on tSNE plot
category_boundary <- c(-Inf, -60, -40, -20, 0, Inf);
category_per_cell <- sapply(c(prob_P_given_geneprofile_aggr_T, prob_P_given_geneprofile_aggr_UT), function(x) {
  return(which(x >= exp(category_boundary[1:5]) & x < exp(category_boundary[2:6])))
});
print(data.frame(category=1:5, 
                 counts_treated=hist(category_per_cell[T_cell_id], breaks=0:5, plot=FALSE)$counts,
                 counts_untreated=hist(category_per_cell[UT_cell_id], breaks=0:5, plot=FALSE)$counts));
df_to_plot <- data.frame(TSNE.1=cellranger_tSNE$TSNE.1,
                         TSNE.2=cellranger_tSNE$TSNE.2,
                         category=category_per_cell,
                         cat12=ifelse(category_per_cell<3, as.character(category_per_cell), 'Others'),
                         cat5=ifelse(category_per_cell==5, '5', 'Others'));

p <- ggplot(df_to_plot, aes(x=TSNE.1, y=TSNE.2, col=cat12)) + geom_point(size = 0.5);
ggsave(plot=p, filename='category1to2.pdf', height=3, width=4, useDingbats=FALSE);

p <- ggplot(df_to_plot, aes(x=TSNE.1, y=TSNE.2, col=cat5)) + geom_point(size = 0.5);
ggsave(plot=p, filename='category5.pdf', height=3, width=4, useDingbats=FALSE);

# partition cells into finer bins, Fig 4A
bin_border <- c(-80:0);
bin_category <- sapply(bin_border, function(x) which(x>=category_boundary[1:5] & x<category_boundary[2:6]));
counts_T <- sapply(bin_border, function(x) sum(log(prob_P_given_geneprofile_aggr_T)>=x & log(prob_P_given_geneprofile_aggr_T) < x+1));
counts_UT <- sapply(bin_border, function(x) sum(log(prob_P_given_geneprofile_aggr_UT)>=x & log(prob_P_given_geneprofile_aggr_UT) < x+1));
print(data.frame(left_border=bin_border, counts_treated=counts_T, counts_untreated=counts_UT));

# plot histogram of probabilities (barplot to control bin counting)
pdf('hist_prob_T.pdf', width=4, height=2);
ggplot(data.frame(left=bin_border, counts=counts_T, cat=as.character(bin_category)), 
       aes(x=left, y=counts, fill=cat)) +
  geom_bar(stat="identity") +
  xlab('ln(Probability of Being Proliferative)') +
  ylab('Cell Counts') + coord_cartesian(ylim=c(0,175));
dev.off();
pdf('hist_prob_UT.pdf', width=4, height=2);
ggplot(data.frame(left=bin_border, counts=counts_UT, cat=as.character(bin_category)), 
       aes(x=left, y=counts, fill=cat)) +
  geom_bar(stat="identity") +
  xlab('ln(Probability of Being Proliferative)') +
  ylab('Cell Counts') + coord_cartesian(ylim=c(0,600));
dev.off();

# define subpopulations and compute mean expression
HC_PC_T <- T_cell_id[prob_P_given_geneprofile_aggr_T == 1];
HC_QC_T <- T_cell_id[prob_P_given_geneprofile_aggr_T < exp(-40)];
HC_PC_UT <- UT_cell_id[prob_P_given_geneprofile_aggr_UT == 1];
temp <- Matrix::colSums(aggr_data); normalized_raw_expr <- t(t(aggr_data)*median(temp)/temp); 
get_expr_stat <- function(expr_mat, cell_id) {
  expr_mat <- expr_mat[, cell_id];
  mean_val <- Matrix::rowMeans(expr_mat);
  sem_val <- sqrt(Matrix::rowMeans(expr_mat**2) - mean_val**2)/sqrt(length(cell_id));
  return(list(mean_val, sem_val));
}
temp <- get_expr_stat(normalized_raw_expr, HC_PC_T); mean_expr_HC_PC_T <- temp[[1]]; sem_expr_HC_PC_T <- temp[[2]];
temp <- get_expr_stat(normalized_raw_expr, HC_QC_T); mean_expr_HC_QC_T <- temp[[1]]; sem_expr_HC_QC_T <- temp[[2]];
temp <- get_expr_stat(normalized_raw_expr, HC_PC_UT); mean_expr_HC_PC_UT <- temp[[1]]; sem_expr_HC_PC_UT <- temp[[2]];
temp <- get_expr_stat(normalized_raw_expr, T_cell_id); mean_expr_allT <- temp[[1]]; sem_expr_allT <- temp[[2]];
temp <- get_expr_stat(normalized_raw_expr, UT_cell_id); mean_expr_allUT <- temp[[1]]; sem_expr_allUT <- temp[[2]];

## Step 3. Differential expression analysis
# create Seurat object
aggr_obj <- CreateSeuratObject(counts=aggr_data, project='aggr_obj', min.cells=3, min.features=200);
aggr_obj <- NormalizeData(aggr_obj, normalization.method ='LogNormalize', scale.factor=10000);
aggr_obj <- FindVariableFeatures(aggr_obj, selection.method = "vst", nfeatures = 2000);
aggr_obj <- ScaleData(aggr_obj, features = rownames(aggr_obj));

# PCA
aggr_obj <- RunPCA(aggr_obj, features = VariableFeatures(object = aggr_obj));
num_selected_PC <- 15; # select the first 15 PCs

# clustering
aggr_obj <- FindNeighbors(aggr_obj, dims = 1:num_selected_PC);
aggr_obj <- FindClusters(aggr_obj);
aggr_obj <- RunTSNE(aggr_obj, dims = 1:num_selected_PC);

# differential expression analysis
aggr_obj <- SetIdent(object = aggr_obj, cells=HC_PC_T, value='HC_PC_T');
aggr_obj <- SetIdent(object = aggr_obj, cells=HC_QC_T, value='HC_QC_T');
aggr_obj <- SetIdent(object = aggr_obj, cells=HC_PC_UT, value='HC_PC_UT');
aggr_obj <- SetIdent(object = aggr_obj, cells=setdiff(1:ncol(aggr_obj), c(HC_PC_T, HC_QC_T, HC_PC_UT)), value='Others');
escapee_vs_nonescapee <- FindMarkers(aggr_obj, ident.1 = "HC_PC_T", ident.2 = "HC_QC_T", min.pct = 0.01);
escapee_vs_ut <- FindMarkers(aggr_obj, ident.1 = "HC_PC_T", ident.2 = "HC_PC_UT", min.pct = 0.01);

# count DEG
EvsNE_upreg <- rownames(escapee_vs_nonescapee)[escapee_vs_nonescapee$avg_logFC > 0 & escapee_vs_nonescapee$p_val_adj < 0.05];
EvsNE_downreg <- rownames(escapee_vs_nonescapee)[escapee_vs_nonescapee$avg_logFC < 0 & escapee_vs_nonescapee$p_val_adj < 0.05];
EvsUT_upreg <- rownames(escapee_vs_ut)[escapee_vs_ut$avg_logFC > 0 & escapee_vs_ut$p_val_adj < 0.05];
EvsUT_downreg <- rownames(escapee_vs_ut)[escapee_vs_ut$avg_logFC < 0 & escapee_vs_ut$p_val_adj < 0.05];
intersect_upreg <- intersect(EvsNE_upreg, EvsUT_upreg);
intersect_downreg <- intersect(EvsNE_downreg, EvsUT_downreg);

## Step 4. Outputs of DEG
# generate output tables for differentially expressed genes
generate_deg_table <- function(query_genes, filename, 
                               mean_expr_HC_PC_T, mean_expr_HC_QC_T, mean_expr_HC_PC_UT, 
                               escapee_vs_nonescapee, escapee_vs_ut) {
  query_genes <- query_genes[query_genes %in% names(mean_expr_HC_PC_T)];
  data_to_output <- data.frame(
    gene.name=query_genes,
    mean.expr.HC_PC_T=mean_expr_HC_PC_T[query_genes],
    mean.expr.HC_QC_T=mean_expr_HC_QC_T[query_genes],
    mean.expr.HC_PC_UT=mean_expr_HC_PC_UT[query_genes],
    logfold.change.HC_PC_T.vs.HC_QC_T=escapee_vs_nonescapee[query_genes, 'avg_logFC'],
    logfold_change.HC_PC_T.vs.HC_PC_UT=escapee_vs_ut[query_genes, 'avg_logFC'],
    pval.HC_PC_T.vs.HC_QC_T=escapee_vs_nonescapee[query_genes, 'p_val_adj'],
    pval.HC_PC_T.vs.HC_PC_UT=escapee_vs_ut[query_genes, 'p_val_adj']
  );
  write.csv(data_to_output, file=filename, row.names=FALSE, quote=FALSE);
}
generate_deg_table(intersect_upreg, 'intersect_upreg_genes.csv',
                   mean_expr_HC_PC_T, mean_expr_HC_QC_T, mean_expr_HC_PC_UT, escapee_vs_nonescapee, escapee_vs_ut);
generate_deg_table(intersect_downreg, 'intersect_downreg_genes.csv',
                   mean_expr_HC_PC_T, mean_expr_HC_QC_T, mean_expr_HC_PC_UT, escapee_vs_nonescapee, escapee_vs_ut);

# generate table for iRegulon (cytoscape)
write.csv(data.frame(intersect_upreg, intersect_upreg), file='for_iRegulon.csv', quote=FALSE, row.names = FALSE);

# plot genes of interest on tSNE
plot_gene_expr <- function(gene_name, cellranger_tSNE, normalized_raw_expr) {
  expr <- normalized_raw_expr[gene_name, ];
  df_to_plot <- data.frame(TSNE.1=cellranger_tSNE$TSNE.1,
                           TSNE.2=cellranger_tSNE$TSNE.2,
                           logexpr=log2(expr));
  df_to_plot$logexpr[is.infinite(df_to_plot$logexpr)] <- NA;
  sorted_order <- sort.int(expr, index.return = TRUE)$ix;
  df_to_plot <- df_to_plot[sorted_order, ];
  limits <- c(min(df_to_plot$logexpr, na.rm=TRUE), max(df_to_plot$logexpr, na.rm=TRUE));
  p <- ggplot(df_to_plot, aes(x=TSNE.1, y=TSNE.2, col=logexpr)) + geom_point(size = 0.5) +
          scale_colour_gradient2( low = "blue", mid = "grey", high = "red", 
                                  space = "Lab", limits = limits, midpoint = mean(limits),
                                  na.value = "grey85", guide = "colourbar", aesthetics = "colour");
  ggsave(plot=p, filename=paste('logexpr/', gene_name, '.pdf', sep=''), height=3, width=4, useDingbats=FALSE);
}
genes_to_plot <- c('ATF4', 'CDC42EP1', 'RAB32', 'AXL', 'MITF', 'SOX10', 'NGFR', 'CCNA2', 'CCNB1');
for (i in genes_to_plot) {
  plot_gene_expr(i, cellranger_tSNE, normalized_raw_expr);
}

## Step 5. Melanoma phenotype-switching model
plot_phenotype_switch_model <- function(pheno_switch_gene, 
                                        mean_expr_allUT, mean_expr_allT, mean_expr_HC_PC_T, mean_expr_HC_QC_T,
                                        sem_expr_allUT, sem_expr_allT, sem_expr_HC_PC_T, sem_expr_HC_QC_T,
                                        filename) {
  data_to_plot <- data.frame(gene=rep(pheno_switch_gene, 4),
                             expr=c(mean_expr_allUT[pheno_switch_gene], mean_expr_allT[pheno_switch_gene], mean_expr_HC_PC_T[pheno_switch_gene], mean_expr_HC_QC_T[pheno_switch_gene]),
                             sem=c(sem_expr_allUT[pheno_switch_gene], sem_expr_allT[pheno_switch_gene], sem_expr_HC_PC_T[pheno_switch_gene], sem_expr_HC_QC_T[pheno_switch_gene]),
                             cond=rep(c('1AllUT', '2AllT', '3HC_PC_T', '4HC_QC_T'), each=length(pheno_switch_gene)));
  pdf(filename, width=4, height=1.5);
  print(ggplot(data_to_plot, aes(x=gene, y=expr, fill=cond)) +
          geom_bar(stat="identity", position=position_dodge()) +
          scale_fill_brewer(palette="Paired") +
          geom_errorbar(aes(ymin=expr-sem, ymax=expr+sem), width=.4, size=0.3,
                        position=position_dodge(0.9)));
  dev.off();
}
plot_phenotype_switch_model(c('MITF','AXL','SOX10','NGFR'), 
                            mean_expr_allUT, mean_expr_allT, mean_expr_HC_PC_T, mean_expr_HC_QC_T,
                            sem_expr_allUT, sem_expr_allT, sem_expr_HC_PC_T, sem_expr_HC_QC_T,
                            'phenotype_switch_gene.pdf');

# Table S3
query_genes <- c('MITF', 'AXL', 'SOX10', 'NGFR', 'EGFR', 'WNT5A', 'JUN', 'PDGFRB');
data_to_output <- data.frame(
  gene.name=query_genes,
  mean.expr.untreated=mean_expr_allUT[query_genes],
  mean.expr.treated=mean_expr_allT[query_genes],
  mean.expr.HC_PC_T=mean_expr_HC_PC_T[query_genes],
  mean.expr.HC_QC_T=mean_expr_HC_QC_T[query_genes]
);
write.csv(data_to_output, file='TableS3.csv', row.names=FALSE, quote=FALSE);

# save image of workspace
save.image(file='Result.RData');
save(EvsNE_upreg, EvsNE_downreg, EvsUT_upreg, EvsUT_downreg, intersect_upreg, intersect_downreg, file='gene_lists.RData');

# session info
sessionInfo();
