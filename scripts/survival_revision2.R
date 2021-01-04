# rm(list=ls()), restart R, then execute rmarkdown::render('survival_revision2.R') to generate the html report
# set up working space, load library
setwd('/Users/chti4479/scRNAseq_chen/');
library(survival); library(survminer); library(ggplot2);
set.seed(1);

# oncolnc patients: exactly the patients with sex, age, days, status, and expr levels. so should use them

# load data
oncolnc_info <- read.delim('data/oncolnc_patient_id');
oncolnc_info$Patient <- as.character(oncolnc_info$Patient);
oncolnc_info$Days <- as.numeric(as.character(oncolnc_info$Days));
oncolnc_info$Status <- as.character(oncolnc_info$Status);

patient_info <- read.delim('data/nationwidechildrens.org_clinical_patient_skcm.txt', header = TRUE)[-c(1:2), ];
patient_info$bcr_patient_barcode <- as.character(patient_info$bcr_patient_barcode);
patient_info$gender <- as.character(patient_info$gender);
patient_info$age_at_diagnosis <- as.numeric(as.character(patient_info$age_at_diagnosis));

expr_data <- read.delim('data/data_RNA_Seq_v2_expression_median.txt');
expr_data <- expr_data[, !endsWith(colnames(expr_data), '.07')];
expr_data <- expr_data[, !(colnames(expr_data) %in% c('TCGA.ER.A19T.06', 'TCGA.ER.A2NF.06'))];

# extract information
patient_gender <- array(NA, nrow(oncolnc_info));
patient_age <- array(NA, nrow(oncolnc_info));
patient_expr_col <- array(NA, nrow(oncolnc_info));
patient_stage <- array(NA, nrow(oncolnc_info));
for (i in 1:nrow(oncolnc_info)) {
  id <- which(patient_info$bcr_patient_barcode ==oncolnc_info$Patient[i]);
  patient_gender[i] <- patient_info$gender[id];
  patient_age[i] <- patient_info$age_at_diagnosis[id];
  
  id <- which(startsWith(colnames(expr_data), gsub('-', '.', oncolnc_info$Patient[i])));
  patient_expr_col[i] <- id;
  if (as.numeric(substr(colnames(expr_data)[id], 15, 15)) == 1) {
    patient_stage[i] <- 'Primary';
  } else {
    patient_stage[i] <- 'Metastatic';
  }
}
expr_data <- expr_data[, c(1, 2, patient_expr_col)];
oncolnc_info <- data.frame(Patient=oncolnc_info$Patient,
                           Days=oncolnc_info$Days,
                           Status=ifelse(oncolnc_info$Status=='Alive', 0, 1),
                           Gender=as.factor(patient_gender),
                           Age=patient_age,
                           Stage=as.factor(patient_stage));

# calculate cox coefficient for all genes
cox_coef <- matrix(NA, nrow=nrow(expr_data), 3); cox_coef[, 1] <- expr_data[, 2]; 
for (i in 1:nrow(expr_data)) {
  data <- cbind(oncolnc_info, data.frame(Gene=scale(as.numeric(expr_data[i, -c(1:2)]))));
  res.cox <- tryCatch(summary(coxph(Surv(Days, Status) ~ Gender + Age + Stage + Gene, data=data))$coefficients,
                      error=function(e) matrix(NA, nrow=4, ncol=5));
  cox_coef[i, 2] <- res.cox[4, 1]; # positive: bad survival; negative: good survival
  cox_coef[i, 3] <- res.cox[4, 5];
}

# all genes:
print(mean(cox_coef[, 3] < 0.05 & cox_coef[, 2] > 0, na.rm=TRUE)); # fraction of genes with bad survival
print(mean(cox_coef[, 3] < 0.05 & cox_coef[, 2] < 0, na.rm=TRUE)); # fraction of genes with good survival

# baseline KM
all_pval <- array(NA, nrow(expr_data)); all_median_high <- array(NA, nrow(expr_data)); all_median_low <- array(NA, nrow(expr_data));
num_sample_include <- floor(nrow(oncolnc_info)/4);
for (i in 1:nrow(expr_data)) {
  if ((i %% 100) == 0) {
    print(i);
  }
  low_id <- sort.int(as.numeric(expr_data[i, -c(1,2)]), decreasing = FALSE, index.return = TRUE)$ix[1:num_sample_include];
  high_id <- sort.int(as.numeric(expr_data[i, -c(1,2)]), decreasing = TRUE, index.return = TRUE)$ix[1:num_sample_include];
  data <- data.frame(time=oncolnc_info$Days[c(low_id, high_id)],
                     vital_status=oncolnc_info$Status[c(low_id, high_id)],
                     expr_level=c(rep('Low', num_sample_include), rep('High', num_sample_include)));
  sfit <- surv_fit(Surv(time, vital_status)~expr_level, data=data);
  all_pval[i] <- surv_pvalue(sfit)$pval;
  temp <- surv_median(sfit)$median; all_median_high[i] <- temp[1]; all_median_low[i] <- temp[2];
}

# calculate null probability
pos_prob <- mean(cox_coef[, 3] < 0.05 & cox_coef[, 2] > 0 & all_pval < 0.05 & all_median_high < all_median_low, na.rm=TRUE);
neg_prob <- mean(cox_coef[, 3] < 0.05 & cox_coef[, 2] < 0, na.rm=TRUE);
print(pos_prob); 
print(neg_prob);
null_prob <- sum(sapply(5:38, function(x) {
  choose(38, x) * (pos_prob**x) * ((1-pos_prob-neg_prob)**(38-x));
}));
print(null_prob);

save(cox_coef, all_pval, all_median_high, all_median_low, pos_prob, neg_prob, null_prob, file='survival_revision2.RData');

sessionInfo();
