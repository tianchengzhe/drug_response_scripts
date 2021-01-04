# rm(list=ls()), restart R, then execute rmarkdown::render('gene_overlap.R') to generate the html report
# set up working space, load library
setwd('/Users/chti4479/scRNAseq_chen/');
library(GeneOverlap);
set.seed(1);
load(file='gene_lists.RData');

## define gene signatures
# Tirosh et al 2016 Science
tirosh_mitf <- c('MITF', 'TYR', 'PMEL', 'PLP1', 'GPR143', 'MLANA', 'STX7', 'IRF4', 'ERBB3', 'CDH1', 'GPNMB', 'IGSF11', 'SLC24A5', 'SLC45A2', 'RAP2B', 'ASAH1', 'MYO10', 'GRN', 'DOCK10', 'ACSL3', 'SORT1', 'QPCT', 'S100B', 'MYC', 'LZTS1', 'GYG2', 'SDCBP', 'LOXL4', 'ETV5', 'C1orf85', 'HMCN1', 'OSTM1', 'ALDH7A1', 'FOSB', 'RAB38', 'ELOVL2', 'MLPH', 'PLK2', 'CHL1', 'RDH11', 'LINC00473', 'RELL1', 'C21orf91', 'SCAMP3', 'SGK3', 'ABCB5', 'SLC7A5', 'SIRPA', 'WDR91', 'PIGS', 'CYP27A1', 'TM7SF3', 'PTPRZ1', 'CNDP2', 'CTSK', 'BNC2', 'TOB1', 'CELF2', 'ROPN1', 'TMEM98', 'CTSA', 'LIMA1', 'CD99', 'IGSF8', 'FDFT1', 'CPNE3', 'SLC35B4', 'EIF3E', 'TNFRSF14', 'VAT1', 'HPS5', 'CDK2', 'CAPN3', 'SUSD5', 'ADSL', 'PIGY', 'PON2', 'SLC19A1', 'KLF6', 'MAGED1', 'ERGIC3', 'PIR', 'SLC25A5', 'JUN', 'ARPC1B', 'SLC19A2', 'AKR7A2', 'HPGD', 'TBC1D7', 'TFAP2A', 'PTPLAD1', 'SNCA', 'GNPTAB', 'DNAJA4', 'APOE', 'MTMR2', 'ATP6V1B2', 'C16orf62', 'EXOSC4', 'STAM');
tirosh_axl <- c('ANGPTL4', 'FSTL3', 'GPC1', 'TMSB10', 'SH3BGRL3', 'PLAUR', 'NGFR', 'SEC14L2', 'FOSL1', 'SERPINE1', 'IGFBP3', 'TNFRSF12A', 'GBE1', 'AXL', 'PHLDA2', 'MAP1B', 'GEM', 'SLC22A4', 'TYMP', 'TREM1', 'RIN1', 'S100A4', 'COL6A2', 'FAM46A', 'CITED1', 'S100A10', 'UCN2', 'SPHK1', 'TRIML2', 'S100A6', 'TMEM45A', 'CDKN1A', 'UBE2C', 'ERO1L', 'SLC16A6', 'CHI3L1', 'FN1', 'S100A16', 'CRIP1', 'SLC25A37', 'LCN2', 'ENO2', 'PFKFB4', 'SLC16A3', 'DBNDD2', 'LOXL2', 'CFB', 'CADM1', 'LTBP3', 'CD109', 'AIM2', 'TCN1', 'STRA6', 'C9orf89', 'DDR1', 'TBC1D8', 'METTL7B', 'GADD45A', 'UPP1', 'SPATA13', 'GLRX', 'PPFIBP1', 'PMAIP1', 'COL6A1', 'JMJD6', 'CIB1', 'HPCAL1', 'MT2A', 'ZCCHC6', 'IL8', 'TRIM47', 'SESN2', 'PVRL2', 'DRAP1', 'MTHFD2', 'SDC4', 'NNMT', 'PPL', 'TIMP1', 'RHOC', 'GNB2', 'PDXK', 'CTNNA1', 'CD52', 'SLC2A1', 'BACH1', 'ARHGEF2', 'UBE2J1', 'CD82', 'ZYX', 'P4HA2', 'PEA15', 'GLRX2', 'HAPLN3', 'RAB36', 'SOD2', 'ESYT2', 'IL18BP', 'FGFRL1', 'PLEC');

# Rambow et al 2018 Cell
rambow_pigment <- c('APOE', 'EDNRB', 'FABP7', 'GPR143', 'KIT', 'MLANA', 'MLPH', 'PMEL', 'RAB27A', 'SLC24A5', 'SLC45A2', 'SNAI2', 'TRPM1', 'TYR', 'TYRP1');
rambow_smc <- c('NDUFA4L2', 'PGM1', 'B3GNT2', 'TK1', 'GAPDH', 'PRDX6', 'AMD1', 'ACSL3', 'DHCR24', 'LDHA', 'RRM1', 'TYR', 'NDUFA4', 'PKM', 'RPE', 'PHGDH', 'BAAT', 'TUSC3', 'ENPP1', 'ACAT2', 'PGK1', 'ACSS2', 'ALDH1A3', 'TYMS', 'RPN2', 'ALDH1A1', 'LDHB');
rambow_invasion <- c('ADM', 'ANGPTL4', 'AXL', 'BCAT1', 'BGN', 'CCL2', 'CDH13', 'CDH2', 'COL1A1', 'COL1A2', 'COL3A1', 'CYSLTR2', 'DDAH1', 'DLC1', 'DLX1', 'EDNRA', 'ERRFI1', 'FABP4', 'FGF1', 'FOSL2', 'GPC3', 'IGFBP5', 'IGFBP6', 'IL13RA2', 'LMO4', 'LOX', 'LOXL2', 'MGP', 'NDNF', 'NES', 'NR2F1', 'PDGFRB', 'PLXDC1', 'PRDX1', 'PTGER4', 'RGS16', 'RGS5', 'SH2B3', 'SLIT2', 'SOX4', 'SPRY2', 'TGFBI', 'TGM2', 'TMSB4X', 'TNC', 'UNC5B', 'VCAN', 'VEGFA', 'VSNL1');
rambow_ncsc <- c('A2M', 'ADAMTS4', 'ADGB', 'ANXA1', 'AQP1', 'ATP1A2', 'ATP1B2', 'CNN3', 'CADM1', 'COL1A1', 'COL4A1', 'GFRA1', 'GFRA2', 'GFRA3', 'IGF1', 'IL1RAP', 'ITGA1', 'ITGA6', 'L1CAM', 'LAMC1', 'MATN2', 'MPZ', 'NGFR', 'NLGN3', 'NRXN1', 'PDGFB', 'PLAT', 'PRIMA1', 'RSPO3', 'S100A4', 'SEMA3B', 'SLC22A17', 'SLITRK6', 'SYT11', 'THBS2', 'TMEM176B', 'VCAN');

# tsoi et al 2018 Cancer Cell
tsoi_undiff <- c('AJUBA', 'TOR4A', 'MARCH4', 'ZDHHC2', 'ZNF467', 'ZNF185', 'ZIC2', 'VASN', 'UCP2', 'GALNT6', 'TNFAIP2', 'TNFSF18', 'TMEM40', 'TMEM200A', 'TMEM184A', 'TBL1X', 'TRERF1', 'TOX', 'TBC1D2', 'SFN', 'SAMD12', 'SAMD11', 'SOX9', 'SLC8A1', 'SLC38A4', 'SLC16A14', 'SCN5A', 'SCNN1A', 'SH3RF2', 'SERPINB7', 'SLPI', 'SECTM1', 'RUNX2', 'ARHGAP29', 'REN', 'PAWR', 'PSG9', 'PSG5', 'PSG4', 'PBX1', 'PLAGL1', 'PHLDB2', 'PLEKHA6', 'PDGFC', 'PLAU', 'PKP2', 'PLAC8', 'PADI3', 'PITX1', 'NUAK1', 'NTNG1', 'NMT2', 'MYEOV', 'MICAL2', 'MGST1', 'MECOM', 'LYPD6B', 'LAMA5', 'KISS1', 'KRT86', 'KRT81', 'KRT80', 'KRT8', 'KRT7', 'KRT18', 'JUP', 'IL7R', 'IL4R', 'IRS1', 'IGFN1', 'HES7', 'GDA', 'GLIS2', 'GATA2', 'GPRC5C', 'GPRC5A', 'FMNL1', 'FOXA1', 'FLNC', 'FERMT1', 'FAT4', 'FAM196B', 'ELFN2', 'EGFR', 'DSE', 'DMBT1', 'DIO2', 'DOCK2', 'CYP2S1', 'CRIM1', 'CDK15', 'CORO6', 'COLEC10', 'CCDC88C', 'CCDC69', 'F3', 'F2RL1', 'CLU', 'CDYL2', 'CITED2', 'CARD11', 'CPA4', 'CREB3L1', 'CNN1', 'CALB2', 'CDH4', 'BTBD11', 'BDNF', 'BASP1', 'BNC1', 'ATP8B1', 'ABCG2', 'ARMC4', 'ANKRD1', 'AR', 'AMIGO2', 'ADAMTSL1', 'ACSL5');
tsoi_undiff_ncsc <- c('VIT', 'VIPR1', 'VEGFC', 'TWIST2', 'TNFRSF12A', 'TPM1', 'TPBG', 'TLE4', 'TOX2', 'TLR4', 'THSD4', 'STX1A', 'SYT1', 'SYNPO', 'STRA6', 'STC2', 'SPRED3', 'SPOCD1', 'SPOCK1', 'SLC2A1', 'SLC16A2', 'SLC14A1', 'SLC12A8', 'SMAGP', 'SLIT2', 'SDK1', 'STAC', 'SLFN11', 'S100A2', 'ROBO4', 'RAB27B', 'PKIA', 'PRSS23', 'PAPPA', 'PRDM1', 'KCNMA1', 'KCNN4', 'PODXL', 'PDGFRB', 'PLAUR', 'PXDN', 'PTX3', 'NMNAT2', 'NRP1', 'NGEF', 'NEGR1', 'NRG1', 'NTN4', 'MT2A', 'MT1E', 'MPP4', 'LOXL2', 'LDOC1', 'LAMB3', 'JUN', 'IL31RA', 'IL11', 'IL1B', 'ITGA3', 'ITGA2', 'IGFBP6', 'ID1', 'INHBA', 'HRH1', 'GAS6', 'GLIPR1', 'GFRA1', 'GATA3', 'GPR176', 'FZD2', 'FJX1', 'FOSL1', 'FOXF1', 'FBLIM1', 'FLNB', 'FAM83G', 'FAM20C', 'FAM171A1', 'FAM155A', 'ERRFI1', 'EFNB2', 'DPYD', 'DKK1', 'DOCK5', 'CYR61', 'CLMP', 'COL13A1', 'COL12A1', 'COL5A1', 'F2RL2', 'C16orf45', 'C15orf52', 'C12orf75', 'CD163L1', 'CAV1', 'CARD10', 'CLCF1', 'CDH13', 'BMP2', 'AXL', 'ABCC3', 'ARNTL2', 'ANTXR2', 'ANXA1', 'AKR1C3', 'ARL4C');
tsoi_ncsc <- c('PXYLP1', 'CXCL8', 'CEMIP', 'TCAF2', 'ZNF469', 'WNT5A', 'TMEM47', 'TMEM171', 'TGFBI', 'TGFA', 'TFAP2C', 'TSPAN13', 'SQRDL', 'SULF1', 'ST8SIA5', 'SOX2', 'SLC24A3', 'SLITRK6', 'SHISA2', 'SH3PXD2A', 'SERTAD4', 'STK32B', 'SEMA3B', 'SFRP1', 'S100A6', 'RAMP1', 'PMEPA1', 'PCSK5', 'PHLDA2', 'PLA2G7', 'OPRD1', 'NTM', 'NRXN3', 'NES', 'MUC5B', 'MAP1LC3A', 'LRRC15', 'KIAA1755', 'ITGB8', 'IER3', 'HHEX', 'GDNF', 'GLI2', 'FOXC2', 'FLT1', 'FAT3', 'FEZ1', 'FAM135B', 'EHF', 'EML1', 'DRD2', 'DEPDC7', 'CYB5R2', 'CSRP2', 'CCL2', 'CADM3', 'CADM1', 'CD96', 'CTSS', 'CHST2', 'CHST1', 'CACNA2D3', 'BST1', 'ABCA6', 'ANGPTL4', 'AIM2');
tsoi_ncsc_tran <- c('SPRY4', 'SORCS1', 'SLC35F1', 'SERPINA5', 'RFTN2', 'PCDH1', 'PTPRZ1', 'PRICKLE2', 'OLIG2', 'LOXL4', 'LOXL3', 'LGI4', 'LAMA4', 'GAS7', 'GRIK2', 'FREM2', 'FREM1', 'EPHB3', 'CRIP2', 'COL4A1', 'CADM4', 'BAALC', 'ABCA8', 'AGMO', 'ALDH1A3');
tsoi_tran <- c('XYLT1', 'TSPAN7', 'SOD3', 'SCRG1', 'SORL1', 'SEMA3E', 'SELENBP1', 'RNASE1', 'RAPGEF4', 'PCDH7', 'PRSS33', 'PCSK6', 'PLBD1', 'NELL1', 'NPR1', 'MCAM', 'MMP15', 'MAMDC2', 'LSAMP', 'LRRTM4', 'GDF11', 'FXYD3', 'EBF3', 'COL11A2', 'COL9A1', 'CX3CL1', 'BCHE', 'ANO4', 'ALDH1A1');
tsoi_tran_mela <- c('ADGRG1', 'MOB3B', 'SEPT4', 'TUBB4A', 'UBAP1L', 'ZNF704', 'WFDC1', 'VGF', 'VAT1', 'GALNT3', 'UGT2B7', 'TYRP1', 'TYR', 'TTYH2', 'TMC6', 'TMCC2', 'TBC1D7', 'TBC1D16', 'STXBP6', 'ST8SIA1', 'ST3GAL6', 'SOX6', 'SLC5A4', 'SLC45A2', 'SLC27A3', 'SLC24A5', 'SIRPA', 'SCUBE2', 'STK32A', 'RLBP1', 'RENBP', 'RRAGD', 'RASSF3', 'RAP1GAP', 'RAB38', 'QDPR', 'P2RX7', 'PRR5', 'PMEL', 'PLXNC1', 'PLEKHH1', 'PLA1A', 'PDE3B', 'PHACTR1', 'PPARGC1A', 'PMP2', 'PI15', 'OGDHL', 'NRG3', 'NKAIN4', 'ASAH1', 'NAT8L', 'GNPTAB', 'MYO10', 'MBP', 'MCC', 'MITF', 'MFAP3L', 'LDB3', 'LRGUK', 'LGI3', 'LINGO1', 'LGALS3', 'LAMC3', 'LAMA1', 'KLF15', 'KAZN', 'IRX6', 'IRF4', 'INPP4B', 'ID4', 'IGSF11', 'HAS2', 'HPS4', 'GREB1', 'GHR', 'GDF15', 'GAB2', 'GPM6B', 'GPNMB', 'GYPC', 'GYG2', 'GAPDHS', 'GJB1', 'GPRC5B', 'FMN1', 'FCGR2B', 'FCER1G', 'FAM69C', 'FAM167B', 'ESRP1', 'DUSP15', 'DSTYK', 'DCT', 'D4S234E', 'DAPK1', 'CDK5R1', 'CELF2', 'CTTNBP2', 'CHCHD6', 'CHCHD10', 'C11orf96', 'CHN2', 'CHL1', 'CITED1', 'CARD14', 'CPN1', 'CA14', 'CAPN3', 'MERTK', 'BCAS3', 'BEST1', 'BCL2A1', 'BIRC7', 'ATP6V0A4', 'ATP10A', 'APOE', 'APOC1', 'ASB2', 'ANK2', 'ADRBK2', 'ADCY1', 'ACP5', 'PFKFB2', 'HTR2B');
tsoi_mela <- c('CCDC171', 'CFAP61', 'ZDHHC11B', 'VEPH1', 'TNFRSF14', 'TDRD3', 'TPPP', 'TRIM63', 'TRPM1', 'TTC39A', 'TSPAN10', 'SLC7A8', 'SLC16A6', 'SLAMF7', 'SEMA6A', 'RUNX3', 'RNF144B', 'RNLS', 'RGS12', 'PYCARD', 'PRUNE2', 'PRKCB', 'PRDM7', 'KCNAB2', 'OCA2', 'NR4A3', 'NAV2', 'MYO1D', 'MAPK4', 'MAT1A', 'MLANA', 'LXN', 'KCP', 'IL16', 'IL12RB2', 'HSD17B14', 'HMOX1', 'H2AFJ', 'GOLGA7B', 'QPCT', 'GFOD1', 'GPR143', 'FYB', 'FAM83H', 'FAM174B', 'EPHA5', 'ENTHD1', 'DNAJA4', 'DENND2D', 'C2orf88', 'CCL18', 'CEACAM1', 'CAPG', 'CDH3', 'CDH1', 'ATP6V0D2', 'ABCD1', 'ABCB5', 'APOLD1', 'ANKRD30B', 'ADCY2', 'ADAM23');

# organize gene lists 
all_external_gene_lists <- list(tirosh_mitf, tirosh_axl, 
                                rambow_pigment, rambow_smc, rambow_invasion, rambow_ncsc,
                                tsoi_undiff, tsoi_undiff_ncsc, tsoi_ncsc, tsoi_ncsc_tran,
                                tsoi_tran, tsoi_tran_mela, tsoi_mela);
num_external_genes <- sapply(all_external_gene_lists, length);
name_external_genes <- c('tirosh_mitf', 'tirosh_axl', 
                         'rambow_pigment', 'rambow_smc', 'rambow_invasion', 'rambow_ncsc',
                         'tsoi_undiff', 'tsoi_undiff_ncsc', 'tsoi_ncsc', 'tsoi_ncsc_tran',
                         'tsoi_tran', 'tsoi_tran_mela', 'tsoi_mela');
print(data.frame(gene_list=name_external_genes, gene_counts=num_external_genes));

## organize gene lists in this work
all_curr_gene_lists <- list(EvsNE_upreg, EvsNE_downreg, c(EvsNE_upreg, EvsNE_downreg), intersect_upreg);
num_curr_genes <- sapply(all_curr_gene_lists, length);
name_curr_genes <- c('up_escapee', 'up_nonescapee', 'up_either', 'up_both');
print(data.frame(gene_list=name_curr_genes, gene_counts=num_curr_genes));

## define gene overlap function
disp_gene_overlap <- function(gene_sig1, gene_sig2, print_msg=NULL) {
  res <- testGeneOverlap(newGeneOverlap(gene_sig1, gene_sig2, spec = "hg19.gene"));
  if (!is.null(print_msg)) {
    print(print_msg);
    print(paste('Intersect Genes: ', paste(res@intersection, collapse=', '), sep=''));
    print(paste('pval=', res@pval, sep=''));
  } else {
    return(res@pval);
  }
}

## Tirosh, need details to draw venn diagram
disp_gene_overlap(EvsNE_upreg, tirosh_axl, 'E vs AXL');
disp_gene_overlap(EvsNE_upreg, tirosh_mitf, 'E vs MITF');
disp_gene_overlap(EvsNE_downreg, tirosh_axl, 'NE vs AXL');
disp_gene_overlap(EvsNE_downreg, tirosh_mitf, 'NE vs MITF');

## Others, only need p-val
all_pval <- matrix(NA, length(all_external_gene_lists), length(all_curr_gene_lists));
rownames(all_pval) <- name_external_genes;
colnames(all_pval) <- name_curr_genes;
for (i in 1:length(all_external_gene_lists)) {
  for (j in 1:length(all_curr_gene_lists)) {
    all_pval[i, j] <- disp_gene_overlap(all_external_gene_lists[[i]], all_curr_gene_lists[[j]]);
  }
}
print(all_pval);

sessionInfo();
