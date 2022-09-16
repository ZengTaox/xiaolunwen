rm(list = ls())
library(ROCR)
library(ggplot2)
library(limma)
options("expressions" = 20000)
memory.limit(size = 8000000)

alldata <- readRDS('./all_sample.Rds')

patient_data <- alldata[["patient_data"]]
log2_TPM <- alldata$log2_TPM

v <- apply(log2_TPM, 1, mean)
log2_TPM <- log2_TPM[v > 1,]

noTR <- patient_data$type != 'TR' & !is.na(patient_data$type)
noTR_patient_data <- patient_data[noTR,]
noTR_log2_TPM <- log2_TPM[,noTR]

noTR_log2_TPM <- as.data.frame(noTR_log2_TPM)
length(noTR_patient_data$type)
dim(noTR_log2_TPM)

group <- noTR_patient_data$type
mm <- model.matrix( ~ 0 + group)
fit <- lmFit(noTR_log2_TPM, mm)
head(coef(fit))
contr <- makeContrasts(groupTP - groupNT, levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table <- topTable(tmp, sort.by = "P", n = Inf)
# write.csv(top.table, file = './output/TCGA_volcano.csv')
# ===================================== cox ====================================
library("survival")
library("survminer")
library(limma)
library(ROCR)

data0 <- noTR_log2_TPM[,group == 'TP']
data0 <- apply(data0, 1,scale)
dim(data0)
data0 <- as.data.frame(data0)

noNT_noTR_patient_data <- noTR_patient_data[group == 'TP',]
data0$time <- noNT_noTR_patient_data$days
data0$status <- 0
data0$status[noNT_noTR_patient_data$status == 'Dead'] <- 1
data0 <- na.omit(data0)

data0$status
data0$time
length(data0)
dim(data0)
colnames(data0)[1:19166]
sv_gl <- list()
for(i in 1:19166){
  datat <- data0[,c(i,19167,19168)]
  res.cox <- coxph(Surv(time, status) ~ ., data =  datat)
  tmp0 <- cox.zph(res.cox)
  global_zph_p <- cox.zph(res.cox)$table[dim(cox.zph(res.cox)$table)[1],
                                         dim(cox.zph(res.cox)$table)[2]]
  tmp <- summary(res.cox)
  coef <- tmp$coefficients[1]
  exp_coef <- tmp$coefficients[2]
  se <- tmp$coefficients[3]
  p_sv <- tmp$coefficients[5]
  
  # rms::vif(res.cox)
  
  C_index <- tmp[["concordance"]][1]
  C_index_se <- tmp[["concordance"]][2]
  
  sv_gl[[colnames(data0)[i]]] <- c(p_sv,global_zph_p, coef,
                                   exp_coef, se, 
                                   C_index,C_index_se)
}

sv_gld <- t(list2DF(sv_gl))
colnames(sv_gld) <- c('p_sv','zph_p','coef','exp_coef',
                      'se','C_index', 'C_index_se')
sv_gld <- as.data.frame(sv_gld)
# saveRDS(sv_gld,file = './output/sv_cox.Rds')
# ============================== 开始建模 ========================================
rm(list = ls())
library(ggplot2)
library(ggsci)
library(glmnet)
library(survival)
library(ROCR)
library(tidyverse)
library(ComplexHeatmap)
library(clusterProfiler)
library("survival")
library("survminer")
library(timeROC)
library(limma)

alldata <- readRDS('all_sample.Rds')

patient_data <- alldata[["patient_data"]]
log2_TPM <- alldata$log2_TPM

v <- apply(log2_TPM, 1, mean)
log2_TPM <- log2_TPM[v > 1,]

noTR <- patient_data$type != 'TR' & !is.na(patient_data$type)
noTR_patient_data <- patient_data[noTR,]
noTR_log2_TPM <- log2_TPM[,noTR]
noTR_log2_TPM <- as.data.frame(noTR_log2_TPM)

data0 <- noTR_log2_TPM[,noTR_patient_data$type == 'TP']
data0 <- apply(data0, 1,scale)
data0 <- as.data.frame(data0)

noNT_noTR_patient_data <- noTR_patient_data[noTR_patient_data$type == 'TP',]
data0$time <- noNT_noTR_patient_data$days
data0$status <- 0
data0$status[noNT_noTR_patient_data$status == 'Dead'] <- 1
data0$barcode <- noNT_noTR_patient_data$patient
data0 <- na.omit(data0)

sv_gld <- readRDS('./sv_cox.Rds')
# write.csv(sv_gld,file = 'singleCox_svgld.csv')
top.table <- read.csv('./TCGA_volcano.csv')
rownames(top.table) <- top.table$X

gndff <- rownames(top.table)[top.table$adj.P.Val < 0.05 & 
                               abs(top.table$logFC) > 0.5]
select <- intersect(gndff, 
                    rownames(sv_gld)[sv_gld$zph_p > 0.4 & sv_gld$p_sv < 0.05 & 
                                       sv_gld$C_index > 0.58]
)
select <- select[-grep(select, pattern = '(^AC[0-9]|^AL[0-9])', 
                       value = F, 
                       ignore.case = T)]
top.table$`-logP` <- -log10(top.table$adj.P.Val)
top.table$group <- 'No-Diff'
top.table$group[top.table$X %in% select] <- 'Diff_Prognostic'
top.table$group[top.table$X %in% gndff & (!top.table$X %in% select)] <- 'Diff_Non_Prognostic'

# vo <- ggplot(top.table, aes(x = logFC, y = `-logP`, color = group,
#                             group = group)) +
#   geom_point() +
#   geom_point(top.table[top.table$group == 'Diff_Prognostic',],
#              mapping = aes(x = logFC, y = `-logP`, color = group, group = group)) +
#   theme_bw() + scale_color_jco()
# vo
# ggsave(vo, filename = 'volcano.pdf',width = 7, height = 7)
# ============================= COX_glmn =======================================
data0 <- data0[data0[,'time'] > 0,]
cox_data <- data0[,select]

train_y <- data0[,c(19167,19168)]
colnames(train_y) <- c('time', 'status')

cox_data <- as.matrix(cox_data)
train_y <- as.matrix(train_y)

cox_data <- cox_data[train_y[,'time'] > 0,]

train_y <- train_y[train_y[,'time'] > 0,]

fit <- cv.glmnet(cox_data, train_y, family = "cox")
# plot(fit)

# pdf('模型超参数.pdf')
# plot(fit)
# plot(fit, xvar = "dev", label = TRUE)
# plot(fit, xvar = "lambda", label = TRUE)
# plot(fit, xvar = "norm", label = TRUE)
# dev.off()

# L0 <- fit$lambda.min
L0 <- 2.7^(-3.5)
# L0 <- 0.056
L0
# 0.06
fit0 <- glmnet(cox_data, train_y, family = "cox", lambda = L0)
t <- coef(fit0)
rownames(t)[t[,1] != 0]
# [1] "PLEK2"     "PTPRH"    
# [3] "OGFRP1"    "CHRNA5"   
# [5] "CBFA2T3"   "SMIM15"   
# [7] "AVEN"      "MELTF"    
# [9] "KRT8"      "RGS20"    
# [11] "FAM207A"   "SOWAHC"   
# [13] "ELF5"      "LSP1P4"   
# [15] "C20orf197" "C11orf16" 
# [17] "DNALI1"
cox_data0 <- cox_data[,rownames(t)[t[,1] != 0]]
lung <- as.data.frame(cbind(cox_data0, train_y))
res.cox <- coxph(Surv(time, status) ~ ., data =  lung)
# cox.zph(res.cox)$table
# summary(res.cox)
# rms::vif(res.cox)
# PLEK2     PTPRH    OGFRP1    CHRNA5   CBFA2T3    SMIM15      AVEN     MELTF 
# 1.924750  1.470266  1.193339  1.309820  1.279034  1.252861  1.493640  1.815486 
# KRT8     RGS20   FAM207A    SOWAHC      ELF5    LSP1P4 C20orf197  C11orf16 
# 1.312085  1.559300  1.277998  1.332380  1.279083  1.232585  1.413043  1.344663 
# DNALI1 
# 1.594311 
# > summary(res.cox)
# Call:
#   coxph(formula = Surv(time, status) ~ ., data = lung)
# 
# n= 520, number of events= 187 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)  
# PLEK2      0.08023   1.08354  0.10568  0.759   0.4477  
# PTPRH      0.08992   1.09409  0.08929  1.007   0.3139  
# OGFRP1     0.15723   1.17027  0.07651  2.055   0.0399 *
#   CHRNA5     0.19302   1.21291  0.08414  2.294   0.0218 *
#   CBFA2T3   -0.03516   0.96545  0.08641 -0.407   0.6841  
# SMIM15     0.07812   1.08125  0.07996  0.977   0.3286  
# AVEN       0.10300   1.10849  0.08190  1.258   0.2085  
# MELTF      0.12072   1.12831  0.10183  1.185   0.2358  
# KRT8       0.06052   1.06239  0.09133  0.663   0.5076  
# RGS20      0.03030   1.03076  0.08213  0.369   0.7122  
# FAM207A    0.12458   1.13268  0.09576  1.301   0.1933  
# SOWAHC     0.11292   1.11955  0.08147  1.386   0.1657  
# ELF5      -0.08243   0.92088  0.08825 -0.934   0.3503  
# LSP1P4     0.14295   1.15368  0.07854  1.820   0.0687 .
# C20orf197 -0.07075   0.93170  0.09542 -0.741   0.4585  
# C11orf16  -0.11591   0.89056  0.09411 -1.232   0.2181  
# DNALI1    -0.03104   0.96944  0.08447 -0.367   0.7133  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# PLEK2        1.0835     0.9229    0.8808     1.333
# PTPRH        1.0941     0.9140    0.9184     1.303
# OGFRP1       1.1703     0.8545    1.0073     1.360
# CHRNA5       1.2129     0.8245    1.0285     1.430
# CBFA2T3      0.9654     1.0358    0.8150     1.144
# SMIM15       1.0813     0.9249    0.9244     1.265
# AVEN         1.1085     0.9021    0.9441     1.302
# MELTF        1.1283     0.8863    0.9242     1.378
# KRT8         1.0624     0.9413    0.8883     1.271
# RGS20        1.0308     0.9702    0.8775     1.211
# FAM207A      1.1327     0.8829    0.9388     1.367
# SOWAHC       1.1195     0.8932    0.9543     1.313
# ELF5         0.9209     1.0859    0.7746     1.095
# LSP1P4       1.1537     0.8668    0.9891     1.346
# C20orf197    0.9317     1.0733    0.7728     1.123
# C11orf16     0.8906     1.1229    0.7405     1.071
# DNALI1       0.9694     1.0315    0.8215     1.144
# 
# Concordance= 0.679  (se = 0.023 )
# Likelihood ratio test= 88.73  on 17 df,   p=1e-11
# Wald test            = 86.91  on 17 df,   p=2e-11
# Score (logrank) test = 92.53  on 17 df,   p=2e-12

mld <- names(rms::vif(res.cox))
sv_gld$Upper <- sv_gld$coef + 1.96 * sv_gld$se
sv_gld$Lower <- sv_gld$coef - 1.96 * sv_gld$se
sv_gld$Varnames <- rownames(sv_gld)

sv_gld0 <- sv_gld[sv_gld$Varnames %in% mld,]
p <- ggplot(sv_gld0, aes(x = coef, y = Varnames,
                         group = Varnames, col = p_sv)) # 不同形状shape= Factor

forsv <- p + geom_point(size=3) +
  
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.3,size=2) +
  
  scale_x_continuous(limits= c(-0.8, 0.8), breaks= seq(-1, 1, 0.2)) +
  
  geom_vline(aes(xintercept = 0)) +
  
  xlab('Coef') + ylab(' ') + theme_bw() + scale_color_viridis_c()
forsv
# ggsave(forsv, filename = 'sv.pdf',width = 5,height = 5)
# saveRDS(res.cox, file = 'res.cox_model.Rds')
res.cox <- readRDS('res.cox_model.Rds')
lac <- lung
ROCR.simple <- predict(res.cox, lac)
names(ROCR.simple) <- data0$barcode
# saveRDS(ROCR.simple, file = 'TCGA_HR.Rds')
# saveRDS(lung, file = 'TCGA_modle_lung.Rds')

mx <- t(cox_data0)
od <- order(ROCR.simple, decreasing = F)
mx <- mx[,od]
ROCR.simple0 <- ROCR.simple[od]
vvv <- rep('high', times = length(ROCR.simple))
vvv[ROCR.simple0 < 0] <- 'low'
colnames(mx) <- NULL
an <- HeatmapAnnotation(group = vvv,
                        col = list(
                          group = c('low' = "blue",
                                    "high" = 'red')
                        ))
# pdf('TCGA_heatmap.pdf', width = 7,height = 5)
# Heatmap(mx, cluster_columns = F, top_annotation = an)
# dev.off()
rownames(mx)
# [1] "PLEK2"     "PTPRH"     "OGFRP1"    "CHRNA5"   
# [5] "CBFA2T3"   "SMIM15"    "AVEN"      "MELTF"    
# [9] "KRT8"      "RGS20"     "FAM207A"   "SOWAHC"   
# [13] "ELF5"      "LSP1P4"    "C20orf197" "C11orf16" 
# [17] "DNALI1" 

ddd <- data.frame(HR = ROCR.simple,
                  time = train_y[,'time'],
                  status = train_y[,'status'])
ddd <- ddd[order(ddd$HR, decreasing = F),]
ddd$rank <- 1:dim(ddd)[1]
ddd$group <- rep('high', times = length(ROCR.simple))
ddd$group[ddd$HR < 0] <- 'low'

ppp <- ggplot(ddd, aes(x = rank, y = HR, color = group)) + 
  geom_point() + theme_bw() + scale_color_manual(
    values = c('red','blue')
  )
ppp
# ggsave(ppp, filename = 'rank_HR.pdf', width = 7,height = 7)

ddd$status <- factor(ddd$status, levels = c(0,1))
ddd$stat <- 'Alive'
ddd$stat[ddd$status == 1] <- 'Dead'
ppp <- ggplot(ddd, aes(x = rank, y = time, color = stat)) + 
  geom_point() + theme_bw() + scale_color_manual(
    values = c('blue','red')
  )
ppp
# ggsave(ppp, filename = 'rank_Dead.pdf', width = 7, height = 7)

ddd$status <- as.numeric(as.character(ddd$status))

fit <- survfit(Surv(time, status) ~ group, data = ddd)
res.sum <- surv_summary(fit)
# pdf('sv1.pdf')
# ggsurvplot(fit,
#            pval = TRUE, conf.int = F,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            # linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("red", "blue"),
#            xlim = c(0, 5000)
# )
# ggsurvplot(fit, 
#            conf.int = F, # 增加置信区间
#            fun = "cumhaz",
#            palette = c("red", "blue"),
#            xlim = c(0, 5000)) # 绘制累计风险曲线
# dev.off()

# PCA
# mxpca <- prcomp(t(mx))
# ddd <- data.frame(
#   HR = ROCR.simple,
#   PC1 = mxpca$x[,1],
#   PC2 = mxpca$x[,2]
# )
# ddd$lable <- 'high'
# ddd$lable[ddd$HR < 0] <- 'low'
# ppca <- ggplot(ddd,aes(x = PC1, y = PC2, color = lable)) + geom_point() + theme_bw()
# ppca
# ggsave(ppca, filename = 'TCGA_PCA.pdf')
alldata <- readRDS('./all_sample.Rds')
patient_data <- alldata[["patient_data"]]
log2_TPM <- alldata$log2_TPM
v <- apply(log2_TPM, 1, mean)
log2_TPM <- log2_TPM[v > 1,]
noTR <- patient_data$type != 'TR' & !is.na(patient_data$type)
noTR_patient_data <- patient_data[noTR,]
noTR_log2_TPM <- log2_TPM[,noTR]
noTR_log2_TPM <- as.data.frame(noTR_log2_TPM)

data0 <- t(noTR_log2_TPM[,noTR_patient_data$type == 'TP'])
data0 <- as.data.frame(data0)
dim(data0)
noNT_noTR_patient_data <- noTR_patient_data[noTR_patient_data$type == 'TP',]
data0$time <- noNT_noTR_patient_data$days
data0$status <- 0
data0$status[noNT_noTR_patient_data$status == 'Dead'] <- 1
data0 <- na.omit(data0)
data0 <- data0[data0[,'time'] > 0,]
genes <- t(data0[,1:19166])
dim(genes)

group <- rep('high', times = length(ROCR.simple))
group[ROCR.simple < 0] <- 'low'
mm <- model.matrix( ~ 0 + group)
fit <- lmFit(genes, mm)
head(coef(fit))

contr <- makeContrasts(grouphigh - grouplow, 
                       levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_HL <- topTable(tmp, sort.by = "P", n = Inf)
top.table_HL$`-logP` <- - log10(top.table_HL$P.Value)
top.table_HL$group <- 'Sig-UP'
top.table_HL$group[top.table_HL$logFC < 0.5 & 
                     top.table_HL$logFC > -0.5 |
                     top.table_HL$adj.P.Val > 0.05] <- 'No-Sig'
top.table_HL$group[top.table_HL$logFC < -0.5 & 
                     top.table_HL$adj.P.Val < 0.05] <- 'Sig-DN'


set <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt')
v <- top.table_HL$logFC
names(v) <- rownames(top.table_HL)
v <- sort(v, decreasing = T)
gseaHL <- GSEA(v, TERM2GENE = set)
gseaHLd <- gseaHL@result
gseaHLd <- gseaHLd[order(gseaHLd$NES, decreasing = F),]
gseaHLd$ID <- factor(gseaHLd$ID, levels = gseaHLd$ID)
gseaHLd$color <- 'orange'
gseaHLd$color[gseaHLd$NES < 0] <- 'green'

library(viridis)
# write.csv(gseaHLd,file = 'TCAG_GSEA.csv')
gseaHLd <- read.csv('TCAG_GSEA.csv')
gseaHLd$ID <- factor(gseaHLd$ID, levels = gseaHLd$ID)
gseaHLd$absNES <- abs(gseaHLd$NES)
gseap <- ggplot(gseaHLd, aes(x = ID, y = NES,
                             fill = absNES)) +
  geom_bar(position = 'dodge', stat = 'identity', color = 'black') +
  coord_flip() + theme_bw() + 
  scale_fill_viridis(option="magma",
                     direction = -1,
                     begin = 0.35,
                     end = 0.9)
gseap
ggsave(gseap, filename = 'gseapviridis.pdf', width = 10, height = 7)

top.table_HL$group <- factor(top.table_HL$group, levels = c('Sig-DN',
                                                            'Sig-UP',
                                                            'No-Sig'))
vo <- ggplot(top.table_HL, aes(x = logFC, y = `-logP`, color = group, 
                               group = group)) + 
  geom_point() +
  theme_bw() + scale_color_jco()
vo
# ggsave(vo, filename = 'HL_volcano.pdf',width = 7, height = 7)

upg <- rownames(top.table_HL)[top.table_HL$group == 'Sig-UP']
dng <- rownames(top.table_HL)[top.table_HL$group == 'Sig-DN']

listd <- list()
setGObp <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt')
upgebp <- enricher(upg,TERM2GENE = setGObp)
dngebp <- enricher(dng,TERM2GENE = setGObp)
listd[['bp']] <- list(
  upgebp,
  dngebp)

setGOcc <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.cc.v7.4.symbols.gmt')
upgecc <- enricher(upg,TERM2GENE = setGOcc)
dngecc <- enricher(dng,TERM2GENE = setGOcc)
listd[['cc']] <- list(
  upgecc,
  dngecc)

setGOmf <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.mf.v7.4.symbols.gmt')
upgemf <- enricher(upg,TERM2GENE = setGOmf)
dngemf <- enricher(dng,TERM2GENE = setGOmf)
listd[['mf']] <- list(
  upgemf,
  dngemf)

setKEGG <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c2.cp.kegg.v7.4.symbols.gmt')
upgesetKEGG <- enricher(upg,TERM2GENE = setKEGG)
dngesetKEGG <- enricher(dng,TERM2GENE = setKEGG)
listd[['setKEGG']] <- list(
  upgesetKEGG,
  dngesetKEGG)
# saveRDS(listd, file = 'output/HL_enrich.Rds')
View(listd)

m=0
nl <- c('bp','cc','mf','kegg')
for(i in names(listd)){
  m = m + 1
  fn1 <- paste0(nl[m],'up.pdf')
  fn2 <- paste0(nl[m],'dn.pdf')
  
  f1 <- dotplot(listd[[i]][[1]])
  f2 <- dotplot(listd[[i]][[2]])
  
  # ggsave(f1, filename = fn1, width = 13, height = 7)
  # ggsave(f2, filename = fn2, width = 13, height = 7)
}
names(listd)

# fn1 <- paste0('bp','dn.pdf')
# f1 <- dotplot(listd[['bp']][[2]])
# ggsave(f1, filename = fn1, width = 20, height = 7,limitsize = FALSE)
# 
# 
status <- train_y[,2]
time <- train_y[,1]
ROC.DSST <- timeROC(
  T=time,
  delta = status,
  marker = ROCR.simple,
  cause=1,
  weighting="marginal",
  times = c(365, 365 *3, 365*5, 365*10),
  iid=TRUE,
  ROC = T)

# plot(ROC.DSST,time=365)
# plot(ROC.DSST,time=365 *3,add=TRUE,col="orange")
# plot(ROC.DSST,time=365*5,add=TRUE,col="blue")
# 
# plot(ROC.DSST,time=365*7,add=TRUE,col="green") 
# plot(ROC.DSST,time=365*10,add=TRUE,col="grey50") 
# 
# legend("bottomright",c("Y-5","Y-3","Y-10"),
#        col=c("blue","orange","grey50"),
#        lty=1,lwd=2)
new_dat <- train_y
new_dat <- as.data.frame(new_dat)
new_dat$fp <- ROCR.simple
result <- with(new_dat, ROC.DSST)
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(365, 365 *3,365*5, 365*10)),
                            each = nrow(result$TP)))

p <- ggplot() + 
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", "#66C2A5","#BC3C29FF"),
                     labels = paste0("AUC of ",c(1,3,5, 10),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
p
# ggsave(p, filename = 'ROC.pdf',width = 7, height = 7)
pdf('cox_p.pdf')
ggsurvplot(survfit(res.cox), palette = "#2E9FDF",
           ggtheme = theme_minimal(), data =  lung)
dev.off()


ROC.DSST0 <- timeROC(
  T=time,
  delta = status,
  marker = ROCR.simple,
  cause=1,
  weighting="marginal",
  times = c(365, 365 *3, 365*5,365*7, 365*10),
  iid=TRUE,
  ROC = T)
saveRDS(ROC.DSST0, file = 'timeROC.Rds')

pd <- data.frame(
  time = c(365, 365 *3, 365*5,365*7, 365*10),
  survProb = ROC.DSST0$survProb,
  AUC = ROC.DSST0$AUC
)

pAUC <- ggplot(pd, aes(x = time, y = AUC)) + geom_point(size = 5) +
  ylim(0.5,1) + geom_line(size = 3) + theme_bw()
# ggsave(pAUC, filename = 'pAUC.pdf')
psvp <- ggplot(pd, aes(x = time, y = survProb)) + geom_point(size = 5) +
  ylim(0.1,1) + geom_line(size = 3) + theme_bw()                    
# ggsave(psvp, filename = 'pAUC.pdf')
# =============================== 韦恩图 ========================================
gndff <- rownames(top.table)[top.table$adj.P.Val < 0.05 & abs(top.table$logFC) > 0.7]
svgen <- rownames(sv_gld)[sv_gld$zph_p > 0.3 & sv_gld$p_sv < 0.05 & 
                            sv_gld$C_index > 0.58]
library(ggvenn)
a <- list(`DEGs` = gndff,
          `Prognostic Gene` = svgen
)
pdf('ven.pdf')
ggvenn(a, c("DEGs", "Prognostic Gene")) 
dev.off()
# --------------------------- 临床数据 --------------------------------------
library(tidyverse)
alldata <- readRDS('./all_sample.Rds')
patient_data <- alldata[["patient_data"]]
log2_TPM <- alldata$log2_TPM
v <- apply(log2_TPM, 1, mean)
log2_TPM <- log2_TPM[v > 1,]
noTR <- patient_data$type != 'TR' & !is.na(patient_data$type)
noTR_patient_data <- patient_data[noTR,]
noTR_log2_TPM <- log2_TPM[,noTR]
noTR_log2_TPM <- as.data.frame(noTR_log2_TPM)

data0 <- t(noTR_log2_TPM[,noTR_patient_data$type == 'TP'])
data0 <- as.data.frame(data0)
dim(data0)
noNT_noTR_patient_data <- noTR_patient_data[noTR_patient_data$type == 'TP',]
data0$time <- noNT_noTR_patient_data$days
data0$status <- 0
data0$status[noNT_noTR_patient_data$status == 'Dead'] <- 1
data0$barcode <- noNT_noTR_patient_data$patient
# data0$age <- noNT_noTR_patient_data$age_at_index
# data0$gender <- noNT_noTR_patient_data$gender
# data0$T_stage <- noNT_noTR_patient_data$ajcc_pathologic_t
# data0$N_stage <- noNT_noTR_patient_data$ajcc_pathologic_n
# data0$M_stage <- noNT_noTR_patient_data$ajcc_pathologic_m

data0 <- na.omit(data0)
data0 <- data0[data0[,'time'] > 0,]
noNT_noTR_patient_data <- noNT_noTR_patient_data[
  noNT_noTR_patient_data$patient %in% data0$barcode,]
data0$patient <- noNT_noTR_patient_data$patient
# saveRDS(data0, file = 'TCGA_log2TPM.Rds')
ROCR.simple <- readRDS('TCGA_HR.Rds')
tmp <- cbind(noNT_noTR_patient_data[,c(3,4,7,8,9,10,11,12)],ROCR.simple)
tmp$stage_age <- '>= 60'
tmp$stage_age[tmp$age_at_index < 60] <- '< 60'
tmp$stage_age[is.na(tmp$age_at_index)] <- NA

tmp$stage_T <- str_extract(tmp$ajcc_pathologic_t, pattern = '..') 
tmp$stage_M <- str_extract(tmp$ajcc_pathologic_m, pattern = '..') 

tmp$risk <- 'high'
tmp$risk[tmp$ROCR.simple < 0] <- 'low'

tmp$stage_T0 <- 'T2|3|4'
tmp$stage_T0[tmp$stage_T %in% c('T1')] <- 'T1'
tmp$stage_T0[tmp$stage_T %in% c('TX')] <- 'TX'

tmp$stage_N <- 'N0'
tmp$stage_N[tmp$ajcc_pathologic_n %in% c('N1', 'N2', 'N3')] <- 'N1|2|3'
tmp$stage_N[tmp$ajcc_pathologic_n %in% c('NX')] <- 'NX'
tmp$stage_N[is.na(tmp$ajcc_pathologic_n)] <- NA

tmp$pathologic_stage <- 'Stage1'
tmp$pathologic_stage[tmp$ajcc_pathologic_stage %in% c('Stage II',
                                                      'Stage IIA',
                                                      'Stage IIB',
                                                      'Stage IIIA',
                                                      'Stage IIIB',
                                                      'Stage IV')] <- 'Stage2|3|4'
tmp$pathologic_stage[is.na(tmp$ajcc_pathologic_stage)] <- NA

table(tmp$risk)
# 年龄
table(tmp$risk, tmp$stage_age)
chisq.test(table(tmp$risk, tmp$stage_age))
# 性别
table(tmp$risk, tmp$gender)
chisq.test(
  table(tmp$risk, tmp$gender)
)
# T
table(tmp$risk, tmp$stage_T0)
chisq.test(
  table(tmp$risk, tmp$stage_T0)
)
# M
table(tmp$risk, tmp$stage_M)
chisq.test(
  table(tmp$risk, tmp$stage_M)
)

# N
table(tmp$risk, tmp$stage_N)
chisq.test(
  table(tmp$risk, tmp$stage_N)
)

#Stage
table(tmp$risk, tmp$pathologic_stage)
chisq.test(table(tmp$risk, tmp$pathologic_stage))


tmp[,c(5,6,7)] <- apply(tmp[,c(5,6,7)], 2, function(x) { 
  x[grep(pattern = 'X',x)] <- NA
  x
}
)
tmp[,c(5,6,7)] <- apply(tmp[,c(5,6,7)], 2, function(x) { 
  unlist(str_extract_all(x, pattern = '[0-9]'))
}
)
tmpxxx <- tmp

tmp$main_stage <- NA
tmp$main_stage[tmp$ajcc_pathologic_stage %in% c('Stage I','Stage IA','Stage IB')] <- 'Stage1'
tmp$main_stage[tmp$ajcc_pathologic_stage %in% c('Stage II','Stage IIA','Stage IIB')] <- 'Stage2'
tmp$main_stage[tmp$ajcc_pathologic_stage %in% c('Stage IIIA','Stage IIIB')] <- 'Stage3/4'
tmp$main_stage[tmp$ajcc_pathologic_stage %in% c('Stage IV')] <- 'Stage3/4'

tmp$statu <- 0
tmp$statu[tmp$status == 'Dead'] <- 1

tmp0 <- tmp[,c(2:7,9,17,18)]
tmp0$ajcc_pathologic_t <- paste0('T',tmp0$ajcc_pathologic_t)
tmp0$ajcc_pathologic_n <- paste0('N',tmp0$ajcc_pathologic_n)
tmp0$ajcc_pathologic_t[tmp0$ajcc_pathologic_t %in% c('T3','T4')] <- 'T3/4'

tmp0$ajcc_pathologic_t <- factor(tmp0$ajcc_pathologic_t,
                                 levels = c('T1',
                                            'T2',
                                            'T3/4'))
tmp0$ajcc_pathologic_n <- factor(tmp0$ajcc_pathologic_n,
                                 levels = c('N0',
                                            'N1',
                                            'N2'))
table(tmp0$ajcc_pathologic_n)


colnames(tmp0)
library(survival)
clist <- list()
glc <- list()
for(i in 2:8){
  coxt <- tmp0[,c(1,9,i)]
  res.cox <- coxph(Surv(days, statu) ~ ., data =  coxt)
  cox.zph(res.cox)$table
  res.cox <- summary(res.cox)
  glc[[colnames(tmp0)[i]]] <-  as.data.frame(unlist(res.cox$coefficients))
  clist[[colnames(tmp0)[i]]] <- data.frame(C = res.cox$concordance[1],
                                           se = res.cox$concordance[2])
}
glc <- purrr::reduce(glc, rbind)
clist <- purrr::reduce(clist,rbind)

View(glc)
min(clist[,2])
clist <- as.data.frame(clist)
clist$names <- colnames(tmp0)[2:8]
colnames(clist)
clist <- clist[order(clist$C),]
clist$names <- factor(clist$names,levels = c(clist$names))
ggplot(clist,aes(x = names, y = C, 
                 color = names, group = names)) + 
  geom_segment(aes(x = names, y = 0, 
                   yend = C, xend = names),
               color = 'grey',size = 1) + 
  geom_point(size = 4)+
  ylim(0,0.75)+
  coord_flip() + 
  theme_bw() -> TCGACindex
write.csv(clist, file = 'TCGA_Cindex.csv')
ggsave(TCGACindex, width = 5.5,height = 3.4, 
       filename = 'TCGACindex.pdf')

# library('compareC')
# tmp0 <- na.omit(tmp0)
# tcgaclist <- list()
# for(i in c(2,4,5,6,7,8))
# {
#   compareC(tmp0$days, tmp0$statu, 
#            tmp0[,i], tmp0$ROCR.simple) -> ccc
#   tcgaclist[[i]] <- ccc$pval
# }
# cox.zph(res.cox)$table
# summary(res.cox)
# rms::vif(res.cox)
# mld <- names(rms::vif(res.cox))
glc$Upper <- glc$coef + 1.96 * glc$`se(coef)`
glc$Lower <- glc$coef - 1.96 * glc$`se(coef)`
glc$Varnames <- rownames(glc)
View(glc)
p <- ggplot(glc, aes(x = coef, y = Varnames, 
                     group = Varnames, col = `Pr(>|z|)`)) # 不同形状shape= Factor

forsv <- p + geom_point(size=3) +
  
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.3,size=2) +
  
  scale_x_continuous(limits= c(-0.5, 2), breaks= seq(-2, 2, 0.2)) +
  
  geom_vline(aes(xintercept = 0)) +
  
  xlab('Coef') + ylab(' ') + theme_bw() + scale_color_viridis_c()
forsv
ggsave(forsv, filename = 'single_forest.pdf',width = 7,
       height = 4)
write.csv(glc, file = 'single_TCGA.csv')

tmp00 <- tmp0[,c(1,4,5,6,7,8,9)]
res.cox <- coxph(Surv(days, statu) ~ ., data =  tmp00)
cox.zph(res.cox)$table
res.cox <- summary(res.cox)
res.cox <- res.cox$coefficients
res.cox <- as.data.frame(res.cox)
res.cox$Lower <- res.cox$coef - 1.96 * res.cox$`se(coef)`
res.cox$Upper <- res.cox$coef + 1.96 * res.cox$`se(coef)`
res.cox$Varnames <- rownames(res.cox)


p <- ggplot(res.cox, aes(x = coef, y = Varnames, 
                         group = Varnames, col = `Pr(>|z|)`)) # 不同形状shape= Factor

forsv <- p + geom_point(size=3) +
  
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.3,size=2) +
  
  scale_x_continuous(limits= c(-1, 2), 
                     breaks= seq(-1, 2, 0.5)) +
  
  geom_vline(aes(xintercept = 0)) +
  
  xlab('Coef') + ylab(' ') + theme_bw() + scale_color_viridis_c()
forsv
ggsave(forsv, filename = 'mulit_forest.pdf',
       width = 7,height = 4)
write.csv(res.cox, file = 'multi_TCGA.csv')

tmpxxx <- tmpxxx[,c(5,6,7,8,9)]
tmpxxx <- na.omit(tmpxxx)
tmpxxx$ajcc_pathologic_stage <- gsub(tmpxxx$ajcc_pathologic_stage, pattern = '(A|B)', replacement = "")
tmpxxx$ajcc_pathologic_t <- paste0('T', tmpxxx$ajcc_pathologic_t)
tmpxxx$ajcc_pathologic_n <- paste0('N', tmpxxx$ajcc_pathologic_n)
tmpxxx$ajcc_pathologic_m <- paste0('M', tmpxxx$ajcc_pathologic_m)
library(ggpubr)
mycomparision <- list(c('T1','T2'),
                      c('T1','T3'),
                      c('T1','T4'))
Tp <- ggplot(tmpxxx, aes(x = ajcc_pathologic_t, 
                         y = ROCR.simple, fill = ajcc_pathologic_t))+ 
  geom_boxplot() + theme_bw() + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, label = 'p.signif')
Tp
ggsave(Tp, filename = 'T_TCGA.pdf',width = 7,height = 7)
mycomparision <- list(c('N0','N1'),
                      c('N0','N2'))
Np <- ggplot(tmpxxx, aes(x = ajcc_pathologic_n, 
                         y = ROCR.simple, fill = ajcc_pathologic_n))+ 
  geom_boxplot() + theme_bw() + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, label = 'p.signif')
Np
ggsave(Np, filename = 'N_TCGA.pdf',width = 7,height = 7)
mycomparision <- list(c('M0','M1'))
Mp <- ggplot(tmpxxx, aes(x = ajcc_pathologic_m, 
                         y = ROCR.simple, fill = ajcc_pathologic_m))+ 
  geom_boxplot() + theme_bw() + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, label = 'p.signif')
Mp
ggsave(Mp, filename = 'M_TCGA.pdf',width = 3.5,height = 7)
mycomparision <- list(c('Stage I', 'Stage II'),
                      c('Stage I', 'Stage III'),
                      c('Stage I', 'Stage IV'))
Sp <- ggplot(tmpxxx, aes(x = ajcc_pathologic_stage, 
                         y = ROCR.simple, fill = ajcc_pathologic_stage))+ 
  geom_boxplot() + theme_bw() + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, label = 'p.signif')
Sp
ggsave(Sp, filename = 'S_TCGA.pdf',width = 7,height = 7)

# ===================== PCA ========================
ROCR.simple <- readRDS('TCGA_HR.Rds')
lung <- readRDS('TCGA_modle_lung.Rds')

com1 <- prcomp(lung, center = TRUE,scale. = TRUE)
pc <- com1$x[,1:2]
pc <- as.data.frame(pc)
pc$group <- ROCR.simple

pc$level <- 'high'
pc$level[pc$group < 0] <- 'low'

pdf('PCA_TCGA.pdf')
ggplot(pc, aes(x = PC1, y = PC2, color = level)) + 
  geom_point() + scale_color_manual(
    values = c(
      'high' = 'red',
      'low' = 'blue')
  ) + theme_bw()
dev.off()
# ====================== 大盒子 ========================
rm(list = ls())
res.cox <- readRDS('res.cox_model.Rds')
alldata <- readRDS('./all_sample.Rds')

patient_data <- alldata[["patient_data"]]
log2_TPM <- alldata$log2_TPM

v <- apply(log2_TPM, 1, mean)
log2_TPM <- log2_TPM[v > 1,]

noTR <- patient_data$type != 'TR' & 
  (!is.na(patient_data$type)) & 
  (!is.na(alldata$patient_data$days)) &
  alldata$patient_data$days != 0

noTR_patient_data <- patient_data[noTR,]
noTR_log2_TPM <- log2_TPM[,noTR]
dim(noTR_log2_TPM)
length(noTR_patient_data$type)

d <- as.data.frame(t(noTR_log2_TPM[names(res.cox$coefficients),]))
d$x <- factor(noTR_patient_data$type, levels = c('TP','NT'))

library(reshape2)
d <- melt(d)
GSE10072 <- ggplot(d, aes(x = variable, y = value, color = x)) + 
  geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = 
          element_text(angle = 45, 
                       hjust = 0.5, 
                       vjust = 0.5)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
pdf('TCGA_comfirmBox.pdf',width = 7,height = 7)
GSE10072
dev.off()




# TPM <- read.csv('data/TPM.csv')
# res.cox_model <- readRDS('res.cox_model.Rds')
# 
# tmp <- TPM[TPM$gene_name %in% names(res.cox_model$coefficients),c(5,6,7,8,9,10)]
# rownames(tmp) <- TPM[TPM$gene_name %in% names(res.cox_model$coefficients),2]
# write.csv(tmp, file = 'RNAseq_TPM_modle.csv')
# 
# library(ComplexHeatmap)
# tmp <- log2(1+tmp)
# tmp <- t(apply(tmp, 1, scale))
# colnames(tmp) <- c('C1','N1','NP1','NP2','P1','P2')
# tmp[is.na(tmp)] <- 0
# pdf('RNAseq_TPM_modle.pdf')
# Heatmap(tmp)
# dev.off()

# =========================== 富集分析基因集合准备 =============================
library(GSEABase)
library(GSVA)
library(clusterProfiler)

res.cox_model <- readRDS('res.cox_model.Rds')
ERstress <- getGmt('data/genesets.gmt')
imm <- read.csv('data/mmc3.csv')

driver <- read.csv('D:/GeneSet/Fe/driver.csv')
marker <- read.csv('D:/GeneSet/Fe/marker.csv')
suppressor <- read.csv('D:/GeneSet/Fe/suppressor.csv')

list(
  driver = unique(driver$symbol),
  suppressor = unique(suppressor$symbol)
) -> fe


saveRDS(fe, file = 'fedb.Rds')
saveRDS(ERstress, file = 'ERstressdb.Rds')
saveRDS(imm, file = 'immdb.Rds')
# ==============================================================================
rm(list = ls())
library(GSEABase)
library(GSVA)
library(ggplot2)
library(reshape2)
library(ComplexHeatmap)
library(circlize)
library(clusterProfiler)
library(ggpmisc)

# fe <- readRDS('fedb.Rds')
imm <- readRDS('immdb.Rds')
ERstress <- readRDS('ERstressdb.Rds')
ERstressL <- read.gmt('data/genesets.gmt')
top100_protein <- read.table('ExoCarta_top100_protein_details_5.txt',
                             quote = '',sep = '\t',header = T)
top100_protein <- top100_protein[,1]
top100_protein <- gsub(top100_protein, pattern = ' ',replacement = '')
# tensor <- read.csv('张力annurev.csv')
# tensor <- unique(tensor[,1])

imml <- list()
for(i in unique(imm$Cell.type)){
  imml[[i]] <- imm[imm$Cell.type == i,1]
}

# ======================== 开始画图GEO1 ==========================
GEO1 <- readRDS('GEO1_log2.Rds')
HR_GEO1 <- readRDS('HR_se1p.Rds')
dir.create('GEO1next')
data0 <- as.matrix(GEO1)
BPPARAM = BiocParallel::SnowParam(4)
top100_proteinL <- list(top100_protein = top100_protein)
gset_top100_protein <- gsva(data0, top100_proteinL, method='ssgsea', 
                            parallel.sz = 0, min.sz = 5,
                            BPPARAM = BPPARAM)
# tensorL <- list(tensor = tensor)
# gset_tensor <- gsva(data0, tensorL, method='ssgsea',
#                     parallel.sz = 0, min.sz = 5,
#                     BPPARAM = BPPARAM)

# gset_fe <- gsva(data0, fe, method='ssgsea',
#                  parallel.sz = 0, min.sz = 5,
#                  BPPARAM = BPPARAM)
gset_imm <- gsva(data0, imml, method='ssgsea',
                 parallel.sz = 0, min.sz = 5,
                 BPPARAM = BPPARAM)
gset_ERstress <- gsva(data0, ERstress, method='ssgsea',
                      parallel.sz = 0, min.sz = 5,
                      BPPARAM = BPPARAM)

HR_GEO1_level <- ifelse(HR_GEO1 > 0, 'High', 'Low')
gset_top100_protein_data <- data.frame(t(gset_top100_protein))
gset_top100_protein_data$HR_GEO1_level <- HR_GEO1_level
gset_top100_protein_dataL <- melt(gset_top100_protein_data)

# gset_tensor_data <- data.frame(t(gset_tensor))
# gset_tensor_data$HR_GEO1_level <- HR_GEO1_level
# gset_tensor_dataL <- melt(gset_tensor_data)

# gset_fe_data <- data.frame(t(gset_fe))
# gset_fe_data$HR_GEO1_level <- HR_GEO1_level
# gset_fe_dataL <- melt(gset_fe_data)

gset_imm_data <- data.frame(t(gset_imm))
gset_imm_data$HR_GEO1_level <- HR_GEO1_level
gset_imm_dataL <- melt(gset_imm_data)

gset_ERstress_data <- data.frame(t(gset_ERstress))
gset_ERstress_data$HR_GEO1_level <- HR_GEO1_level
gset_ERstress_dataL <- melt(gset_ERstress_data)



library(ggpubr)
# ----------------------- 盒子图 ----------------------
gset_top100_protein_dataL$HR_GEO1_level <- factor(gset_top100_protein_dataL$HR_GEO1_level,
                                                  levels = unique(gset_top100_protein_dataL$HR_GEO1_level))
top100box <- ggplot(gset_top100_protein_dataL, aes(x = variable, y = value,
                                                   fill = HR_GEO1_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(top100box, filename = 'GEO1next/top100box.pdf', height = 7,width = 8)

gset_tensor_dataL$HR_GEO1_level <- factor(gset_tensor_dataL$HR_GEO1_level,
                                          levels = unique(gset_tensor_dataL$HR_GEO1_level))
tensorbox <- ggplot(gset_tensor_dataL, aes(x = variable, y = value,
                                           fill = HR_GEO1_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(tensorbox, filename = 'GEO1next/tensorbox.pdf', height = 7,width = 8)

gset_fe_dataL$HR_GEO1_level <- factor(gset_fe_dataL$HR_GEO1_level,
                                      levels = unique(gset_fe_dataL$HR_GEO1_level))
febox <- ggplot(gset_fe_dataL, aes(x = variable, y = value,
                                   fill = HR_GEO1_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(febox, filename = 'GEO1next/feubox.pdf', height = 7,width = 8)


gset_imm_dataL$HR_GEO1_level <- factor(gset_imm_dataL$HR_GEO1_level,
                                       levels = unique(gset_imm_dataL$HR_GEO1_level))
immbox <- ggplot(gset_imm_dataL, aes(x = variable, y = value,
                                     fill = HR_GEO1_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(immbox, filename = 'GEO1next/immubox.pdf', height = 7,width = 8)

gset_ERstress_dataL$HR_GEO1_level <- factor(gset_ERstress_dataL$HR_GEO1_level,
                                            levels = unique(gset_ERstress_dataL$HR_GEO1_level))
ERbox <- ggplot(gset_ERstress_dataL, aes(x = variable, y = value,
                                         fill = HR_GEO1_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(ERbox, filename = 'GEO1next/ERbox.pdf', height = 20,width = 10)


# ----------------------- 打分热图 + 基因热图 -----------------------
od <- order(HR_GEO1)
HR_GEO1_level_order <- HR_GEO1_level[od]

mxgene_order <- GEO1[rownames(GEO1) %in% unique(top100_protein),od]
an <- HeatmapAnnotation(group = HR_GEO1_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
mxgene_order_scale <- t(apply(mxgene_order, 1, scale))
pdf('GEO1next/GEO1_heatmap_top100_gene.pdf', width = 7,height = 15)
Heatmap(mxgene_order_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_order <- GEO1[rownames(GEO1) %in% unique(tensor),od]
an <- HeatmapAnnotation(group = HR_GEO1_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
mxgene_order_scale <- t(apply(mxgene_order, 1, scale))
pdf('GEO1next/GEO1_heatmap_tensor_gene.pdf', width = 7,height = 8)
Heatmap(mxgene_order_scale, cluster_columns = F, top_annotation = an)
dev.off()

gset_fe_order <- gset_fe[,od]
mxgene_order <- GEO1[rownames(GEO1) %in% unique(c(fe$driver, fe$suppressor)),od]
an <- HeatmapAnnotation(group = HR_GEO1_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
gset_imm_scale <- t(apply(gset_fe_order, 1, scale))
pdf('GEO1next/GEO1_heatmap_fe_cell.pdf', width = 7,height = 4)
Heatmap(gset_imm_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_scale <- t(apply(mxgene_order, 1, scale))
pdf('GEO1next/GEO1_heatmap_fe_gene.pdf', width = 7,height = 100)
Heatmap(mxgene_scale, cluster_columns = F, top_annotation = an)
dev.off()

gset_imm_order <- gset_imm[,od]
mxgene_order <- GEO1[rownames(GEO1) %in% unique(imm$Metagene),od]
an <- HeatmapAnnotation(group = HR_GEO1_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
gset_imm_scale <- t(apply(gset_imm_order, 1, scale))
pdf('GEO1next/GEO1_heatmap_immu_cell.pdf', width = 7,height = 4)
Heatmap(gset_imm_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_scale <- t(apply(mxgene_order, 1, scale))
pdf('GEO1next/GEO1_heatmap_immu_gene.pdf', width = 7,height = 100)
Heatmap(mxgene_scale, cluster_columns = F, top_annotation = an)
dev.off()

ER_order <- GEO1[rownames(GEO1) %in% unique(ERstressL$gene),od]
an <- HeatmapAnnotation(group = HR_GEO1_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
ER_scale <- t(apply(ER_order, 1, scale))
pdf('GEO1next/GEO1_heatmap_ER_gene.pdf', width = 7,height = 150)
Heatmap(ER_scale, cluster_columns = F, top_annotation = an)
dev.off()


# -------------------------- 通路相关性图 -----------------------
gset_top100_protein_data <- data.frame(t(gset_top100_protein))
gset_top100_protein_cdata <- apply(gset_top100_protein_data, 2, function(x) cor(x, HR_GEO1))
Cor_top100 <- 
  data.frame(cor = gset_top100_protein_cdata, 
             absCor = abs(unlist(gset_top100_protein_cdata)),
             name = colnames(gset_top100_protein_data))
gset_top100_protein_data$HR <- HR_GEO1

dir.create('GEO1next/ssgseafe')
for(i in 1:1){
  ns <- colnames(gset_top100_protein_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_top100_protein_data, 
                aes(x = eval(parse(text = colnames(gset_top100_protein_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) +
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('GEO1next/ssgseafe/',ns,'.pdf'), 
         height = 7, width = 7)
}

# gset_tensor_data <- data.frame(t(gset_tensor))
# gset_tensor_cdata <- apply(gset_tensor_data, 2, function(x) cor(x, HR_GEO1))
# Cor_tensor <- 
#   data.frame(cor = gset_tensor_cdata, 
#              absCor = abs(unlist(gset_tensor_cdata)),
#              name = colnames(gset_tensor_data))
# gset_tensor_data$HR <- HR_GEO1
# 
# for(i in 1:1){
#   ns <- colnames(gset_tensor_data)[i]
#   my.formula <- y ~ x
#   fig <- ggplot(gset_tensor_data, 
#                 aes(x = eval(parse(text = colnames(gset_tensor_data)[i])), y = HR)) + 
#     geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
#     stat_poly_eq(formula = my.formula, 
#                  aes(label = paste(..rr.label.., sep = "~~~")), 
#                  parse = TRUE) +
#     geom_point() + theme_bw() + xlab(ns)
#   ggsave(fig, filename = paste0('GEO1next/ssgseafe/',ns,'.pdf'), 
#          height = 7, width = 7)
# }


# gset_fe_data <- data.frame(t(gset_fe))
# data_gset_fe_cdata <- apply(gset_fe_data, 2, function(x) cor(x, HR_GEO1))
# Cor_fe <- 
#   data.frame(cor = data_gset_fe_cdata, 
#              absCor = abs(unlist(data_gset_fe_cdata)),
#              name = colnames(gset_fe_data))
# gset_fe_data$HR <- HR_GEO1
# 
# dir.create('GEO1next/ssgseafe')
# for(i in 1:2){
#   ns <- colnames(gset_fe_data)[i]
#   my.formula <- y ~ x
#   fig <- ggplot(gset_fe_data, 
#                 aes(x = eval(parse(text = colnames(gset_fe_data)[i])), y = HR)) + 
#     geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
#     stat_poly_eq(formula = my.formula, 
#                  aes(label = paste(..rr.label.., sep = "~~~")), 
#                  parse = TRUE) +
#     geom_point() + theme_bw() + xlab(ns)
#   ggsave(fig, filename = paste0('GEO1next/ssgseafe/',ns,'.pdf'), 
#          height = 7, width = 7)
# }



gset_imm_data <- data.frame(t(gset_imm))
data_gset_imm_cdata <- apply(gset_imm_data, 2, function(x) cor(x, HR_GEO1))
Cor_imm <- 
  data.frame(cor = data_gset_imm_cdata, 
             absCor = abs(unlist(data_gset_imm_cdata)),
             name = colnames(gset_imm_data))
gset_imm_data$HR <- HR_GEO1

dir.create('GEO1next/ssgseaimmu')
for(i in 1:28){
  ns <- colnames(gset_imm_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_imm_data, 
                aes(x = eval(parse(text = colnames(gset_imm_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) +
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('GEO1next/ssgseaimmu/',ns,'.pdf'), 
         height = 7, width = 7)
}




gset_ERstress_data <- data.frame(t(gset_ERstress))
data_gset_ERstress_cdata <- apply(gset_ERstress_data, 2, function(x) cor(x, HR_GEO1))
Cor_ERstress <- 
  data.frame(cor = data_gset_ERstress_cdata, 
             absCor = abs(unlist(data_gset_ERstress_cdata)),
             name = colnames(gset_ERstress_data))
gset_ERstress_data$HR <- HR_GEO1
dir.create('GEO1next/ssgseaERstress')
for(i in 1:23){
  ns <- colnames(gset_ERstress_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_ERstress_data, 
                aes(x = eval(parse(text = colnames(gset_ERstress_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) + 
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('GEO1next/ssgseaERstress/',ns,'.pdf'), 
         height = 7, width = 7)
}


# Cor_fe <- Cor_fe[order(Cor_fe$cor),]
# Cor_fe$name <- factor(Cor_fe$name, levels = Cor_fe$name)
# 
# Cor_feC <- ggplot(Cor_fe, aes(x=name, y=cor, color = absCor)) +
#   geom_point(aes(size = absCor)) + 
#   geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
#   scale_color_gradient2(mid = 'green', high = 'red')+
#   theme_bw() + coord_flip()
# ggsave(Cor_feC, filename ='GEO1next/Cor_feC.pdf', 
#        height = 7, width = 7)

Cor_imm <- Cor_imm[order(Cor_imm$cor),]
Cor_imm$name <- factor(Cor_imm$name, levels = Cor_imm$name)

Cor_immC <- ggplot(Cor_imm, aes(x=name, y=cor, color = absCor)) +
  geom_point(aes(size = absCor)) + 
  geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
  scale_color_gradient2(mid = 'green', high = 'red')+
  theme_bw() + coord_flip()
ggsave(Cor_immC, filename ='GEO1next/Cor_immC.pdf', 
       height = 7, width = 7)

Cor_ERstress <- Cor_ERstress[order(Cor_ERstress$cor),]
Cor_ERstress$name <- factor(Cor_ERstress$name, levels = Cor_ERstress$name)
Cor_ERstress2 <- ggplot(Cor_ERstress, aes(x=name, y=cor, color = absCor)) +
  geom_point(aes(size = absCor)) + 
  geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
  scale_color_gradient2(mid = 'green', high = 'red')+
  theme_bw() + coord_flip()
ggsave(Cor_ERstress2, filename = 'GEO1next/Cor_ERstress2.pdf',
       height = 7, width = 15)
# ======================== 开始画图GEO2 ==========================
GEO2 <- readRDS('GEO2_log2.Rds')
HR_GEO2 <- readRDS('HR_GEO2.Rds')
dir.create('GEO2next')
data0 <- as.matrix(GEO2)
BPPARAM = BiocParallel::SnowParam(4)
top100_proteinL <- list(top100_protein = top100_protein)
gset_top100_protein <- gsva(data0, top100_proteinL, method='ssgsea', 
                            parallel.sz = 0, min.sz = 5,
                            BPPARAM = BPPARAM)
# tensorL <- list(tensor = tensor)
# gset_tensor <- gsva(data0, tensorL, method='ssgsea',
#                     parallel.sz = 0, min.sz = 5,
#                     BPPARAM = BPPARAM)
# gset_fe <- gsva(data0, fe, method='ssgsea',
#                  parallel.sz = 0, min.sz = 5,
#                  BPPARAM = BPPARAM)
gset_imm <- gsva(data0, imml, method='ssgsea',
                 parallel.sz = 0, min.sz = 5,
                 BPPARAM = BPPARAM)
gset_ERstress <- gsva(data0, ERstress, method='ssgsea',
                      parallel.sz = 0, min.sz = 5,
                      BPPARAM = BPPARAM)

HR_GEO2_level <- ifelse(HR_GEO2 > 0, 'High', 'Low')
gset_top100_protein_data <- data.frame(t(gset_top100_protein))
gset_top100_protein_data$HR_GEO2_level <- HR_GEO2_level
gset_top100_protein_dataL <- melt(gset_top100_protein_data)

# gset_tensor_data <- data.frame(t(gset_tensor))
# gset_tensor_data$HR_GEO2_level <- HR_GEO2_level
# gset_tensor_dataL <- melt(gset_tensor_data)

# gset_fe_data <- data.frame(t(gset_fe))
# gset_fe_data$HR_GEO2_level <- HR_GEO2_level
# gset_fe_dataL <- melt(gset_fe_data)


gset_imm_data <- data.frame(t(gset_imm))
gset_imm_data$HR_GEO2_level <- HR_GEO2_level
gset_imm_dataL <- melt(gset_imm_data)

gset_ERstress_data <- data.frame(t(gset_ERstress))
gset_ERstress_data$HR_GEO2_level <- HR_GEO2_level
gset_ERstress_dataL <- melt(gset_ERstress_data)


# ----------------------- 盒子图 ----------------------
gset_top100_protein_dataL$HR_GEO2_level <- factor(gset_top100_protein_dataL$HR_GEO2_level,
                                                  levels = unique(gset_top100_protein_dataL$HR_GEO2_level))
top100box <- ggplot(gset_top100_protein_dataL, aes(x = variable, y = value,
                                                   fill = HR_GEO2_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(top100box, filename = 'GEO2next/top100box.pdf', height = 7,width = 8)

gset_tensor_dataL$HR_GEO2_level <- factor(gset_tensor_dataL$HR_GEO2_level,
                                          levels = unique(gset_tensor_dataL$HR_GEO2_level))
tensorbox <- ggplot(gset_tensor_dataL, aes(x = variable, y = value,
                                           fill = HR_GEO2_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(tensorbox, filename = 'GEO2next/tensorbox.pdf', height = 7,width = 8)

gset_fe_dataL$HR_GEO2_level <- factor(gset_fe_dataL$HR_GEO2_level,
                                      levels = unique(gset_fe_dataL$HR_GEO2_level))
febox <- ggplot(gset_fe_dataL, aes(x = variable, y = value,
                                   fill = HR_GEO2_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(febox, filename = 'GEO2next/immubox.pdf', height = 7,width = 8)

gset_imm_dataL$HR_GEO2_level <- factor(gset_imm_dataL$HR_GEO2_level,
                                       levels = unique(gset_imm_dataL$HR_GEO2_level))
immbox <- ggplot(gset_imm_dataL, aes(x = variable, y = value,
                                     fill = HR_GEO2_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(immbox, filename = 'GEO2next/immubox.pdf', height = 7,width = 8)

gset_ERstress_dataL$HR_GEO2_level <- factor(gset_ERstress_dataL$HR_GEO2_level,
                                            levels = unique(gset_ERstress_dataL$HR_GEO2_level))
ERbox <- ggplot(gset_ERstress_dataL, aes(x = variable, y = value,
                                         fill = HR_GEO2_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(ERbox, filename = 'GEO2next/ERbox.pdf', height = 20,width = 10)

# ----------------------- 打分热图 + 基因热图 -----------------------
od <- order(HR_GEO2)
HR_GEO2_level_order <- HR_GEO2_level[od]

mxgene_order <- GEO2[rownames(GEO2) %in% unique(top100_protein),od]
an <- HeatmapAnnotation(group = HR_GEO2_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
mxgene_order_scale <- t(apply(mxgene_order, 1, scale))
pdf('GEO2next/GEO2_heatmap_top100_gene.pdf', width = 7,height = 15)
Heatmap(mxgene_order_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_order <- GEO2[rownames(GEO2) %in% unique(tensor),od]
an <- HeatmapAnnotation(group = HR_GEO2_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
mxgene_order_scale <- t(apply(mxgene_order, 1, scale))
pdf('GEO2next/GEO2_heatmap_tensor_gene.pdf', width = 7,height = 8)
Heatmap(mxgene_order_scale, cluster_columns = F, top_annotation = an)
dev.off()

gset_fe_order <- gset_fe[,od]
mxgene_order <- GEO2[rownames(GEO2) %in% unique(c(fe$driver,fe$suppressor)),od]
an <- HeatmapAnnotation(group = HR_GEO2_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
gset_fe_scale <- t(apply(gset_fe_order, 1, scale))
pdf('GEO2next/GEO2_heatmap_fe_cell.pdf', width = 7,height = 4)
Heatmap(gset_fe_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_scale <- t(apply(mxgene_order, 1, scale))
pdf('GEO2next/GEO2_heatmap_fe_gene.pdf', width = 7,height = 100)
Heatmap(mxgene_scale, cluster_columns = F, top_annotation = an)
dev.off()


gset_imm_order <- gset_imm[,od]
mxgene_order <- GEO2[rownames(GEO2) %in% unique(imm$Metagene),od]
an <- HeatmapAnnotation(group = HR_GEO2_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
gset_imm_scale <- t(apply(gset_imm_order, 1, scale))
pdf('GEO2next/GEO2_heatmap_immu_cell.pdf', width = 7,height = 4)
Heatmap(gset_imm_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_scale <- t(apply(mxgene_order, 1, scale))
pdf('GEO2next/GEO2_heatmap_immu_gene.pdf', width = 7,height = 100)
Heatmap(mxgene_scale, cluster_columns = F, top_annotation = an)
dev.off()

ER_order <- GEO2[rownames(GEO2) %in% unique(ERstressL$gene),od]
an <- HeatmapAnnotation(group = HR_GEO2_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
ER_scale <- t(apply(ER_order, 1, scale))
pdf('GEO2next/GEO2_heatmap_ER_gene.pdf', width = 7,height = 150)
Heatmap(ER_scale, cluster_columns = F, top_annotation = an)
dev.off()


# -------------------------- 通路相关性图 -----------------------
gset_top100_protein_data <- data.frame(t(gset_top100_protein))
gset_top100_protein_cdata <- apply(gset_top100_protein_data, 2, function(x) cor(x, HR_GEO2))
Cor_top100 <- 
  data.frame(cor = gset_top100_protein_cdata, 
             absCor = abs(unlist(gset_top100_protein_cdata)),
             name = colnames(gset_top100_protein_data))
gset_top100_protein_data$HR <- HR_GEO2

dir.create('GEO2next/ssgseafe')
for(i in 1:1){
  ns <- colnames(gset_top100_protein_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_top100_protein_data, 
                aes(x = eval(parse(text = colnames(gset_top100_protein_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) +
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('GEO2next/ssgseafe/',ns,'.pdf'), 
         height = 7, width = 7)
}

# gset_tensor_data <- data.frame(t(gset_tensor))
# gset_tensor_cdata <- apply(gset_tensor_data, 2, function(x) cor(x, HR_GEO2))
# Cor_tensor <- 
#   data.frame(cor = gset_tensor_cdata, 
#              absCor = abs(unlist(gset_tensor_cdata)),
#              name = colnames(gset_tensor_data))
# gset_tensor_data$HR <- HR_GEO2
# 
# for(i in 1:1){
#   ns <- colnames(gset_tensor_data)[i]
#   my.formula <- y ~ x
#   fig <- ggplot(gset_tensor_data, 
#                 aes(x = eval(parse(text = colnames(gset_tensor_data)[i])), y = HR)) + 
#     geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
#     stat_poly_eq(formula = my.formula, 
#                  aes(label = paste(..rr.label.., sep = "~~~")), 
#                  parse = TRUE) +
#     geom_point() + theme_bw() + xlab(ns)
#   ggsave(fig, filename = paste0('GEO2next/ssgseafe/',ns,'.pdf'), 
#          height = 7, width = 7)
# }

# gset_fe_data <- data.frame(t(gset_fe))
# data_gset_fe_cdata <- apply(gset_fe_data, 2, function(x) cor(x, HR_GEO2))
# Cor_fe <- 
#   data.frame(cor = data_gset_fe_cdata, 
#              absCor = abs(unlist(data_gset_fe_cdata)),
#              name = colnames(gset_fe_data))
# gset_fe_data$HR <- HR_GEO2
# dir.create('GEO2next/ssgseafe')
# for(i in 1:2){
#   ns <- colnames(gset_fe_data)[i]
#   my.formula <- y ~ x
#   fig <- ggplot(gset_fe_data, 
#                 aes(x = eval(parse(text = colnames(gset_fe_data)[i])), y = HR)) + 
#     geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
#     stat_poly_eq(formula = my.formula, 
#                  aes(label = paste(..rr.label.., sep = "~~~")), 
#                  parse = TRUE) +
#     geom_point() + theme_bw() + xlab(ns)
#   ggsave(fig, filename = paste0('GEO2next/ssgseafe/',ns,'.pdf'), 
#          height = 7, width = 7)
# }

gset_imm_data <- data.frame(t(gset_imm))
data_gset_imm_cdata <- apply(gset_imm_data, 2, function(x) cor(x, HR_GEO2))
Cor_imm <- 
  data.frame(cor = data_gset_imm_cdata, 
             absCor = abs(unlist(data_gset_imm_cdata)),
             name = colnames(gset_imm_data))
gset_imm_data$HR <- HR_GEO2
dir.create('GEO2next/ssgseaimmu')
for(i in 1:28){
  ns <- colnames(gset_imm_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_imm_data, 
                aes(x = eval(parse(text = colnames(gset_imm_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) +
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('GEO2next/ssgseaimmu/',ns,'.pdf'), 
         height = 7, width = 7)
}

gset_ERstress_data <- data.frame(t(gset_ERstress))
data_gset_ERstress_cdata <- apply(gset_ERstress_data, 2, function(x) cor(x, HR_GEO2))
Cor_ERstress <- 
  data.frame(cor = data_gset_ERstress_cdata, 
             absCor = abs(unlist(data_gset_ERstress_cdata)),
             name = colnames(gset_ERstress_data))
gset_ERstress_data$HR <- HR_GEO2
dir.create('GEO2next/ssgseaERstress')
for(i in 1:23){
  ns <- colnames(gset_ERstress_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_ERstress_data, 
                aes(x = eval(parse(text = colnames(gset_ERstress_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) + 
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('GEO2next/ssgseaERstress/',ns,'.pdf'), 
         height = 7, width = 7)
}


# Cor_fe <- Cor_fe[order(Cor_fe$cor),]
# Cor_fe$name <- factor(Cor_fe$name, levels = Cor_fe$name)
# Cor_feC <- ggplot(Cor_fe, aes(x=name, y=cor, color = absCor)) +
#   geom_point(aes(size = absCor)) + 
#   geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
#   scale_color_gradient2(mid = 'green', high = 'red')+
#   theme_bw() + coord_flip()
# ggsave(Cor_feC, filename ='GEO2next/Cor_feC.pdf', 
#        height = 7, width = 7)

Cor_imm <- Cor_imm[order(Cor_imm$cor),]
Cor_imm$name <- factor(Cor_imm$name, levels = Cor_imm$name)

Cor_immC <- ggplot(Cor_imm, aes(x=name, y=cor, color = absCor)) +
  geom_point(aes(size = absCor)) + 
  geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
  scale_color_gradient2(mid = 'green', high = 'red')+
  theme_bw() + coord_flip()
ggsave(Cor_immC, filename ='GEO2next/Cor_immC.pdf', 
       height = 7, width = 7)

Cor_ERstress <- Cor_ERstress[order(Cor_ERstress$cor),]
Cor_ERstress$name <- factor(Cor_ERstress$name, levels = Cor_ERstress$name)
Cor_ERstress2 <- ggplot(Cor_ERstress, aes(x=name, y=cor, color = absCor)) +
  geom_point(aes(size = absCor)) + 
  geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
  scale_color_gradient2(mid = 'green', high = 'red')+
  theme_bw() + coord_flip()
ggsave(Cor_ERstress2, filename = 'GEO2next/Cor_ERstress2.pdf',
       height = 7, width = 15)
# ======================== 开始画图TCGA ==========================
TCGA <- readRDS('TCGA_log2TPM.Rds')
TCGA_HR <- readRDS('TCGA_HR.Rds')
dir.create('TCGAnext')
TCGA <- apply(TCGA, 2, function(x) as.numeric(x))
TCGA <- na.omit(t(TCGA))
data0 <- TCGA
BPPARAM = BiocParallel::SnowParam(4)
top100_proteinL <- list(top100_protein = top100_protein)
gset_top100_protein <- gsva(data0, top100_proteinL, method='ssgsea', 
                            parallel.sz = 0, min.sz = 5,
                            BPPARAM = BPPARAM)
# tensorL <- list(tensor = tensor)
# gset_tensor <- gsva(data0, tensorL, method='ssgsea',
#                     parallel.sz = 0, min.sz = 5,
#                     BPPARAM = BPPARAM)
# gset_fe <- gsva(data0, fe, method='ssgsea',
#                  parallel.sz = 0, min.sz = 5,
#                  BPPARAM = BPPARAM)
gset_imm <- gsva(data0, imml, method='ssgsea',
                 parallel.sz = 0, min.sz = 5,
                 BPPARAM = BPPARAM)
gset_ERstress <- gsva(data0, ERstress, method='ssgsea',
                      parallel.sz = 0, min.sz = 5,
                      BPPARAM = BPPARAM)

HR_TCGA_level <- ifelse(TCGA_HR > 0, 'High', 'Low')
gset_top100_protein_data <- data.frame(t(gset_top100_protein))
gset_top100_protein_data$HR_TCGA_level <- HR_TCGA_level
gset_top100_protein_dataL <- melt(gset_top100_protein_data)

# gset_tensor_data <- data.frame(t(gset_tensor))
# gset_tensor_data$HR_TCGA_level <- HR_TCGA_level
# gset_tensor_dataL <- melt(gset_tensor_data)
# 
# gset_fe_data <- data.frame(t(gset_fe))
# gset_fe_data$HR_TCGA_level <- HR_TCGA_level
# gset_fe_dataL <- melt(gset_fe_data)

gset_imm_data <- data.frame(t(gset_imm))
gset_imm_data$HR_TCGA_level <- HR_TCGA_level
gset_imm_dataL <- melt(gset_imm_data)

gset_ERstress_data <- data.frame(t(gset_ERstress))
gset_ERstress_data$HR_TCGA_level <- HR_TCGA_level
gset_ERstress_dataL <- melt(gset_ERstress_data)


# ----------------------- 盒子图 ----------------------
gset_top100_protein_dataL$HR_TCGA_level <- factor(gset_top100_protein_dataL$HR_TCGA_level,
                                                  levels = unique(gset_top100_protein_dataL$HR_TCGA_level))
library(ggpubr)
top100box <- ggplot(gset_top100_protein_dataL, aes(x = variable, y = value,
                                                   fill = HR_TCGA_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(top100box, filename = 'TCGAnext/top100box.pdf', height = 7,width = 8)

# gset_tensor_dataL$HR_TCGA_level <- factor(gset_tensor_dataL$HR_TCGA_level,
#                                           levels = unique(gset_tensor_dataL$HR_TCGA_level))
# tensorbox <- ggplot(gset_tensor_dataL, aes(x = variable, y = value,
#                                            fill = HR_TCGA_level)) + 
#   geom_boxplot() +   
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
#                      method.args = list(alternative = "two.sided"))
# ggsave(tensorbox, filename = 'TCGAnext/tensorbox.pdf', height = 7,width = 8)

# gset_fe_dataL$HR_TCGA_level <- factor(gset_fe_dataL$HR_TCGA_level,
#                                        levels = unique(gset_fe_dataL$HR_TCGA_level))
# febox <- ggplot(gset_fe_dataL, aes(x = variable, y = value,
#                                      fill = HR_TCGA_level)) + 
#   geom_boxplot() +   
#   theme_bw() + 
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
#                      method.args = list(alternative = "two.sided"))
# ggsave(febox, filename = 'TCGAnext/febox.pdf', height = 7,width = 8)


gset_imm_dataL$HR_TCGA_level <- factor(gset_imm_dataL$HR_TCGA_level,
                                       levels = unique(gset_imm_dataL$HR_TCGA_level))
immbox <- ggplot(gset_imm_dataL, aes(x = variable, y = value,
                                     fill = HR_TCGA_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
immbox
ggsave(immbox, filename = 'TCGAnext/immubox.pdf', height = 7,width = 8)

gset_ERstress_dataL$HR_TCGA_level <- factor(gset_ERstress_dataL$HR_TCGA_level,
                                            levels = unique(gset_ERstress_dataL$HR_TCGA_level))
ERbox <- ggplot(gset_ERstress_dataL, aes(x = variable, y = value,
                                         fill = HR_TCGA_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))


ggsave(ERbox, filename = 'TCGAnext/ERbox.pdf', height = 20,width = 10)

# ----------------------- 打分热图 + 基因热图 -----------------------
od <- order(TCGA_HR)
HR_TCGA_level_order <- HR_TCGA_level[od]

mxgene_order <- TCGA[rownames(TCGA) %in% unique(top100_protein),od]
an <- HeatmapAnnotation(group = HR_TCGA_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
mxgene_order_scale <- t(apply(mxgene_order, 1, scale))
pdf('TCGAnext/TCGA_heatmap_top100_gene.pdf', width = 7,height = 15)
Heatmap(mxgene_order_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_order <- TCGA[rownames(TCGA) %in% unique(tensor),od]
an <- HeatmapAnnotation(group = HR_TCGA_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
mxgene_order_scale <- t(apply(mxgene_order, 1, scale))
pdf('TCGAnext/TCGA_heatmap_tensor_gene.pdf', width = 7,height = 8)
Heatmap(mxgene_order_scale, cluster_columns = F, top_annotation = an)
dev.off()

od <- order(TCGA_HR)
HR_TCGA_level_order <- HR_TCGA_level[od]


gset_fe_order <- gset_fe[,od]
mxgene_order <- TCGA[rownames(TCGA) %in% unique(c(fe$driver,fe$suppressor)),od]
an <- HeatmapAnnotation(group = HR_TCGA_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
gset_fe_scale <- t(apply(gset_fe_order, 1, scale))
pdf('TCGAnext/TCGA_heatmap_fe_cell.pdf', width = 7,height = 4)
Heatmap(gset_fe_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_scale <- t(apply(mxgene_order, 1, scale))
pdf('TCGAnext/TCGA_heatmap_fe_gene.pdf', width = 7,height = 100)
Heatmap(mxgene_scale, cluster_columns = F, top_annotation = an)
dev.off()


gset_imm_order <- gset_imm[,od]
mxgene_order <- TCGA[rownames(TCGA) %in% unique(imm$Metagene),od]
an <- HeatmapAnnotation(group = HR_TCGA_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
gset_imm_scale <- t(apply(gset_imm_order, 1, scale))
pdf('TCGAnext/TCGA_heatmap_immu_cell.pdf', width = 7,height = 4)
Heatmap(gset_imm_scale, cluster_columns = F, top_annotation = an)
dev.off()

mxgene_scale <- t(apply(mxgene_order, 1, scale))
pdf('TCGAnext/TCGA_heatmap_immu_gene.pdf', width = 7,height = 100)
Heatmap(mxgene_scale, cluster_columns = F, top_annotation = an)
dev.off()

ER_order <- TCGA[rownames(TCGA) %in% unique(ERstressL$gene),od]
an <- HeatmapAnnotation(group = HR_TCGA_level_order,
                        col = list(
                          group = c('Low' = "blue",
                                    "High" = 'red')
                        ))
ER_scale <- t(apply(ER_order, 1, scale))
pdf('TCGAnext/TCGA_heatmap_ER_gene.pdf', width = 7,height = 150)
Heatmap(ER_scale, cluster_columns = F, top_annotation = an)
dev.off()


# -------------------------- 通路相关性图 -----------------------
gset_top100_protein_data <- data.frame(t(gset_top100_protein))
gset_top100_protein_cdata <- apply(gset_top100_protein_data, 2, function(x) cor(x, TCGA_HR))
Cor_top100 <- 
  data.frame(cor = gset_top100_protein_cdata, 
             absCor = abs(unlist(gset_top100_protein_cdata)),
             name = colnames(gset_top100_protein_data))
gset_top100_protein_data$HR <- TCGA_HR

dir.create('TCGAnext/ssgseafe')
for(i in 1:1){
  ns <- colnames(gset_top100_protein_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_top100_protein_data, 
                aes(x = eval(parse(text = colnames(gset_top100_protein_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) +
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('TCGAnext/ssgseafe/',ns,'.pdf'), 
         height = 7, width = 7)
}

# gset_tensor_data <- data.frame(t(gset_tensor))
# gset_tensor_cdata <- apply(gset_tensor_data, 2, function(x) cor(x, TCGA_HR))
# Cor_tensor <- 
#   data.frame(cor = gset_tensor_cdata, 
#              absCor = abs(unlist(gset_tensor_cdata)),
#              name = colnames(gset_tensor_data))
# gset_tensor_data$HR <- TCGA_HR
# 
# for(i in 1:1){
#   ns <- colnames(gset_tensor_data)[i]
#   my.formula <- y ~ x
#   fig <- ggplot(gset_tensor_data, 
#                 aes(x = eval(parse(text = colnames(gset_tensor_data)[i])), y = HR)) + 
#     geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
#     stat_poly_eq(formula = my.formula, 
#                  aes(label = paste(..rr.label.., sep = "~~~")), 
#                  parse = TRUE) +
#     geom_point() + theme_bw() + xlab(ns)
#   ggsave(fig, filename = paste0('TCGAnext/ssgseafe/',ns,'.pdf'), 
#          height = 7, width = 7)
# }


# gset_fe_data <- data.frame(t(gset_fe))
# data_gset_fe_cdata <- apply(gset_fe_data, 2, function(x) cor(x, TCGA_HR))
# Cor_fe <- 
#   data.frame(cor = data_gset_fe_cdata, 
#              absCor = abs(unlist(data_gset_fe_cdata)),
#              name = colnames(gset_fe_data))
# gset_fe_data$HR <- TCGA_HR
# dir.create('TCGAnext/ssgseafe')
# for(i in 1:2){
#   ns <- colnames(gset_fe_data)[i]
#   my.formula <- y ~ x
#   fig <- ggplot(gset_fe_data, 
#                 aes(x = eval(parse(text = colnames(gset_fe_data)[i])), y = HR)) + 
#     geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
#     stat_poly_eq(formula = my.formula, 
#                  aes(label = paste(..rr.label.., sep = "~~~")), 
#                  parse = TRUE) +
#     geom_point() + theme_bw() + xlab(ns)
#   ggsave(fig, filename = paste0('TCGAnext/ssgseafe/',ns,'.pdf'), 
#          height = 7, width = 7)
# }

gset_imm_data <- data.frame(t(gset_imm))
data_gset_imm_cdata <- apply(gset_imm_data, 2, function(x) cor(x, TCGA_HR))
Cor_imm <- 
  data.frame(cor = data_gset_imm_cdata, 
             absCor = abs(unlist(data_gset_imm_cdata)),
             name = colnames(gset_imm_data))
gset_imm_data$HR <- TCGA_HR
dir.create('TCGAnext/ssgseaimmu')
for(i in 1:28){
  ns <- colnames(gset_imm_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_imm_data, 
                aes(x = eval(parse(text = colnames(gset_imm_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) +
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('TCGAnext/ssgseaimmu/',ns,'.pdf'), 
         height = 7, width = 7)
}

gset_ERstress_data <- data.frame(t(gset_ERstress))
data_gset_ERstress_cdata <- apply(gset_ERstress_data, 2, function(x) cor(x, TCGA_HR))
Cor_ERstress <- 
  data.frame(cor = data_gset_ERstress_cdata, 
             absCor = abs(unlist(data_gset_ERstress_cdata)),
             name = colnames(gset_ERstress_data))
gset_ERstress_data$HR <- TCGA_HR
dir.create('TCGAnext/ssgseaERstress')
for(i in 1:23){
  ns <- colnames(gset_ERstress_data)[i]
  my.formula <- y ~ x
  fig <- ggplot(gset_ERstress_data, 
                aes(x = eval(parse(text = colnames(gset_ERstress_data)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) + 
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('TCGAnext/ssgseaERstress/',ns,'.pdf'), 
         height = 7, width = 7)
}



# Cor_fe <- Cor_fe[order(Cor_fe$cor),]
# Cor_fe$name <- factor(Cor_fe$name, levels = Cor_fe$name)
# 
# Cor_feC <- ggplot(Cor_fe, aes(x=name, y=cor, color = absCor)) +
#   geom_point(aes(size = absCor)) + 
#   geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
#   scale_color_gradient2(mid = 'green', high = 'red')+
#   theme_bw() + coord_flip()
# ggsave(Cor_feC, filename ='TCGAnext/Cor_feC.pdf', 
#        height = 7, width = 7)

Cor_imm <- Cor_imm[order(Cor_imm$cor),]
Cor_imm$name <- factor(Cor_imm$name, levels = Cor_imm$name)

Cor_immC <- ggplot(Cor_imm, aes(x=name, y=cor, color = absCor)) +
  geom_point(aes(size = absCor)) + 
  geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
  scale_color_gradient2(mid = 'green', high = 'red')+
  theme_bw() + coord_flip()
ggsave(Cor_immC, filename ='TCGAnext/Cor_immC.pdf', 
       height = 7, width = 7)

Cor_ERstress <- Cor_ERstress[order(Cor_ERstress$cor),]
Cor_ERstress$name <- factor(Cor_ERstress$name, levels = Cor_ERstress$name)
Cor_ERstress2 <- ggplot(Cor_ERstress, aes(x=name, y=cor, color = absCor)) +
  geom_point(aes(size = absCor)) + 
  geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
  scale_color_gradient2(mid = 'green', high = 'red')+
  theme_bw() + coord_flip()
ggsave(Cor_ERstress2, filename = 'TCGAnext/Cor_ERstress2.pdf',
       height = 7, width = 15)
# ======================== TCGA suppler ==========================
rm(list = ls())
TCGA_HR <- readRDS('TCGA_HR.Rds')
immue <- read.csv('data/immsuppler.csv',stringsAsFactors = F)
interpatient <- intersect(immue$TCGA.Participant.Barcode,names(TCGA_HR))
TCGA_HR <- TCGA_HR[interpatient]
rownames(immue) <- immue$TCGA.Participant.Barcode
immue <- immue[interpatient,]
immue <- immue[,-1]
immue_scale <- apply(immue, 2, scale)
immue_scale <- as.data.frame(immue_scale)
dir.create('TCGAimmsuppler')
HR_TCGA_level <- ifelse(TCGA_HR > 0, 'High', 'Low')

gset_imm_data <- immue_scale
gset_imm_data$HR_TCGA_level <- HR_TCGA_level
gset_imm_dataL <- melt(gset_imm_data)
gset_imm_dataL$HR_TCGA_level <- factor(gset_imm_dataL$HR_TCGA_level,
                                       levels = c('High','Low'))
immbox <- ggplot(gset_imm_dataL, aes(x = variable, y = value,
                                     fill = HR_TCGA_level)) + 
  geom_boxplot() +   
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
# ggsave(immbox, filename = 'TCGAimmsuppler/immubox.pdf', height = 7,width = 15)


data_gset_imm_cdata <- apply(immue_scale, 2, function(x) cor.test(x, 
                                                                  TCGA_HR,
                                                                  method = "spearman"
)$estimate)
Cor_imm <- 
  data.frame(cor = data_gset_imm_cdata, 
             absCor = abs(unlist(data_gset_imm_cdata)),
             name = colnames(immue_scale))
immue$HR <- TCGA_HR
dir.create('TCGAimmsuppler/ssgseaimmu')
for(i in 1:56){
  ns <- colnames(immue)[i]
  my.formula <- y ~ x
  fig <- ggplot(immue, 
                aes(x = eval(parse(text = colnames(immue)[i])), y = HR)) + 
    geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
    stat_poly_eq(formula = my.formula, 
                 aes(label = paste(..p.value.label.., sep = "~~~")), 
                 parse = TRUE) +
    geom_point() + theme_bw() + xlab(ns)
  ggsave(fig, filename = paste0('TCGAimmsuppler/ssgseaimmu/',ns,'.pdf'), 
         height = 7, width = 7)
}

Cor_imm <- Cor_imm[order(Cor_imm$cor),]
Cor_imm$name <- factor(Cor_imm$name, levels = Cor_imm$name)
Cor_imm0 <- ggplot(Cor_imm, aes(x=name, y=cor, color = absCor)) +
  geom_point(aes(size = absCor)) + 
  geom_segment( aes(x=name, xend=name, y=0, yend=cor, color = absCor)) + 
  scale_color_gradient2(mid = 'green', high = 'red')+
  theme_bw() + coord_flip()
ggsave(Cor_imm0, filename = 'TCGAimmsuppler/Cor_imm.pdf',
       height = 7, width = 7)

# ============================== TCGA mutation =================================
rm(list = ls())
library(maftools)
library(tidyverse)
library(scales)
library(ggpubr)
# library(TCGAbiolinks)
# setwd('data/')
# maf <- GDCquery_Maf("LUAD", pipelines = "muse")
# saveRDS(maf, file = 'maf.Rds')
# setwd('..')
TCGA_HR <- readRDS('TCGA_HR.Rds')
mf <- readRDS('data/maf.Rds')

hl <- names(which(TCGA_HR > 0))
ll <- names(which(TCGA_HR <= 0))

pm <- str_extract(mf$Tumor_Sample_Barcode, pattern = 'TCGA-..-....')
Hpm <- mf[pm %in% hl,]
Lpm <- mf[pm %in% ll,]

Hpm <- read.maf(Hpm)
Lpm <- read.maf(Lpm)
dir.create('TCGA_mutation')
dir.create('TCGA_mutation/high_risk')

# pdf('TCGA_mutation/high_risk/summary.pdf')
# plotmafSummary(maf = Hpm, rmOutlier = F, 
#                addStat = 'median', 
#                dashboard = TRUE, 
#                titvRaw = FALSE)
# dev.off()
# pdf('TCGA_mutation/high_risk/oncoplot.pdf', width = 10, height = 10)
# oncoplot(maf = Hpm, top=50)
# dev.off()
# 
# laml.titv = titv(maf = Hpm, plot = FALSE, useSyn = TRUE)
# 
# pdf('TCGA_mutation/high_risk/TiTv.pdf')
# plotTiTv(res = laml.titv)
# dev.off()
# 
# pdf('TCGA_mutation/high_risk/somaticInteractions.pdf')
# somaticInteractions(maf = Hpm, top = 25, pvalue = c(0.05, 0.1))
# dev.off()
# 
# pdf('TCGA_mutation/high_risk/drugInteractions.pdf')
# dgi = drugInteractions(maf = Hpm, fontSize = 0.75)
# dev.off()
# 
# pdf('TCGA_mutation/high_risk/OncogenicPathways.pdf')
# OncogenicPathways(maf = Hpm)
# dev.off()
# 
# 
# 
# dir.create('TCGA_mutation/low_risk')
# pdf('TCGA_mutation/low_risk/summary.pdf')
# plotmafSummary(maf = Lpm, rmOutlier = F, 
#                addStat = 'median', 
#                dashboard = TRUE, 
#                titvRaw = FALSE)
# dev.off()
# pdf('TCGA_mutation/low_risk/oncoplot.pdf', width = 10, height = 10)
# oncoplot(maf = Lpm, top=50)
# dev.off()
# 
# laml.titv = titv(maf = Lpm, plot = FALSE, useSyn = TRUE)
# 
# pdf('TCGA_mutation/low_risk/TiTv.pdf')
# plotTiTv(res = laml.titv)
# dev.off()
# 
# pdf('TCGA_mutation/low_risk/somaticInteractions.pdf')
# somaticInteractions(maf = Lpm, top = 25, pvalue = c(0.05, 0.1))
# dev.off()
# 
# pdf('TCGA_mutation/low_risk/drugInteractions.pdf')
# dgi = drugInteractions(maf = Lpm, fontSize = 0.75)
# dev.off()
# 
# pdf('TCGA_mutation/low_risk/OncogenicPathways.pdf')
# OncogenicPathways(maf = Lpm)
# dev.off()
# 
# pdf('TCGA_mutation/RTK-RAS_p53_patheay.pdf')
# PlotOncogenicPathways(maf = Hpm, pathways = "RTK-RAS")
# PlotOncogenicPathways(maf = Lpm, pathways = "RTK-RAS")
# PlotOncogenicPathways(maf = Hpm, pathways = "TP53")
# PlotOncogenicPathways(maf = Lpm, pathways = "TP53")
# dev.off()
# 
# pdf('TCGA_mutation/p53.pdf')
# lollipopPlot(
#   maf = Hpm,
#   gene = 'TP53',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot2(
#   maf = Lpm,
#   gene = 'TP53',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot2(Hpm, Lpm,'TP53')
# dev.off()
# 
# pdf('TCGA_mutation/Kras.pdf')
# lollipopPlot(
#   maf = Hpm,
#   gene = 'KRAS',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot(
#   maf = Lpm,
#   gene = 'KRAS',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot2(Hpm, Lpm,'KRAS')
# dev.off()
# 
# pdf('TCGA_mutation/CSMD3.pdf')
# lollipopPlot(
#   maf = Hpm,
#   gene = 'CSMD3',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot(
#   maf = Lpm,
#   gene = 'CSMD3',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot2(Hpm, Lpm,'CSMD3')
# dev.off()
# pdf('TCGA_mutation/MUC16.pdf')
# lollipopPlot(
#   maf = Hpm,
#   gene = 'MUC16',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot(
#   maf = Lpm,
#   gene = 'MUC16',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot2(Hpm, Lpm,'MUC16')
# dev.off()
# pdf('TCGA_mutation/RYR2.pdf')
# lollipopPlot(
#   maf = Hpm,
#   gene = 'RYR2',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot(
#   maf = Lpm,
#   gene = 'RYR2',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot2(Hpm, Lpm,'RYR2')
# dev.off()
# pdf('TCGA_mutation/LRP1B.pdf')
# lollipopPlot(
#   maf = Hpm,
#   gene = 'LRP1B',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot(
#   maf = Lpm,
#   gene = 'LRP1B',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot2(Hpm, Lpm,'LRP1B')
# dev.off()
# pdf('TCGA_mutation/RB1.pdf')
# lollipopPlot(
#   maf = Hpm,
#   gene = 'KRAS',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot(
#   maf = Lpm,
#   gene = 'KRAS',
#   showMutationRate = TRUE,
#   printCount = T,
#   defaultYaxis = T,
#   cBioPortal = T
# )
# lollipopPlot2(Hpm, Lpm,'KRAS')
# dev.off()

pointL <- c('G12A','G12C','G12V','G12D','G12S','G13C','G13D')
numbL <- c('6','24','10','4','1','1','1')
names(numbL) <- pointL

pointH <- c('G12V','G12C','G12D','G12A','G12S','G13C')
numbH <- c('11','19','6','7','3','2','1')
names(numbH) <- pointH

pointt <- unique(c(pointL, pointH))
numbL <- numbL[pointt]
numbH <- numbH[pointt]

numbT <- data.frame(numbL = as.numeric(numbL),
                    numbH = as.numeric(numbH))
numbT[is.na(numbT)] <- 0
rownames(numbT) <- pointt

fisher.test(numbT)
binom.test(numbT[5,1],rowSums(numbT)[5])
sum(colSums(numbT)) -> tttt

pt.vs.rt <- mafCompare(m1 = Hpm, m2 = Lpm, 
                       m1Name = 'High',
                       m2Name = 'Low',
                       pathways = T,
                       minMut = 5)
View(pt.vs.rt)
pt.vs.rt$results$Hugo_Symbol -> pw
# write.csv(pt.vs.rt$results,file = 'pathwaycompair.csv')
pt.vs.rtg <- mafCompare(m1 = Hpm, m2 = Lpm, 
                        m1Name = 'High',
                        m2Name = 'Low',
                        pathways = F,
                        minMut = 5)
# write.csv(pt.vs.rtg$results,file = 'genescompair.csv')

pw
# pdf('oncopathwayNOTCH.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[1]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[1]))
# dev.off()
# pdf('oncopathwayTP53.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[2]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[2]))
# dev.off()
# pdf('oncopathwayNRF2.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[3]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[3]))
# dev.off()
# pdf('oncopathwayRTK-RAS.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[4]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[4]))
# dev.off()
# pdf('oncopathwayWNT.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[5]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[5]))
# dev.off()
# pdf('oncopathwayMYC.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[6]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[6]))
# dev.off()
# pdf('oncopathwayHippo.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[7]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[7]))
# dev.off()
# pdf('oncopathwayPI3K.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[8]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[8]))
# dev.off()
# pdf('oncopathwayCell_Cycle.pdf')
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[9]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[9]))
# dev.off()
# pdf('oncopathwayTGF-Beta.pdf',width = 3,height = 3)
# print(PlotOncogenicPathways(maf = Hpm, pathways = pw[10]))
# print(PlotOncogenicPathways(maf = Lpm, pathways = pw[10]))
# dev.off()
# pdf('kras_p53_all.pdf',width = 7, height = 3)
# oncostrip(maf = Hpm, genes = c('TP53','KRAS'))
# oncostrip(maf = Lpm, genes = c('TP53','KRAS'))
# dev.off()






Hpmx <- mf[pm %in% hl & 
             mf$Hugo_Symbol %in% c('TP53','KRAS') &
             (mf$Hugo_Symbol == 'TP53' | 
                str_detect(string = mf$HGVSp_Short,pattern = "p.G12")),]
Lpmx <- mf[pm %in% ll & 
             mf$Hugo_Symbol %in% c('TP53','KRAS') &
             (mf$Hugo_Symbol == 'TP53' | 
                str_detect(string = mf$HGVSp_Short,pattern = "p.G12")),]
Hpmx <- read.maf(Hpmx)
Lpmx <- read.maf(Lpmx)
# pdf('KRAS_hotspot.pdf')
# lollipopPlot2(Hpmx, Lpmx,'KRAS')
# dev.off()

# pdf('kras_p53_hotspot.pdf',width = 7, height = 3)
# oncoplot(maf = Hpmx, top=50)
# oncoplot(maf = Lpmx, top=50)
# dev.off()

tmpd <- data.frame(
  H = as.numeric(Hpmx@gene.summary[1,]),
  L = as.numeric(Lpmx@gene.summary[1,]))
tmpd <- tmpd[-1,]
rownames(tmpd) <- names(Lpmx@gene.summary)[2:7]
pl53 <- list()
for(i in 1:6){
  binom.test(tmpd[i,1],rowSums(tmpd)[i]) -> tmp0
  pl53[[i]] <- tmp0$p.value
}
names(pl53) <- names(Lpmx@gene.summary)[2:7]
pl53 <- unlist(pl53)
tmpd$names <- rownames(tmpd)
tmpd$p <- pl53
tmpd$p <- signif(tmpd$p,3)
library(reshape2)
tmpdl <- melt(tmpd,c('names', 'p'))
pl53txt <- data.frame(p = pl53,
                      names = names(pl53))

ggplot(tmpdl,aes(x = names, y = value, 
                 group = variable, 
                 color = variable,
                 shape = variable)) +
  geom_point(size = 3, alpha = 0.9) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  ylim(0,150) + geom_text(aes(x = names, 
                              y = 150, 
                              label = p),
                          color = 'black',size = 2) -> p1
tmpdl <- tmpdl[order(tmpdl$names),]
tmpdl[tmpdl$variable == 'L',4]/(tmpdl %>% 
                                  group_by(names) %>% 
                                  summarise(sum(value)) %>% 
                                  pull(`sum(value)`)) -> pct


data.frame(x = unique(tmpdl$names),
           High = rep(1,times = 6),
           Low = pct) -> d0
p2 <- ggplot(d0) + geom_bar(
  mapping = aes(x = x, y = High,
                fill = hue_pal()(2)[2]),
  stat = 'identity',
  color = 'black') + 
  geom_bar(
    mapping = aes(x = x, y = Low,
                  fill = hue_pal()(2)[1]),
    stat = 'identity',
    color = 'black') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))+ 
  geom_hline(
    yintercept = 0.5,
    color = 'red',
    linetype = 2,
    size = 1)
# ==================== 折现 =============================
# saveRDS(list(p1,p2),file = 'p53pointline.Rds')
fp <- readRDS('p53pointline.Rds')
library(cowplot)
library(scales)

pdf('p53pointline.pdf',height = 3,width = 4)
fp[[1]]
fp[[2]]
dev.off()

tmpd <- data.frame(
  H = as.numeric(Hpmx@gene.summary[2,]),
  L = as.numeric(Lpmx@gene.summary[2,]))
tmpd <- tmpd[-1,]
plkras <- list()
for(i in 1:6){
  if(rowSums(tmpd)[i] != 0) binom.test(tmpd[i,1],rowSums(tmpd)[i]) -> tmp0
  plkras[[i]] <- tmp0$p.value
}
names(plkras) <- names(Lpmx@gene.summary)[2:7]
plkras <- unlist(plkras)
tmpd$names <- rownames(tmpd)
tmpd$p <- plkras
tmpd$p <- signif(tmpd$p,3)
library(reshape2)
tmpd$names <- names(plkras)
tmpdl <- melt(tmpd,c('names', 'p'))

ggplot(tmpdl,aes(x = names, y = value, 
                 group = variable, 
                 color = variable,
                 shape = variable)) +
  geom_point(size = 3, alpha = 0.9) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  ylim(0,100) + geom_text(aes(x = names, 
                              y = 100, 
                              label = p),
                          color = 'black',size = 2) -> p1
tmpdl <- tmpdl[order(tmpdl$names),]
tmpdl[tmpdl$variable == 'L',4]/(tmpdl %>% 
                                  group_by(names) %>% 
                                  summarise(sum(value)) %>% 
                                  pull(`sum(value)`)) -> pct


data.frame(x = unique(tmpdl$names),
           High = rep(1,times = 6),
           Low = pct) -> d0
p2 <- ggplot(d0) + geom_bar(
  mapping = aes(x = x, y = High,
                fill = hue_pal()(2)[2]),
  stat = 'identity',
  color = 'black') + 
  geom_bar(
    mapping = aes(x = x, y = Low,
                  fill = hue_pal()(2)[1]),
    stat = 'identity',
    color = 'black') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))+ 
  geom_hline(
    yintercept = 0.5,
    color = 'red',
    linetype = 2,
    size = 1)
# ==================== 折现 =============================
# saveRDS(list(p1,p2),file = 'KRASpointline.Rds')
fp <- readRDS('KRASpointline.Rds')
pdf('KRASpointline.pdf',height = 3,width = 4)
fp[[1]]
fp[[2]]
dev.off()



distinct(
  data.frame(
    sample = Hpmx@data$Tumor_Sample_Barcode,
    gene = Hpmx@data$Hugo_Symbol
  )
) -> hgene
KRAS <- table(hgene)[,"KRAS"]
TP53 <- table(hgene)[,"TP53"]

table(
  data.frame(KRAS = KRAS,
             TP53 = TP53
  )
) -> tHpmx

distinct(
  data.frame(
    sample = Lpmx@data$Tumor_Sample_Barcode,
    gene = Lpmx@data$Hugo_Symbol
  )
) -> Lgene
KRAS <- table(Lgene)[,"KRAS"]
TP53 <- table(Lgene)[,"TP53"]

table(
  data.frame(KRAS = KRAS,
             TP53 = TP53
  )
) -> tLgene

Highd = melt(tHpmx)
Lowd = melt(tLgene)
cbind(Highd, Low = Lowd[,3]) -> bind
colnames(bind)[3] <- 'High'
plnb <- list()
for(i in 1:4){
  if(rowSums(bind)[i] != 0) binom.test(bind[i,1],rowSums(bind)[i]) -> tmp0
  plnb[[i]] <- tmp0$p.value
}
bind$p <- unlist(plnb)
bind$p <- signif(bind$p,3)
bind$names <- c('KRAS12-/TP53-',
                'KRAS12+/TP53-',
                'KRAS12-/TP53+',
                'KRAS12+/TP53+')
bind <- bind[,3:6]
bind <- melt(bind, c('names','p'))

ggplot(bind,aes(x = names, y = value, 
                group = variable, 
                color = variable,
                shape = variable)) +
  geom_point(size = 3, alpha = 0.9) + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  ylim(0,150) + geom_text(aes(x = names, 
                              y = 150, 
                              label = p),
                          color = 'black',size = 2) -> p1
bind <- bind[order(bind$names),]
bind[bind$variable == 'Low',4]/(bind %>% 
                                  group_by(names) %>% 
                                  summarise(sum(value)) %>% 
                                  pull(`sum(value)`)) -> pct


data.frame(x = unique(bind$names),
           High = rep(1,times = 4),
           Low = pct) -> d0
p2 <- ggplot(d0) + geom_bar(
  mapping = aes(x = x, y = High,
                fill = hue_pal()(2)[2]),
  stat = 'identity',
  color = 'black') + 
  geom_bar(
    mapping = aes(x = x, y = Low,
                  fill = hue_pal()(2)[1]),
    stat = 'identity',
    color = 'black') +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1))+ 
  geom_hline(
    yintercept = 0.5,
    color = 'red',
    linetype = 2,
    size = 1)
# ==================== 折现 =============================
# saveRDS(list(p1,p2),file = 'KRAS12_p53_linepoint.Rds')
pf <- readRDS('KRAS12_p53_linepoint.Rds')
pdf('KRAS12_p53_linepoint.pdf',height = 3,width = 4)
pf[[1]]
pf[[2]]
dev.off()



data.frame(
  H = as.numeric(table(rowSums(table(hgene)))),
  L = as.numeric(table(rowSums(table(Lgene))))
) -> mutn
plnb <- list()
for(i in 1:3){
  if(rowSums(mutn)[i] != 0) binom.test(mutn[i,1],rowSums(mutn)[i]) -> tmp0
  plnb[[i]] <- tmp0$p.value
}
mutn$p <- unlist(plnb)
mutn$names <- 0:2
mutnL <- melt(mutn,c('names', 'p'))
mutnL$p <- signif(mutnL$p,3)
mutnL$names <- as.factor(mutnL$names)
pdf('KRAS_p53_linepoint.pdf',height = 4,width = 4)
ggplot(mutnL,aes(x = names, y = value, 
                 group = variable, color = variable)) +
  geom_point() + geom_line() + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  ylim(0,200) + geom_text(aes(x = names, 
                              y = 200, 
                              label = p),
                          color = 'black',size = 4)
dev.off()


# PlotOncogenicPathways(maf = Hpm, pathways = "RTK-RAS")
# rainfallPlot(maf = Hpm, detectChangePoints = TRUE, pointSize = 0.4)
# lollipopPlot(
#   maf = Hpm,
#   gene = 'TP53',
#   showMutationRate = TRUE,
#   labelPos = 882
# )
# plotProtein(gene = "TP53", refSeqID = "NM_000546")
# laml.mutload = tcgaCompare(maf = Hpm, cohortName = 'Example-LAML', logscale = TRUE, capture_size = 50)
# plotVaf(maf = Hpm)
# laml.sig = oncodrive(maf = Hpm)
# plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)
# laml.pfam = pfamDomains(maf = Hpm, top = 10)
# dnmt3a.dgi = drugInteractions(genes = "DNMT3A", drugs = TRUE)
Hpmx <- mf[pm %in% hl & 
             mf$Hugo_Symbol %in% c('TP53','KRAS') &
             (mf$Hugo_Symbol == 'TP53' | 
                str_detect(string = mf$HGVSp_Short,pattern = "p.G12")),]
Lpmx <- mf[pm %in% ll & 
             mf$Hugo_Symbol %in% c('TP53','KRAS') &
             (mf$Hugo_Symbol == 'TP53' | 
                str_detect(string = mf$HGVSp_Short,pattern = "p.G12")),]
Hpmx <- read.maf(Hpmx)
Lpmx <- read.maf(Lpmx)
# =============== 比较 ===============
pt.vs.rt <- mafCompare(m1 = Hpm, m2 = Lpm, m1Name = 'High_risk', 
                       m2Name = 'Low_risk', minMut = 5)
print(pt.vs.rt)
pdf('TCGA_mutation/forestPlot.pdf',height = 10)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.01)
dev.off()
pt.vs.rt[["results"]]-> cmp
cmp[cmp$Hugo_Symbol %in% c('CSMD3',
                           'MUC16',
                           'RYR2',
                           'LRP1B'),]



pdf('TCGA_mutation/coOncoplot.pdf')
coOncoplot(m1 = Hpm, m2 = Lpm,
           m1Name = 'High_risk', m2Name = 'Low_risk', 
           removeNonMutated = TRUE)
dev.off()

pdf('TCGA_mutation/coBarplot.pdf')
coBarplot(m1 = Hpm, m2 = Lpm, m1Name = "High_risk", 
          m2Name = "Low_risk")
dev.off()

mf0 <- read.maf(mf)
t0 <- data.frame(tmb = tmb(mf0)[,3], 
                 barcode = unlist(str_extract(mf0@clinical.data$Tumor_Sample_Barcode, 
                                              pattern = '....-..-....')
                 )
)

t0$HR <- TCGA_HR[t0$barcode]
t0 <- na.omit(t0)
t0$HR_level <- ifelse(t0$HR > 0, 'High','Low')
mycomparision <- list(c('High','Low'))
bx <- ggplot(t0, aes(x = HR_level, y = total_perMB, fill = HR_level)) + 
  geom_boxplot() + theme_bw() + 
  stat_compare_means(method = 'wilcox.test', comparisons = mycomparision, label = 'p.signif')
ggsave(bx, filename = 'TCGA_mutation/TMB_BOX.pdf')

my.formula <- y ~ x
pp <- ggplot(t0, aes(x = HR, y = total_perMB))+ 
  geom_smooth(method = "lm", se=T, color="black", formula = my.formula) +
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..rr.label.., sep = "~~~")), 
               parse = TRUE) + 
  geom_point() + theme_bw() 
ggsave(pp, filename = 'TCGA_mutation/TMB_point.pdf')



# > summary(res.cox)
# Call:
#   coxph(formula = Surv(time, status) ~ ., data = lung)
# 
# n= 520, number of events= 187 
# 
# coef exp(coef) se(coef)      z Pr(>|z|)  
# PLEK2      0.08023   1.08354  0.10568  0.759   0.4477  
# PTPRH      0.08992   1.09409  0.08929  1.007   0.3139  
# OGFRP1     0.15723   1.17027  0.07651  2.055   0.0399 *
#   CHRNA5     0.19302   1.21291  0.08414  2.294   0.0218 *
#   CBFA2T3   -0.03516   0.96545  0.08641 -0.407   0.6841  
# SMIM15     0.07812   1.08125  0.07996  0.977   0.3286  
# AVEN       0.10300   1.10849  0.08190  1.258   0.2085  
# MELTF      0.12072   1.12831  0.10183  1.185   0.2358  
# KRT8       0.06052   1.06239  0.09133  0.663   0.5076  
# RGS20      0.03030   1.03076  0.08213  0.369   0.7122  
# FAM207A    0.12458   1.13268  0.09576  1.301   0.1933  
# SOWAHC     0.11292   1.11955  0.08147  1.386   0.1657  
# ELF5      -0.08243   0.92088  0.08825 -0.934   0.3503  
# LSP1P4     0.14295   1.15368  0.07854  1.820   0.0687 .
# C20orf197 -0.07075   0.93170  0.09542 -0.741   0.4585  
# C11orf16  -0.11591   0.89056  0.09411 -1.232   0.2181  
# DNALI1    -0.03104   0.96944  0.08447 -0.367   0.7133  
# ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# exp(coef) exp(-coef) lower .95 upper .95
# PLEK2        1.0835     0.9229    0.8808     1.333
# PTPRH        1.0941     0.9140    0.9184     1.303
# OGFRP1       1.1703     0.8545    1.0073     1.360
# CHRNA5       1.2129     0.8245    1.0285     1.430
# CBFA2T3      0.9654     1.0358    0.8150     1.144
# SMIM15       1.0813     0.9249    0.9244     1.265
# AVEN         1.1085     0.9021    0.9441     1.302
# MELTF        1.1283     0.8863    0.9242     1.378
# KRT8         1.0624     0.9413    0.8883     1.271
# RGS20        1.0308     0.9702    0.8775     1.211
# FAM207A      1.1327     0.8829    0.9388     1.367
# SOWAHC       1.1195     0.8932    0.9543     1.313
# ELF5         0.9209     1.0859    0.7746     1.095
# LSP1P4       1.1537     0.8668    0.9891     1.346
# C20orf197    0.9317     1.0733    0.7728     1.123
# C11orf16     0.8906     1.1229    0.7405     1.071
# DNALI1       0.9694     1.0315    0.8215     1.144
# 
# Concordance= 0.679  (se = 0.023 )
# Likelihood ratio test= 88.73  on 17 df,   p=1e-11
# Wald test            = 86.91  on 17 df,   p=2e-11
# Score (logrank) test = 92.53  on 17 df,   p=2e-12
# ============================
# library(GEOquery)
# library(tidyverse)
# v1 <- c('GSE31210', 'GSE50081','GSE68465','GSE72094',
#         'GSE10072','GSE32863','GSE43458')
# 
# dir.create('GEO/')
# dir.create('GEO/matrix')
# 
# setwd('GEO/matrix')
# 
# for(i in v1){
#     t1 <- getGEO(i, destdir ='.')
#     fn <- paste0(i,'.Rds')
#     saveRDS(t1,fn)
# }
# fn <- paste0(v1,'.Rds')
# 
# for(i in fn[2:7]){
#     ti <- readRDS(i)
#     dn <- unlist(str_extract_all(i, pattern = 'GSE[0-9]+'))
#     dir.create(dn)
#     fj <- ti[[1]]$`geo_accession`
#     setwd(paste0('./', dn))
#     for(j in fj){
#         getGEOSuppFiles(j, makeDirectory = F)
#     }
#     setwd('..')
# }
# i='GSE32863.Rds'
# ti <- readRDS(i)
# dn <- unlist(str_extract_all(i, pattern = 'GSE[0-9]+'))
# dir.create(dn)
# fj <- ti[[1]]$`geo_accession`
# setwd(paste0('./', dn))
# for(j in fj){
#     getGEOSuppFiles(j, makeDirectory = F)
# }
# setwd('..')
# ==========================================
# [1] "PLEK2"     "PTPRH"     "OGFRP1"    "CHRNA5"   
# [5] "CBFA2T3"   "SMIM15"    "AVEN"      "MELTF"    
# [9] "KRT8"      "RGS20"     "FAM207A"   "SOWAHC"   
# [13] "ELF5"      "LSP1P4"    "C20orf197" "C11orf16" 

res.cox <- readRDS(file = 'res.cox_model.Rds')
sv <- c('GSE31210', 'GSE50081','GSE68465','GSE72094')
deff <- c('GSE10072','GSE32863','GSE43458')

library('affy')
library(GEOquery) 
library(tidyverse)
library(hgu133plus2cdf)
library(hgu133plus2.db)
library(annotate)
library(ggsci)

cre <- paste0('./GEO/matrix/')
setwd(cre)
# ========================== GSE31210 ============================================
id <- 'GSE31210.Rds'
temp <- readRDS(id)
metadata <- pData(temp[[1]])
expr0 <- exprs(temp[[1]])
gene.symbols <- getSYMBOL(rownames(expr0), "hgu133plus2")

res.cox$coefficients[!names(res.cox$coefficients) %in% gene.symbols]
nona <- !is.na(rownames(expr0))

expr0 <- as.data.frame(expr0[nona,])
gene.symbols <- gene.symbols[nona]

listE <- split(expr0, gene.symbols)
listE <- lapply(listE, function(x) apply(x, 2, mean))
gn <- names(listE)
se <- purrr::reduce(listE, rbind)
rownames(se) <- gn
GSE31210 <- list(se = se,
                 metadata = metadata)
saveRDS(GSE31210, file = 'GSE31210_exprs_metadata.Rds')

# ====================================  GSE50081 ===============================
id <- 'GSE50081.Rds'
temp <- readRDS(id)
metadata <- pData(temp[[1]])
expr0 <- exprs(temp[[1]])
gene.symbols <- getSYMBOL(rownames(expr0), "hgu133plus2")

res.cox$coefficients[!names(res.cox$coefficients) %in% gene.symbols]
nona <- !is.na(rownames(expr0))

expr0 <- as.data.frame(expr0[nona,])
gene.symbols <- gene.symbols[nona]

listE <- split(expr0, gene.symbols)
listE <- lapply(listE, function(x) apply(x, 2, mean))
gn <- names(listE)
se <- purrr::reduce(listE, rbind)
rownames(se) <- gn
GSE50081 <- list(se = se,
                 metadata = metadata)
saveRDS(GSE50081, file = 'GSE50081_exprs_metadata.Rds')
# ==================================== GSE68465 ================================
library(hgu133a.db)
id <- 'GSE68465.Rds'
temp <- readRDS(id)
metadata <- pData(temp[[1]])
expr0 <- exprs(temp[[1]])
gene.symbols <- getSYMBOL(rownames(expr0), "hgu133a")

res.cox$coefficients[!names(res.cox$coefficients) %in% gene.symbols]
# 5
nona <- !is.na(rownames(expr0))

expr0 <- as.data.frame(expr0[nona,])
gene.symbols <- gene.symbols[nona]

listE <- split(expr0, gene.symbols)
listE <- lapply(listE, function(x) apply(x, 2, mean))
gn <- names(listE)
se <- purrr::reduce(listE, rbind)
rownames(se) <- gn
GSE68465 <- list(se = se,
                 metadata = metadata)
saveRDS(GSE68465, file = 'GSE68465_exprs_metadata.Rds')
# =============================== GSE72094 ========================================
id <- 'GSE72094.Rds'
temp <- readRDS(id)
metadata <- pData(temp[[1]])
expr0 <- exprs(temp[[1]])
gene.symbols <- 
  temp[["GSE72094_series_matrix.txt.gz"]]@featureData@data[["GeneSymbol"]]

res.cox$coefficients[!names(res.cox$coefficients) %in% gene.symbols]
# 5
nona <- !is.na(rownames(expr0))

expr0 <- as.data.frame(expr0[nona,])
gene.symbols <- gene.symbols[nona]

listE <- split(expr0, gene.symbols)
listE <- lapply(listE, function(x) apply(x, 2, mean))
gn <- names(listE)
se <- purrr::reduce(listE, rbind)
rownames(se) <- gn
GSE72094 <- list(se = se,
                 metadata = metadata)
saveRDS(GSE72094, file = 'GSE72094_exprs_metadata.Rds')

# ========================= GSE10072 ==========================================
library(hgu133a.db)
id <- 'GSE10072.Rds'
temp <- readRDS(id)
metadata <- pData(temp[[1]])
expr0 <- exprs(temp[[1]])
gene.symbols <- getSYMBOL(rownames(expr0), "hgu133a")

res.cox$coefficients[!names(res.cox$coefficients) %in% gene.symbols]
# 5
nona <- !is.na(rownames(expr0))

expr0 <- as.data.frame(expr0[nona,])
gene.symbols <- gene.symbols[nona]

listE <- split(expr0, gene.symbols)
listE <- lapply(listE, function(x) apply(x, 2, mean))
gn <- names(listE)
se <- purrr::reduce(listE, rbind)
rownames(se) <- gn
GSE10072 <- list(se = se,
                 metadata = metadata)
saveRDS(GSE10072, file = 'GSE10072_exprs_metadata.Rds')
# ===================== GSE32863 =======================================
id <- 'GSE32863.Rds'
temp <- readRDS(id)
metadata <- pData(temp[[1]])
expr0 <- exprs(temp[[1]])
gene.symbols <- temp[["GSE32863_series_matrix.txt.gz"]]@featureData@data[["Symbol"]]

res.cox$coefficients[!names(res.cox$coefficients) %in% gene.symbols]
# 7
nona <- !is.na(rownames(expr0))

expr0 <- as.data.frame(expr0[nona,])
gene.symbols <- gene.symbols[nona]

listE <- split(expr0, gene.symbols)
listE <- lapply(listE, function(x) apply(x, 2, mean))
gn <- names(listE)
se <- purrr::reduce(listE, rbind)
rownames(se) <- gn
GSE32863 <- list(se = se,
                 metadata = metadata)
saveRDS(GSE32863, file = 'GSE32863_exprs_metadata.Rds')
# ====================================== GSE43458 ==================================
library(hugene10sttranscriptcluster.db)
id <- 'GSE43458.Rds'
temp <- readRDS(id)
metadata <- pData(temp[[1]])
expr0 <- exprs(temp[[1]])
gn <- temp[["GSE43458_series_matrix.txt.gz"]]@featureData@data$gene_assignment != '---'
expr0 <- expr0[gn,]
gene.symbols <- temp[["GSE43458_series_matrix.txt.gz"]]@featureData@data$gene_assignment
gene.symbols <- gene.symbols[gn]

gene.symbols <- strsplit(gene.symbols, split = '//')
gene.symbols <- unlist(lapply(gene.symbols, function(x) x[2]))
gene.symbols <- gsub(gene.symbols, pattern = ' ',replacement = '')
res.cox$coefficients[!names(res.cox$coefficients) %in% gene.symbols]
expr0 <- as.data.frame(expr0)
listE <- split(expr0, gene.symbols)
listE <- lapply(listE, function(x) apply(x, 2, mean))
gn <- names(listE)
se <- purrr::reduce(listE, rbind)
rownames(se) <- gn
GSE43458 <- list(se = se,
                 metadata = metadata)
saveRDS(GSE43458, file = 'GSE43458_exprs_metadata.Rds')
# ============================= DEG ==================================================
library(limma)
library(reshape2)
library(ggpubr)

cre <- paste0('./GEO/matrix/')
setwd(cre)
deff <- paste0(deff, '.Rds')
de1 <- readRDS('GSE10072_exprs_metadata.Rds')
table(de1$metadata$source_name_ch1)
dim(de1$se)
names(res.cox$coefficients)[names(res.cox$coefficients) %in% rownames(de1$se)]
d <- data.frame(x = de1$metadata$source_name_ch1,
                PLEK2 = de1$se['PLEK2',],
                PTPRH = de1$se['PTPRH',],
                CHRNA5 = de1$se['CHRNA5',],
                CBFA2T3 = de1$se['CBFA2T3',],
                AVEN = de1$se['AVEN',],
                MELTF = de1$se['MELTF',],
                KRT8 = de1$se['KRT8',],
                RGS20 = de1$se['RGS20',],
                SOWAHC = de1$se['SOWAHC',],
                ELF5 = de1$se['ELF5',],
                C11orf16 = de1$se['C11orf16',],
                DNALI1 = de1$se['DNALI1',])
d <- melt(d)
GSE10072 <- ggplot(d, aes(x = variable, y = value, color = x)) + 
  geom_boxplot() +  
  theme_bw() +
  theme(axis.text.x = 
          element_text(angle = 45, 
                       hjust = 0.5, 
                       vjust = 0.5)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided"))
ggsave(GSE10072, filename = 'GSE10072_deff.pdf', width = 8)


de2 <- readRDS('GSE32863_exprs_metadata.Rds')
table(de2$metadata$source_name_ch1)
table(de2$metadata$characteristics_ch1.1)
de2bc <- str_extract(de2$metadata$characteristics_ch1.1, pattern = '[0-2]')
de2$se <- removeBatchEffect(de2$se, batch = de2bc)
names(res.cox$coefficients)[names(res.cox$coefficients) %in% rownames(de2$se)]
d <- data.frame(x = de2$metadata$source_name_ch1,
                PLEK2 = de2$se['PLEK2',],
                PTPRH = de2$se['PTPRH',],
                CHRNA5 = de2$se['CHRNA5',],
                CBFA2T3 = de2$se['CBFA2T3',],
                AVEN = de2$se['AVEN',],
                KRT8 = de2$se['KRT8',],
                RGS20 = de2$se['RGS20',],
                ELF5 = de2$se['ELF5',],
                C11orf16 = de2$se['C11orf16',],
                DNALI1 = de2$se['DNALI1',])
d$x <- factor(d$x, levels = c("Lung adenocarcinoma", "Adjacent non-tumor lung"))
d <- melt(d)
GSE32863 <- ggplot(d, aes(x = variable, y = value, color = x)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = 
          element_text(angle = 45, 
                       hjust = 0.5, 
                       vjust = 0.5)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided")) 

ggsave(GSE32863, filename = 'GSE32863_deff.pdf', width = 8)

de3 <- readRDS('GSE43458_exprs_metadata.Rds')
dim(de3$se)
table(de3$metadata$characteristics_ch1.1)
names(res.cox$coefficients)[names(res.cox$coefficients) %in% rownames(de3$se)]
d <- data.frame(x = de3$metadata$characteristics_ch1.1,
                PLEK2 = de3$se['PLEK2',],
                PTPRH = de3$se['PTPRH',],
                OGFRP1 = de3$se['OGFRP1',],
                CHRNA5 = de3$se['CHRNA5',],
                CBFA2T3 = de3$se['CBFA2T3',],
                SMIM15 = de3$se['SMIM15',],
                AVEN = de3$se['AVEN',],
                KRT8 = de3$se['KRT8',],
                RGS20 = de3$se['RGS20',],
                FAM207A = de3$se['FAM207A',],
                SOWAHC = de3$se['SOWAHC',],
                ELF5 = de3$se['ELF5',],
                C20orf197 = de3$se['C20orf197',],
                C11orf16 = de3$se['C11orf16',],
                DNALI1 = de3$se['DNALI1',])
d <- melt(d)
GSE43458 <- ggplot(d, aes(x = variable, y = value, color = x)) + 
  geom_boxplot() + theme_bw() +
  theme(axis.text.x = 
          element_text(angle = 45, 
                       hjust = 0.5, 
                       vjust = 0.5)) + 
  stat_compare_means(aes(label =..p.signif..), method = "wilcox.test", 
                     method.args = list(alternative = "two.sided")) 

ggsave(GSE43458, filename = 'GSE43458_deff.pdf', width = 8)


# ================================= sv =================================================
rm(list = ls())
library(tidyverse)
library("survival")
library("survminer")
library(limma)
library(timeROC)

res.cox <- readRDS(file = './res.cox_model.Rds')
sv <- c('GSE31210', 'GSE50081')
sv1 <- readRDS('GEO/matrix/GSE31210_exprs_metadata.Rds')
sv2 <- readRDS('GEO/matrix/GSE50081_exprs_metadata.Rds')

metadata1 <- read.csv('GEO/matrix/sv1.csv')
metadata2 <- read.csv('GEO/matrix/sv2.csv')

days_metadata1 <- str_extract(metadata1$characteristics_ch1.15,
                              pattern = '[0-9]+')
case_metadata1 <- str_extract(metadata1$characteristics_ch1.14, 
                              pattern = '(?<=death: ).*')
exclud_metadata1 <-  str_extract(metadata1$characteristics_ch1.17, 
                                 pattern = '(?<=exclude for prognosis analysis due to incomplete resection or adjuvant therapy: ).*')
se1 <- sv1$se[,metadata1$X]
se1 <- log2(1 + se1)
# saveRDS(se1, file = 'GEO1_log2.Rds')
sv1d <- data.frame(time = days_metadata1,
                   cas = case_metadata1,
                   exclud = exclud_metadata1,
                   GSM = metadata1$X)
f1 <- sv1d$exclud == 'none'
se1 <- se1[,f1]
sv1d <- sv1d[f1,]
sv1d$status <- 0
sv1d$status[sv1d$cas == 'dead'] <- 1
sv1_data <- sv1d[,c(1,4,5)]

days_metadata2 <- floor(as.numeric(str_extract(metadata2$characteristics_ch1.8,
                                               pattern = '(?<=survival time: ).+')
) * 365)
case_metadata2 <- str_extract(metadata2$characteristics_ch1.9, 
                              pattern = '(?<=status: ).*')
type_metadata2 <- metadata2$characteristics_ch1.1
se2 <- sv2$se[,metadata2$X]
# saveRDS(se2, file = 'GEO2_log2.Rds')
sv2d <- data.frame(time = days_metadata2,
                   cas = case_metadata2,
                   type = type_metadata2,
                   GSM = metadata2$X)

f1 <- sv2d$type == 'histology: adenocarcinoma'
sv2d <- sv2d[f1,]
se2 <- se2[,f1]
sv2d$status <- 0
sv2d$status[sv2d$cas == 'dead'] <- 1
sv2_data <- sv2d[,c(1,4,5)]


se1 <- as.data.frame(se1)
se2 <- as.data.frame(se2)

gn <- names(res.cox$coefficients)[names(res.cox$coefficients) %in% rownames(se1)]
names(res.cox$coefficients)[!names(res.cox$coefficients) %in% rownames(se1)]
# =============================== ste1 =========================================
se1_data <- as.data.frame(apply(t(se1[gn,]), 2, scale))
se1_data$C20orf197 <- 0
rownames(se1_data) <- sv1_data[,2]
se1p <- predict(res.cox, se1_data)
# saveRDS(se1p, file = 'HR_se1p.Rds')
# ===================== PCA ========================
com1 <- prcomp(se1_data[,1:16], center = TRUE,scale. = TRUE)
pc <- com1$x[,1:2]
pc <- as.data.frame(pc)
pc$group <- se1p
pc$level <- 'high'
pc$level[pc$group < 0] <- 'low'
# pdf('PCA_GEO1.pdf')
# ggplot(pc, aes(x = PC1, y = PC2, color = level)) + 
#   geom_point() + scale_color_manual(
#     values = c(
#       'high' = 'red',
#       'low' = 'blue')
#   ) + theme_bw()
# dev.off()


status <- sv1_data[,3]
time <- as.numeric(sv1_data[,1])
ROC.DSST <- timeROC(
  T=time,
  delta = status,
  marker = se1p,
  cause=1,
  weighting="marginal",
  times = c(365, 365 *3, 365*5),
  iid=TRUE,
  ROC = T)

new_dat <- sv1_data
new_dat <- as.data.frame(sv1_data)
new_dat$fp <- se1p
result <- with(new_dat, ROC.DSST)
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(365, 365 *3, 365*5)),
                            each = nrow(result$TP)))

# p1 <- ggplot() + 
#     geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
#     scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", 
#                                               "#66C2A5"),
#                        labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
#                                        format(round(result$AUC,2),nsmall = 2)))+
#     geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
#     theme_bw()+
#     theme(panel.grid = element_blank(),
#           legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
#           legend.position = c(0.765,0.125))+
#     scale_x_continuous(expand = c(0.005,0.005))+
#     scale_y_continuous(expand = c(0.005,0.005))+
#     labs(x = "1 - Specificity",
#          y = "Sensitivity")+
#     coord_fixed()
# p1
# ggsave(p1, filename = 'GEO_set1.pdf',width = 7,height = 7)


names(se1p) <- sv1_data[,2]
mx <- t(se1_data[,1:16])
od <- order(se1p, decreasing = F)
mx <- mx[,od]
ROCR.simple0 <- se1p[od]
vvv <- rep('high', times = length(se1p))
vvv[ROCR.simple0 < 0] <- 'low'
colnames(mx) <- NULL
# an <- HeatmapAnnotation(group = vvv,
#                         col = list(
#                             group = c('low' = "blue",
#                                       "high" = 'red')
#                         ))
# pdf('GEO1_heatmap.pdf', width = 7,height = 4)
# Heatmap(mx, cluster_columns = F, top_annotation = an)
# dev.off()

ddd <- data.frame(HR = se1p,
                  time = sv1_data[,'time'],
                  status = sv1_data[,'status'])
ddd$time <- as.numeric(ddd$time)
ddd <- ddd[order(ddd$HR, decreasing = F),]
ddd$rank <- 1:dim(ddd)[1]
ddd$group <- rep('high', times = length(sv1_data))
ddd$group[ddd$HR < 0] <- 'low'

# ppp <- ggplot(ddd, aes(x = rank, y = HR, color = group)) + 
#     geom_point() + theme_bw() + scale_color_manual(
#         values = c('red','blue')
#     )
# ppp
# ggsave(ppp, filename = 'GEO1_rank_HR.pdf', width = 7,height = 7)

# ddd$status <- factor(ddd$status, levels = c(0,1))
# ddd$stat <- 'Alive'
# ddd$stat[ddd$status == 1] <- 'Dead'
# ppp <- ggplot(ddd, aes(x = rank, y = time, color = stat)) + 
#     geom_point() + theme_bw() + scale_color_manual(
#         values = c('blue','red')
#     )
# ppp
# ggsave(ppp, filename = 'GEO1_rank_Dead.pdf', width = 7, height = 7)

ddd$status <- as.numeric(as.character(ddd$status))

fit <- survfit(Surv(time, status) ~ group, data = ddd)
# pdf('GEO1_sv1.pdf')
# ggsurvplot(fit,
#            pval = TRUE, conf.int = F,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            # linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("red", "blue"),
#            xlim = c(0, 3000)
# )
# ggsurvplot(fit,
#            conf.int = F, # 增加置信区间
#            fun = "cumhaz",
#            palette = c("red", "blue"),
#            xlim = c(0, 3000)) # 绘制累计风险曲线
# dev.off()

group <- rep('high', times = length(se1p))
group[se1p < 0] <- 'low'
mm <- model.matrix( ~ 0 + group)
fit <- lmFit(se1, mm)
head(coef(fit))

contr <- makeContrasts(grouphigh - grouplow, 
                       levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_HL <- topTable(tmp, sort.by = "P", n = Inf)
top.table_HL$`-logP` <- - log10(top.table_HL$P.Value)
top.table_HL$group <- 'Sig-UP'
top.table_HL$group[top.table_HL$logFC < 0.5 & 
                     top.table_HL$logFC > -0.5 |
                     top.table_HL$adj.P.Val > 0.05] <- 'No-Sig'
top.table_HL$group[top.table_HL$logFC < -0.5 & 
                     top.table_HL$adj.P.Val < 0.05] <- 'Sig-DN'

# library(clusterProfiler)
# set <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt')
# v <- top.table_HL$logFC
# names(v) <- rownames(top.table_HL)
# v <- sort(v, decreasing = T)
# gseaHL <- GSEA(v, TERM2GENE = set)
# gseaHLd <- gseaHL@result
# gseaHLd <- gseaHLd[order(gseaHLd$NES, decreasing = F),]
# gseaHLd$ID <- factor(gseaHLd$ID, levels = gseaHLd$ID)
# gseaHLd$color <- 'orange'
# gseaHLd$color[gseaHLd$NES < 0] <- 'green'
# 
library(viridis)
gseaHLd$absNES <- abs(gseaHLd$NES)
# write.csv(gseaHLd,file = 'GSEA_GEO1.csv')
gseaHLd <- read.csv('GSEA_GEO1.csv')
gseaHLd$ID <- factor(gseaHLd$ID, levels = gseaHLd$ID)
gseap <- ggplot(gseaHLd, aes(x = ID, y = NES,
                             fill = absNES)) +
  geom_bar(position = 'dodge', stat = 'identity', 
           color = 'black') +
  coord_flip() + theme_bw()  + 
  scale_fill_viridis(option="magma",
                     direction = -1,
                     begin = 0.35,
                     end = 0.9)
gseap
ggsave(gseap, filename = 'GEO1_gseapviridis.pdf', width = 10, height = 7)

# top.table_HL$group <- factor(top.table_HL$group, levels = c('Sig-DN',
#                                                             'Sig-UP',
#                                                             'No-Sig'))
# vo <- ggplot(top.table_HL, aes(x = logFC, y = `-logP`, color = group, 
#                                group = group)) + 
#     geom_point() +
#     theme_bw() + scale_color_jco()
# vo
# ggsave(vo, filename = 'GEO1_HL_volcano.pdf',width = 7, height = 7)

upg <- rownames(top.table_HL)[top.table_HL$group == 'Sig-UP']
dng <- rownames(top.table_HL)[top.table_HL$group == 'Sig-DN']

# listd <- list()
# setGObp <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt')
# upgebp <- enricher(upg,TERM2GENE = setGObp)
# dngebp <- enricher(dng,TERM2GENE = setGObp)
# listd[['bp']] <- list(
#     upgebp,
#     dngebp)
# 
# setGOcc <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.cc.v7.4.symbols.gmt')
# upgecc <- enricher(upg,TERM2GENE = setGOcc)
# dngecc <- enricher(dng,TERM2GENE = setGOcc)
# listd[['cc']] <- list(
#     upgecc,
#     dngecc)
# 
# setGOmf <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.mf.v7.4.symbols.gmt')
# upgemf <- enricher(upg,TERM2GENE = setGOmf)
# dngemf <- enricher(dng,TERM2GENE = setGOmf)
# listd[['mf']] <- list(
#     upgemf,
#     dngemf)
# 
# setKEGG <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c2.cp.kegg.v7.4.symbols.gmt')
# upgesetKEGG <- enricher(upg,TERM2GENE = setKEGG)
# dngesetKEGG <- enricher(dng,TERM2GENE = setKEGG)
# listd[['setKEGG']] <- list(
#     upgesetKEGG,
#     dngesetKEGG)
# saveRDS(listd, file = 'output/GEO1_HL_enrich.Rds')
# View(listd)

# m=0
# nl <- c('bp','cc','mf','kegg')
# for(i in names(listd)){
#     m = m + 1
#     fn1 <- paste0(nl[m],'up.pdf')
#     fn2 <- paste0(nl[m],'dn.pdf')
#     
#     f1 <- dotplot(listd[[i]][[1]])
#     f2 <- dotplot(listd[[i]][[2]])
#     
#     ggsave(f1, filename = fn1, width = 13, height = 7)
#     ggsave(f2, filename = fn2, width = 13, height = 7)
# }
# ======================= 临床森林图 =====================
st_sv1 <- read.csv('GEO/matrix/st_sv1.csv')

st_sv1 <- st_sv1[st_sv1$X %in% names(se1p),]
rownames(st_sv1) <- st_sv1$X
st_sv1 <- st_sv1[names(se1p),]
st_sv1$HR <- se1p
st_sv1$level <- 'High'
st_sv1$level[st_sv1$HR < 0] <- 'low'
st_sv1$age <- ' < 60'
st_sv1$age[st_sv1$age..years..ch1 >= 60] <- ' >= 60'


t1 <- table(st_sv1$gene.alteration.status.ch1,st_sv1$level)
chisq.test(t1)

t2 <- table(st_sv1$smoking.status.ch1,st_sv1$level)
chisq.test(t2)

t3 <- table(st_sv1$gender.ch1,st_sv1$level)
chisq.test(t3)

t4 <- table(st_sv1$myc.ch1,st_sv1$level)
chisq.test(t4)

# t5 <- table(st_sv1$pathological.stage.ch1,st_sv1$level)
# chisq.test(t5)

t6 <- table(st_sv1$pstage.iorii.ch1,st_sv1$level)
chisq.test(t6)

t7 <- table(st_sv1$relapse.ch1,st_sv1$level)
chisq.test(t7)

t8 <- table(st_sv1$age,st_sv1$level)
chisq.test(t8)

View(st_sv1)
colnames(st_sv1)
st_sv10 <- st_sv1[,c(4,6,2,8,14,16,17)]

# st_sv10$pstage.iorii.ch1 <- ifelse(st_sv10$pstage.iorii.ch1 == 'II', 2, 1)
st_sv10$pstage.iorii.ch1 <- factor(st_sv10$pstage.iorii.ch1,levels = c("I",
                                                                       "II"))
colnames(st_sv10)[c(1,2)] <- c('days','statu') 
st_sv10$statu <- ifelse(st_sv10$statu == 'alive',0,1)
# st_sv10$myc.ch1 <- ifelse(st_sv10$myc.ch1 == 'Low',0,1)
# st_sv10$smoking.status.ch1 <- ifelse(
#   st_sv10$smoking.status.ch1 == 'Ever-smoker',1,0)
st_sv10$smoking.status.ch1 <- factor(st_sv10$smoking.status.ch1,
                                     levels = c('Never-smoker',
                                                'Ever-smoker'))
# st_sv10$gene.alteration.status.ch1 <- factor(st_sv10$gene.alteration.status.ch1,
#                                              levels = c("EGFR/KRAS/ALK -",
#                                                         "ALK-fusion +",
#                                                         "KRAS mutation +",
#                                                         "EGFR mutation +"))
colnames(st_sv10)
glc <- list()
cl1 <- list()
for(i in 3:7){
  coxt <- st_sv10[,c(1,2,i)]
  res.cox <- coxph(Surv(days, statu) ~ ., data =  coxt)
  cox.zph(res.cox)$table
  res.cox <- summary(res.cox)
  glc[[colnames(st_sv10)[i]]] <-  as.data.frame(unlist(res.cox$coefficients))
  cl1[[colnames(st_sv10)[i]]] <- res.cox$concordance
}
glc <- purrr::reduce(glc, rbind)
cl1 <- purrr::reduce(cl1,rbind)
rownames(cl1) <- colnames(st_sv10)[3:7]

cl1 <- as.data.frame(cl1)
cl1$names <- rownames(cl1)
colnames(cl1)
cl1 <- cl1[order(cl1$C),]
cl1$names <- factor(cl1$names,levels = c(cl1$names))
ggplot(cl1,aes(x = names, y = C, 
               color = names, group = names)) + 
  geom_segment(aes(x = names, y = 0, 
                   yend = C, xend = names),
               color = 'grey',size = 1) + 
  geom_point(size = 4)+
  ylim(0,0.75)+
  coord_flip() + 
  theme_bw() -> GEO1Cindex
write.csv(cl1, file = 'GEO1Cindex.csv')
ggsave(GEO1Cindex, width = 6,height = 3, 
       filename = 'GEO1Cindex.pdf')

glc$Upper <- glc$coef + 1.96 * glc$`se(coef)`
glc$Lower <- glc$coef - 1.96 * glc$`se(coef)`
glc$Varnames <- rownames(glc)
View(glc)
p <- ggplot(glc, aes(x = coef, y = Varnames, 
                     group = Varnames, col = `Pr(>|z|)`)) # 不同形状shape= Factor

forsv <- p + geom_point(size=3) +
  
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.3,size=2) +
  
  scale_x_continuous(limits= c(-3, 3.5), breaks= seq(-3, 3.5, 1)) +
  
  geom_vline(aes(xintercept = 0)) +
  
  xlab('Coef') + ylab(' ') + theme_bw() + scale_color_viridis_c()
ggsave(forsv, filename = 'GEO1_single_forest.pdf', width=7,height = 2.5)
write.csv(glc, file = 'single_GEO1.csv')

tmp00 <- st_sv10[,c(1,2,5,7)]
res.cox <- coxph(Surv(days, statu) ~ ., data =  tmp00)
cox.zph(res.cox)$table
res.cox <- summary(res.cox)
res.cox <- res.cox$coefficients
res.cox <- as.data.frame(res.cox)
res.cox$Lower <- res.cox$coef - 1.96 * res.cox$`se(coef)`
res.cox$Upper <- res.cox$coef + 1.96 * res.cox$`se(coef)`
res.cox$Varnames <- rownames(res.cox)


p <- ggplot(res.cox, aes(x = coef, y = Varnames, 
                         group = Varnames, col = `Pr(>|z|)`)) # 不同形状shape= Factor

forsv <- p + geom_point(size=3) +
  
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.3,size=2) +
  
  scale_x_continuous(limits= c(-4, 4), breaks= seq(-2, 3.5, 1)) +
  
  geom_vline(aes(xintercept = 0)) +
  
  xlab('Coef') + ylab(' ') + theme_bw() + scale_color_viridis_c()
ggsave(forsv, filename = 'GEO1_mulit_forest.pdf',width = 7,height = 1.5)
write.csv(res.cox, file = 'multi_GEO1.csv')

st_sv1$Stage <- paste0('Stage', st_sv1$pstage.iorii.ch1)
mycomparision <- list(c('StageII', 'StageI'))
Sp1 <- ggplot(st_sv1, aes(x = Stage, 
                          y = HR, fill = Stage))+ 
  geom_boxplot() + theme_bw() + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, 
                     label = 'p.signif')
Sp1
ggsave(Sp1, filename = 'Stage_GEO1.pdf',width = 3.5,height = 7)
colnames(st_sv1)
st_sv1$gene.alteration.status.ch1 <- factor(st_sv1$gene.alteration.status.ch1,
                                            levels = c("EGFR/KRAS/ALK -",
                                                       "ALK-fusion +",
                                                       "KRAS mutation +",
                                                       "EGFR mutation +"))
View(st_sv10)
mycomparision <- list(c("EGFR/KRAS/ALK -", 'ALK-fusion +'),
                      c("EGFR/KRAS/ALK -", 'KRAS mutation +'),
                      c("EGFR/KRAS/ALK -", 'EGFR mutation +'))
altp1 <- ggplot(st_sv1, aes(x = gene.alteration.status.ch1, 
                            y = HR, fill = gene.alteration.status.ch1))+ 
  geom_boxplot() + theme_bw()+
  theme(axis.text.x = 
          element_text(angle = 45, 
                       hjust = 0.5, 
                       vjust = 0.5)) + 
  ylim(-2,4.5) +
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, 
                     aes(label = 'p.signif')) 
ggsave(altp1, filename = 'altp1_GEO1.pdf',width = 7.2,height = 7)
# ========================= set2 ==============================================
se2_data <- as.data.frame(apply(t(se2[gn,]), 2, scale))
se2_data$C20orf197 <- 0
rownames(se2_data) <- colnames(se2)
res.cox <- readRDS(file = 'res.cox_model.Rds')
se2p <- predict(res.cox, se2_data)
# saveRDS(se2p, file = 'HR_GEO2.Rds')
# ===================== PCA ========================
com1 <- prcomp(se2_data[,1:16], center = TRUE,scale. = TRUE)
pc <- com1$x[,1:2]
pc <- as.data.frame(pc)
pc$group <- se2p
pc$level <- 'high'
pc$level[pc$group < 0] <- 'low'
# pdf('PCA_GEO2.pdf')
# ggplot(pc, aes(x = PC1, y = PC2, color = level)) + 
#   geom_point() + scale_color_manual(
#     values = c(
#       'high' = 'red',
#       'low' = 'blue')
#   ) + theme_bw()
# dev.off()


status <- sv2_data[,3]
time <- sv2_data[,1]
time <- as.numeric(time)
ROC.DSST <- timeROC(
  T=time,
  delta = status,
  marker = se2p,
  cause=1,
  weighting="marginal",
  times = c(365, 365 *3, 365*5),
  iid=TRUE,
  ROC = T)

new_dat <- sv2_data
new_dat <- as.data.frame(sv2_data)
new_dat$fp <- se2p
result <- with(new_dat, ROC.DSST)
dat = data.frame(fpr = as.numeric(result$FP),
                 tpr = as.numeric(result$TP),
                 time = rep(as.factor(c(365, 365 *3, 365*5)),
                            each = nrow(result$TP)))

# p2 <- ggplot() + 
#     geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
#     scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", 
#                                               "#66C2A5"),
#                        labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
#                                        format(round(result$AUC,2),nsmall = 2)))+
#     geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
#     theme_bw()+
#     theme(panel.grid = element_blank(),
#           legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
#           legend.position = c(0.765,0.125))+
#     scale_x_continuous(expand = c(0.005,0.005))+
#     scale_y_continuous(expand = c(0.005,0.005))+
#     labs(x = "1 - Specificity",
#          y = "Sensitivity")+
#     coord_fixed()
# p2
# ggsave(p2, filename = 'GEO_set2.pdf',width = 7,height = 7)

names(se2p) <- sv2_data[,2]
mx <- t(se2_data[,1:16])
od <- order(se2p, decreasing = F)
mx <- mx[,od]
ROCR.simple0 <- se2p[od]
vvv <- rep('high', times = length(se2p))
vvv[ROCR.simple0 < 0] <- 'low'
colnames(mx) <- NULL
# an <- HeatmapAnnotation(group = vvv,
#                         col = list(
#                             group = c('low' = "blue",
#                                       "high" = 'red')
#                         ))
# pdf('GEO2_heatmap.pdf', width = 7,height = 4)
# Heatmap(mx, cluster_columns = F, top_annotation = an)
# dev.off()

ddd <- data.frame(HR = se2p,
                  time = sv2_data[,'time'],
                  status = sv2_data[,'status'])
ddd$time <- as.numeric(ddd$time)
ddd <- ddd[order(ddd$HR, decreasing = F),]
ddd$rank <- 1:dim(ddd)[1]
ddd$group <- rep('high', times = length(sv2_data$time))
ddd$group[ddd$HR < 0] <- 'low'

# ppp <- ggplot(ddd, aes(x = rank, y = HR, color = group)) + 
#     geom_point() + theme_bw() + scale_color_manual(
#         values = c('red','blue')
#     )
# ppp
# ggsave(ppp, filename = 'GEO2_rank_HR.pdf', width = 7,height = 7)

# ddd$status <- factor(ddd$status, levels = c(0,1))
# ddd$stat <- 'Alive'
# ddd$stat[ddd$status == 1] <- 'Dead'
# ppp <- ggplot(ddd, aes(x = rank, y = time, color = stat)) + 
#     geom_point() + theme_bw() + scale_color_manual(
#         values = c('blue','red')
#     )
# ppp
# ggsave(ppp, filename = 'GEO2_rank_Dead.pdf', width = 7, height = 7)

# ddd$status <- as.numeric(as.character(ddd$status))
# 
# fit <- survfit(Surv(time, status) ~ group, data = ddd)
# pdf('GEO2_sv1.pdf')
# ggsurvplot(fit,
#            pval = TRUE, conf.int = F,
#            risk.table = TRUE, # Add risk table
#            risk.table.col = "strata", # Change risk table color by groups
#            # linetype = "strata", # Change line type by groups
#            surv.median.line = "hv", # Specify median survival
#            ggtheme = theme_bw(), # Change ggplot2 theme
#            palette = c("red", "blue"),
#            xlim = c(0, 3000)
# )
# ggsurvplot(fit,
#            conf.int = F, # 增加置信区间
#            fun = "cumhaz",
#            palette = c("red", "blue"),
#            xlim = c(0, 3000)) # 绘制累计风险曲线
# dev.off()


group <- rep('high', times = length(se2p))
group[se2p < 0] <- 'low'
mm <- model.matrix( ~ 0 + group)
fit <- lmFit(se2, mm)
head(coef(fit))

contr <- makeContrasts(grouphigh - grouplow, 
                       levels = colnames(coef(fit)))
tmp <- contrasts.fit(fit, contr)
tmp <- eBayes(tmp)
top.table_HL <- topTable(tmp, sort.by = "P", n = Inf)
top.table_HL$`-logP` <- - log10(top.table_HL$P.Value)
top.table_HL$group <- 'Sig-UP'
top.table_HL$group[top.table_HL$logFC < 0.5 & 
                     top.table_HL$logFC > -0.5 |
                     top.table_HL$adj.P.Val > 0.05] <- 'No-Sig'
top.table_HL$group[top.table_HL$logFC < -0.5 & 
                     top.table_HL$adj.P.Val < 0.05] <- 'Sig-DN'


# set <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt')
# v <- top.table_HL$logFC
# names(v) <- rownames(top.table_HL)
# v <- sort(v, decreasing = T)
# gseaHL <- GSEA(v, TERM2GENE = set)
# gseaHLd <- gseaHL@result
# gseaHLd <- gseaHLd[order(gseaHLd$NES, decreasing = F),]
# gseaHLd$ID <- factor(gseaHLd$ID, levels = gseaHLd$ID)
# gseaHLd$color <- 'orange'
# gseaHLd$color[gseaHLd$NES < 0] <- 'green'

gseaHLd$absNES <- abs(gseaHLd$NES)
# write.csv(gseaHLd,file = 'GEO2_GSEA.csv')
gseaHLd <- read.csv('GEO2_GSEA.csv')
gseaHLd$ID <- factor(gseaHLd$ID, levels = gseaHLd$ID)
gseap <- ggplot(gseaHLd, aes(x = ID, 
                             y = NES, 
                             fill = absNES)) +
  geom_bar(position = 'dodge', stat = 'identity', color = 'black') +
  coord_flip() + theme_bw()  + 
  scale_fill_viridis(option="magma",
                     direction = -1,
                     begin = 0.35,
                     end = 0.9)
gseap

ggsave(gseap, filename = 'GEO2_gseapviridis.pdf', width = 10, height = 7)

top.table_HL$group <- factor(top.table_HL$group, levels = c('Sig-DN',
                                                            'Sig-UP',
                                                            'No-Sig'))
vo <- ggplot(top.table_HL, aes(x = logFC, y = `-logP`, color = group, 
                               group = group)) + 
  geom_point() +
  theme_bw() + scale_color_jco()
vo
# ggsave(vo, filename = 'GEO2_HL_volcano.pdf',width = 7, height = 7)
# 
# upg <- rownames(top.table_HL)[top.table_HL$group == 'Sig-UP']
# dng <- rownames(top.table_HL)[top.table_HL$group == 'Sig-DN']
# 
# listd <- list()
# setGObp <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt')
# upgebp <- enricher(upg,TERM2GENE = setGObp)
# dngebp <- enricher(dng,TERM2GENE = setGObp)
# listd[['bp']] <- list(
#     upgebp,
#     dngebp)
# 
# setGOcc <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.cc.v7.4.symbols.gmt')
# upgecc <- enricher(upg,TERM2GENE = setGOcc)
# dngecc <- enricher(dng,TERM2GENE = setGOcc)
# listd[['cc']] <- list(
#     upgecc,
#     dngecc)
# 
# setGOmf <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c5.go.mf.v7.4.symbols.gmt')
# upgemf <- enricher(upg,TERM2GENE = setGOmf)
# dngemf <- enricher(dng,TERM2GENE = setGOmf)
# listd[['mf']] <- list(
#     upgemf,
#     dngemf)
# 
# setKEGG <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/c2.cp.kegg.v7.4.symbols.gmt')
# upgesetKEGG <- enricher(upg,TERM2GENE = setKEGG)
# dngesetKEGG <- enricher(dng,TERM2GENE = setKEGG)
# listd[['setKEGG']] <- list(
#     upgesetKEGG,
#     dngesetKEGG)
# saveRDS(listd, file = 'output/GEO2_HL_enrich.Rds')
# View(listd)
# 
# m=0
# nl <- c('bp','cc','mf','kegg')
# for(i in names(listd)){
#     m = m + 1
#     fn1 <- paste0(nl[m],'up.pdf')
#     fn2 <- paste0(nl[m],'dn.pdf')
#     
#     f1 <- dotplot(listd[[i]][[1]])
#     f2 <- dotplot(listd[[i]][[2]])
#     
#     ggsave(f1, filename = fn1, width = 13, height = 7)
#     ggsave(f2, filename = fn2, width = 13, height = 7)
# }
# ==================== 森林图 ====================================
st_sv2 <- read.csv('GEO/matrix/st_sv2.csv')
st_sv2 <- st_sv2[st_sv2$X %in% names(se2p),]
rownames(st_sv2) <- st_sv2$X
st_sv2 <- st_sv2[names(se2p),]
st_sv2$HR <- se2p
st_sv2$level <- 'High'
st_sv2$level[st_sv2$HR < 0] <- 'low'

st_sv2$age <- ' < 60'
st_sv2$age[st_sv2$age.ch1 >= 60] <- ' >= 60'

st_sv2$Stage <- '1'
st_sv2$Stage[st_sv2$Stage.ch1 %in% c('2B','2A')] <- '2'


tt1 <- table(st_sv2$n.stage.ch1,st_sv2$level)
chisq.test(tt1)

tt2 <- table(st_sv2$recurrence.ch1,st_sv2$level)
chisq.test(tt2)

tt3 <- table(st_sv2$Sex.ch1,st_sv2$level)
chisq.test(tt3)

tt4 <- table(st_sv2$smoking.ch1,st_sv2$level)
chisq.test(tt4)

# tt5 <- table(st_sv2$Stage.ch1,st_sv2$level)
# chisq.test(tt5)

tt6 <- table(st_sv2$t.stage.ch1,st_sv2$level)
chisq.test(tt6)

tt7 <- table(st_sv2$age,st_sv2$level)
chisq.test(tt7)

tt8 <- table(st_sv2$Stage,st_sv2$level)
chisq.test(tt8)

View(st_sv2)
colnames(st_sv2)
st_sv2$smoking <- "Ex-smoker/Current"
st_sv2$smoking[st_sv2$smoking.ch1 == 'Never'] <- 'Nevers'
st_sv2$smoking[st_sv2$smoking.ch1 == "Unable to determine"] <- NA
st_sv2$smoking <- factor(st_sv2$smoking, levels = c("Nevers",
                                                    "Ex-smoker/Current"))
unique(st_sv2$smoking)

st_sv20 <- st_sv2[,c(11,12,2,6,8,18,13,14,17)]
st_sv20$t.stage.ch1 <- as.factor(st_sv20$t.stage.ch1)
st_sv20$status.ch1 <- ifelse(st_sv20$status.ch1 == 'dead',1,0)
colnames(st_sv20)[c(1,2)] <- c('statu','years')

View(st_sv20)
glc <- list()
cl2 <- list()
for(i in 3:9){
  coxt <- st_sv20[,c(1,2,i)]
  res.cox <- coxph(Surv(years, statu) ~ ., data =  coxt)
  cox.zph(res.cox)$table
  res.cox <- summary(res.cox)
  glc[[colnames(st_sv20)[i]]] <-  as.data.frame(unlist(res.cox$coefficients))
  cl2[[colnames(st_sv20)[i]]] <- res.cox$concordance
  
}
glc <- purrr::reduce(glc, rbind)
cl2 <- purrr::reduce(cl2,rbind)
rownames(cl2) <- colnames(st_sv20)[3:9]

cl2 <- as.data.frame(cl2)
cl2$names <- rownames(cl2)
colnames(cl2)
cl2 <- cl2[order(cl2$C),]
cl2$names <- factor(cl2$names,levels = c(cl2$names))
ggplot(cl2,aes(x = names, y = C, 
               color = names, group = names)) + 
  geom_segment(aes(x = names, y = 0, 
                   yend = C, xend = names),
               color = 'grey',size = 1) + 
  geom_point(size = 4)+
  ylim(0,0.75)+
  coord_flip() + 
  theme_bw() -> GEO2Cindex
write.csv(cl2, file = 'GEO2Cindex.csv')
ggsave(GEO2Cindex, width = 4.5,height = 3, 
       filename = 'GEO2Cindex.pdf')



rownames(cl2) <- colnames(st_sv20)[3:9]
min(cl2[,2])
pnorm(0.54, mean = 0.64, sd = 0.0374)
pnorm(0.57, mean = 0.64, sd = 0.0374)
pnorm(0.54, mean = 0.64, sd = 0.0374)
pnorm(0.57, mean = 0.64, sd = 0.0374)
pnorm(0.60, mean = 0.64, sd = 0.032)
pnorm(0.59, mean = 0.64, sd = 0.032)

geo2clist <- list()
colnames(st_sv20)
for(i in c(3,4,5,7,9))
{
  compareC(st_sv20$years, st_sv20$statu, 
           st_sv20[,i], st_sv20$HR) -> ccc
  geo2clist[[i]] <- ccc$pval
}



colnames(st_sv20)
status <- sv2_data[,3]
time <- sv2_data[,1]
time <- as.numeric(time)
mk2 <- as.numeric(st_sv20[,7])
ROC.DSST_stage <- timeROC(
  T=time,
  delta = status,
  marker = mk2,
  cause=1,
  weighting="marginal",
  times = c(365, 365 *3, 365*5),
  iid=TRUE,
  ROC = T)


dat = data.frame(fpr = as.numeric(ROC.DSST_stage$FP),
                 tpr = as.numeric(ROC.DSST_stage$TP),
                 time = rep(as.factor(c(365, 365 *3, 365*5)),
                            each = nrow(ROC.DSST_stage$TP)))

p2 <- ggplot() + 
  geom_line(data = dat,aes(x = fpr, y = tpr,color = time),size = 1) + 
  scale_color_manual(name = NULL,values = c("#92C5DE", "#F4A582", 
                                            "#66C2A5"),
                     labels = paste0("AUC of ",c(1,3,5),"-y survival: ",
                                     format(round(result$AUC,2),nsmall = 2)))+
  geom_line(aes(x=c(0,1),y=c(0,1)),color = "grey")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        legend.background = element_rect(linetype = 1, size = 0.2, colour = "black"),
        legend.position = c(0.765,0.125))+
  scale_x_continuous(expand = c(0.005,0.005))+
  scale_y_continuous(expand = c(0.005,0.005))+
  labs(x = "1 - Specificity",
       y = "Sensitivity")+
  coord_fixed()
p2





glc$Upper <- glc$coef + 1.96 * glc$`se(coef)`
glc$Lower <- glc$coef - 1.96 * glc$`se(coef)`
glc$Varnames <- rownames(glc)
View(glc)
glc <- na.omit(glc)
p <- ggplot(glc, aes(x = coef, y = Varnames, 
                     group = Varnames, col = `Pr(>|z|)`)) # 不同形状shape= Factor

forsv <- p + geom_point(size=3) +
  
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.3,size=2) +
  
  scale_x_continuous(limits= c(-2.5, 6), breaks= seq(-2.5, 6, 1)) +
  
  geom_vline(aes(xintercept = 0)) +
  
  xlab('Coef') + ylab(' ') + theme_bw() + scale_color_viridis_c()
forsv
ggsave(forsv, filename = 'GEO2_single_forest.pdf', 
       width=7,height = 3)
write.csv(glc, file = 'single_GEO2.csv')

tmp00 <- st_sv20[,c(1,2,4,7,8,9)]
tmp00$Stage <- as.numeric(tmp00$Stage)
res.cox <- coxph(Surv(years, statu) ~ ., data =  tmp00)
cox.zph(res.cox)$table
res.cox <- summary(res.cox)
res.cox <- res.cox$coefficients
res.cox <- as.data.frame(res.cox)
res.cox$Lower <- res.cox$coef - 1.96 * res.cox$`se(coef)`
res.cox$Upper <- res.cox$coef + 1.96 * res.cox$`se(coef)`
res.cox$Varnames <- rownames(res.cox)
res.cox <- na.omit(res.cox)

p <- ggplot(res.cox, aes(x = coef, y = Varnames, 
                         group = Varnames, col = `Pr(>|z|)`)) # 不同形状shape= Factor

forsv <- p + geom_point(size=3) +
  
  geom_errorbarh(aes(xmax =Upper, xmin = Lower), height = 0.3,size=2) +
  
  scale_x_continuous(limits= c(-4, 5.5), breaks= seq(-4, 5.5, 1)) +
  
  geom_vline(aes(xintercept = 0)) +
  
  xlab('Coef') + ylab(' ') + theme_bw() + scale_color_viridis_c()
forsv
ggsave(forsv, filename = 'GEO2_mulit_forest.pdf',width = 7,height = 2)
write.csv(res.cox, file = 'multi_GEO2.csv')



colnames(st_sv2)
st_sv2$Stage <- paste0('Stage', st_sv2$Stage)
mycomparision <- list(c('Stage2', 'Stage1'))
Stagep2 <- ggplot(st_sv2, aes(x = Stage, 
                              y = HR, fill = Stage))+ 
  geom_boxplot() + theme_bw()+
  theme(axis.text.x = 
          element_text(angle = 45, 
                       hjust = 0.5, 
                       vjust = 0.5)) + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, label = 'p.signif')
ggsave(Stagep2, filename = 'Stage_GEO2.pdf',width = 3.5,height = 7)

st_sv2$t.stage.ch1 <- paste0('T_stage', st_sv2$t.stage.ch1)
mycomparision <- list(c('T_stage1', 'T_stage2'),
                      c('T_stage1', 'T_stage3'))
tp2 <- ggplot(st_sv2, aes(x = t.stage.ch1, 
                          y = HR, fill = t.stage.ch1))+ 
  geom_boxplot() + theme_bw()+
  theme(axis.text.x = 
          element_text(angle = 45, 
                       hjust = 0.5, 
                       vjust = 0.5)) + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, label = 'p.signif')
ggsave(tp2, filename = 'TStage_GEO2.pdf',width = 4.6,height = 7)


st_sv2$n.stage.ch1 <- paste0('N_stage', st_sv2$n.stage.ch1)
mycomparision <- list(c('N_stage1','N_stage0'))
np2 <- ggplot(st_sv2, aes(x = n.stage.ch1, 
                          y = HR, fill = n.stage.ch1))+ 
  geom_boxplot() + theme_bw()+
  theme(axis.text.x = 
          element_text(angle = 45, 
                       hjust = 0.5, 
                       vjust = 0.5)) + 
  stat_compare_means(method = 'wilcox.test', 
                     comparisons = mycomparision, label = 'p.signif')
ggsave(np2, filename = 'NStage_GEO2.pdf',width = 3.5,height = 7)
res.cox <- readRDS('res.cox_model.Rds')
names(res.cox$coefficients)

#  ====== ELMER karyoploteR=====
# library(karyoploteR)
# library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TFBSTools)
# library(rtfbs)
library(motifmatchr)
library(ELMER)
matrixl <- list()
for(i in c('MA0106.1.pfm','MA0106.2.pfm','MA0106.3.pfm')){
  purrr::reduce(
    lapply(
      strsplit(readLines(i)[2:5],'(\\s)'),
      function(x){
        as.integer(x[x!=''])
      } 
    ),rbind) -> m1
  rownames(m1) <- c('A','C','G','T')
  matrixl[[i]] <- m1
}
gl <- list()
for(gene in names(res.cox$coefficients)){
  tssAnnot <- ELMER::getTSS(genome = "hg38",
                            TSS=list(upstream=1000, downstream=1000))
  tssAnnot <- tssAnnot[tssAnnot$external_gene_name == gene]
  
  mtfl <- list()
  for(motif0 in 1:3){
    pfm <- PFMatrix(ID=names(matrixl)[motif0], name="tP53", 
                    profileMatrix = matrixl[[motif0]]
    )
    
    motif_tot <- matchMotifs(pfm,tssAnnot, genome = "hg38")
    motif_tot_m <- motifMatches(motif_tot)
    (motif_tot_m_matrix <- as.matrix(motif_tot_m))
    mtfl[[names(matrixl)[motif0]]] <- motif_tot_m_matrix
    
  }
  gl[[gene]] <- as.data.frame(mtfl)
}
gl0 <- lapply(gl, function(x) do.call(rbind,x))
gl0l <- lapply(gl0, sum)
cn <- lapply(gl0, function(x) dim(x)[1] * dim(x)[2] )
which(unlist(gl0l))
unlist(gl0l) / unlist(cn)
unlist(gl0l)[unlist(gl0l)!=0]
do.call(rbind,lapply(gl, colSums)) -> heatmap0
heatmap0[rowSums(heatmap0) != 0,]  -> heatmap1
colnames(heatmap1) <- gsub(colnames(heatmap1),pattern = '\\.pfm',replacement = '')
library(ComplexHeatmap)
library(circlize)

heatmap1[heatmap1==0] <- NA
heatmap1[!is.na(heatmap1)] <- 1
pdf('p53_motif.pdf')
Heatmap(heatmap1,
        col = colorRamp2(
          c(0, 1.5), 
          c("white","red")
        ),
        na_col = '#6495ED',
        cluster_rows = F,
        cluster_columns = F,
        rect_gp = gpar(col = "black", lwd = 2))

dev.off()





library(hdf5r)
library(dplyr)
library(Seurat)
library(patchwork)
library(harmony)
library(tidyverse)
library(glmGamPoi)
library(Matrix)

file.h5 <- H5File$new('GSE154989_mmLungPlate_fQC_dSp_rawCount.h5', mode = "r+")
ii <- file.h5[['i']]
jj <- file.h5[['j']]
vv <- file.h5[['v']]

A <- sparseMatrix(file.h5[['i']]$read(), 
                  file.h5[['j']]$read(), 
                  x = file.h5[['v']]$read())
dim(A)
smpTable_annot <- read.csv('GSE154989_mmLungPlate_fQC_dZ_annot_smpTable.csv.gz') 
smpTable_QCstat <- read.csv('GSE154989_mmLungPlate_fQC_dZ_QCstat_smpTable.csv.gz') 
geneTable <- read.csv('GSE154989_mmLungPlate_fQC_geneTable.csv.gz') 
smpTable <- read.csv('GSE154989_mmLungPlate_fQC_smpTable.csv.gz') 

colnames(A) <- smpTable$sampleID
rownames(A) <- geneTable$geneID

pbmc <- CreateSeuratObject(counts = A)
pbmc@meta.data <- cbind(pbmc@meta.data,smpTable, smpTable_annot, smpTable_QCstat)
sp <- strsplit(pbmc@meta.data$sampleID,split = '_')
sp <- purrr::reduce(sp, rbind)
colnames(sp) <- c('kp','w','dn','animal','tumor','plate','cell')
sp <- as.data.frame(sp)
pbmc@meta.data <- cbind(pbmc@meta.data,sp)

dim(table(smpTable$plateID, smpTable$mouseID))

# VlnPlot(pbmc, 
#         features = c("expCount", "expMeanHK", 
#                      "expPropHK","expMeanMt","expPropMt",
#                      "PCT_RIBOSOMAL_BASES","PCT_INTERGENIC_BASES"), 
#         ncol = 5)
# pbmc <- NormalizeData(pbmc, 
#                      normalization.method = "LogNormalize", 
#                      scale.factor = 10000)
# pbmc <- ScaleData(pbmc, features = rownames(pbmc))

pbmc <- SCTransform(pbmc, method = "glmGamPoi")
pbmc <- FindVariableFeatures(pbmc, 
                             selection.method = "vst", 
                             nfeatures = 2000)
pbmc <- RunPCA(object = pbmc)
ElbowPlot(pbmc,ndims = 50)
pbmc <- FindNeighbors(pbmc, dims = 1:40)

library(clustree)
sce <- FindClusters(
  object = pbmc,
  resolution = c(seq(.3,1.6,.2))
)
clustree(sce@meta.data, prefix = "SCT_snn_res.")
seq(.3,1.6,.2)
pbmc <- FindClusters(pbmc, resolution = 0.5)
heatmap(table(pbmc$seurat_clusters,pbmc$clusterK12), scale = 'column')
pbmc <- RunICA(pbmc)
pbmc <- RunTSNE(pbmc, dims = 1:40)
pbmc <- RunUMAP(pbmc, dims = 1:40)

# saveRDS(pbmc.markers, file = 'HPSC.markers.Rds')
# saveRDS(pbmc, file = 'HPSC_obj.Rds')

# ================================ analysis ============================
library(Seurat)
library(ComplexHeatmap)
library(tidyverse)

pbmc <- readRDS('HPSC_obj.Rds')

# pbmc_marker <- Seurat::FindAllMarkers(pbmc)
# write.csv(pbmc_marker, file = 'HPSC_marker.csv')
# DimPlot(object = pbmc, reduction = "pca")
# DimPlot(object = pbmc, reduction = "ica")
# pdf('dim_reduce.pdf')
# DimPlot(object = pbmc, reduction = "tsne", label = T)
# DimPlot(object = pbmc, reduction = "umap", label = T)
# dev.off()

# DimPlot(object = pbmc, reduction = "umap", label = F,group.by = 'mouseID') + 
#   NoLegend()
# DimPlot(object = pbmc, reduction = "umap", label = T,group.by = 'kp')
# DimPlot(object = pbmc, reduction = "umap", label = T,group.by = 'w')
# DimPlot(object = pbmc, reduction = "umap", label = T,group.by = 'tumor')
# DimPlot(object = pbmc, reduction = "umap", label = T,group.by = 'tumor')


# FeaturePlot(pbmc, features = c('tdTomato-tdTomato'))
# pdf('gene.big.pdf')
# FeaturePlot(pbmc, features = c('ENSMUSG00000033161.10-Atp1a1'), label = T)
# FeaturePlot(pbmc, features = c('ENSMUSG00000016496.7-Cd274'), label = T)
# FeaturePlot(pbmc, features = c('ENSMUSG00000028382.15-Ptbp3'), label = T)
# FeaturePlot(pbmc, features = c('ENSMUSG00000071552.4-Tigit'), label = T)
# FeaturePlot(pbmc, features = c('ENSMUSG00000026104.14-Stat1'), label = T)
# FeaturePlot(pbmc, features = c('ENSMUSG00000004040.16-Stat3'), label = T)
# 
# # FeaturePlot(pbmc, features = c('ENSMUSG00000024661.6-Fth1'), label = T)
# FeaturePlot(pbmc, features = c('ENSMUSG00000039153.16-Runx2'), label = T)
# FeaturePlot(pbmc, features = c('ENSMUSG00000001496.15-Nkx2-1'), label = T)
# dev.off()

# pdf('violin.big.pdf',height = 15)
# VlnPlot(pbmc, 
#         features = c(
#           "ENSMUSG00000033161.10-Atp1a1",
#           "ENSMUSG00000026104.14-Stat1",
#           "ENSMUSG00000004040.16-Stat3",
#           'ENSMUSG00000039153.16-Runx2',
#           'ENSMUSG00000001496.15-Nkx2-1',
#           'ENSMUSG00000016496.7-Cd274',
#           'ENSMUSG00000071552.4-Tigit',
#           'ENSMUSG00000028382.15-Ptbp3'),
#         ncol = 1)
# dev.off()


pdf('percent2.pdf')
# tpd <- data.frame(
#   seurat_clusters = pbmc$seurat_clusters,
#   tumor = pbmc$tumor)
# ggplot(tpd, aes(x = tumor, fill = seurat_clusters)) + 
#   geom_bar(position = 'fill', color = 'black') + theme_bw()
# 
# htpd <- table(tpd)
# colnames(htpd)
# htpd <- apply(htpd, 1, function(x) x / sum(x))
# Heatmap(htpd,
#         col = circlize::colorRamp2(c(0,1), c("#6495ED","red")))

# twd <- data.frame(
#   week = pbmc$w,
#   tumor = pbmc$tumor)
# ggplot(twd, aes(x = tumor, fill = week)) + 
#   geom_bar(position = 'fill', color = 'black') + theme_bw()
# 
# twdd <- table(twd)
# twdd <- apply(twdd, 1, function(x) x / sum(x))
# Heatmap(twdd,col = circlize::colorRamp2(c(0,1), c("#6495ED","red")))

# seurattwd <- data.frame(
#   seurat_clusters = pbmc$seurat_clusters,
#   week = pbmc$w)
# seurattwd$week <- factor(seurattwd$week, levels = c('0w','2w','4w','12w','18w','20w','30w'))
# ggplot(seurattwd, aes(x = week, fill = seurat_clusters)) + 
#   geom_bar(position = 'fill', color = 'black') + theme_bw()

seurattwdd <- table(seurattwd)
seurattwdd <- apply(seurattwdd, 1, function(x) x / sum(x))
Heatmap(seurattwdd,cluster_rows = F,
        col = circlize::colorRamp2(c(0,1), c("#6495ED","red")))
dev.off()

library(tidyverse)
# markers <- readRDS('HPSC.markers.Rds')
# markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> markers0
# markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> markers2
# pdf('violin.cluster.pdf',height = 60)
# VlnPlot(pbmc, 
#         features = markers2$gene,
#         ncol = 1)
# dev.off()

pdf('feature.pdf', width = 10, height = 10)
DoHeatmap(pbmc, features = markers0$gene, size = 5, angle = -50, hjust=0.8) + NoLegend()
dev.off()

library(ggrepel)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 9,ident.2 = 5, logfc.threshold = 0.25, 
                                only.pos = F)
cluster1.markers$`-log10p` <- -log10(cluster1.markers$p_val_adj)

st <- cluster1.markers[abs(cluster1.markers$avg_log2FC) > 5 &
                         cluster1.markers$`-log10p` > 30,]
pdf('p9vs5.pdf', width = 10,height = 10)
ggplot(cluster1.markers, aes(x = avg_log2FC, y = `-log10p`)) + 
  geom_point() + geom_text_repel(st,mapping = aes(x = avg_log2FC, 
                                                  y = `-log10p`, 
                                                  label = rownames(st),
                                                  color = 'red')) + 
  theme_bw()
dev.off()

cluster75.markers <- FindMarkers(pbmc, ident.1 = 7,ident.2 = 5, 
                                 logfc.threshold = 0.25, 
                                 only.pos = F)
cluster75.markers$`-log10p` <- -log10(cluster75.markers$p_val_adj)
st <- cluster75.markers[abs(cluster75.markers$avg_log2FC) > 5 &
                          cluster75.markers$`-log10p` > 30,]
pdf('p7vs5.pdf', width = 10,height = 10)
ggplot(cluster75.markers, aes(x = avg_log2FC, y = `-log10p`)) + 
  geom_point() + geom_text_repel(st,mapping = aes(x = avg_log2FC, 
                                                  y = `-log10p`, 
                                                  label = rownames(st),
                                                  color = 'red')) + 
  theme_bw()
dev.off()

cluster97.markers <- FindMarkers(pbmc, ident.1 = 9,ident.2 = 7, 
                                 logfc.threshold = 0.25, 
                                 only.pos = F)
cluster97.markers$`-log10p` <- -log10(cluster97.markers$p_val_adj)
st <- cluster97.markers[abs(cluster97.markers$avg_log2FC) > 1 &
                          cluster97.markers$`-log10p` > 30,]
pdf('p9vs7.pdf', width = 15,height = 15)
ggplot(cluster97.markers, aes(x = avg_log2FC, y = `-log10p`)) + 
  geom_point() + geom_text_repel(st,mapping = aes(x = avg_log2FC, 
                                                  y = `-log10p`, 
                                                  label = rownames(st),
                                                  color = 'red')) + 
  theme_bw()
dev.off()

cluster91.markers <- FindMarkers(pbmc, ident.1 = 9,ident.2 = 1, 
                                 logfc.threshold = 0.25, 
                                 only.pos = F)
cluster91.markers$`-log10p` <- -log10(cluster91.markers$p_val_adj)
st <- cluster91.markers[abs(cluster91.markers$avg_log2FC) > 5 &
                          cluster91.markers$`-log10p` > 30,]
pdf('p9vs1.pdf', width = 10,height = 10)
ggplot(cluster91.markers, aes(x = avg_log2FC, y = `-log10p`)) + 
  geom_point() + geom_text_repel(st,mapping = aes(x = avg_log2FC, 
                                                  y = `-log10p`, 
                                                  label = rownames(st),
                                                  color = 'red')) + 
  theme_bw()
dev.off()

# pbmc.markers <- FindAllMarkers(pbmc,
#                                   only.pos = TRUE,
#                                   min.pct = 0.1,
#                                   logfc.threshold = 0.25)

pbmc$seurat_clusters

library(ComplexHeatmap)
pbmc9 <- subset(pbmc, seurat_clusters == 10)
SCT <- pbmc@assays$SCT@data
colnames(SCT) <- NULL

v <- lapply(0:12, function(x){
  d <- data.frame(
    atp1a1 = SCT['ENSMUSG00000033161.10-Atp1a1',],
    stat1 = SCT["ENSMUSG00000026104.14-Stat1",]
  )
  d <- d[rowSums(d) > 1 & pbmc$seurat_clusters == x,]
  cor.test(d$atp1a1,d$stat1,method = 'spearman',exact=FALSE)
})
v

d <- data.frame(
  atp1a1 = SCT['ENSMUSG00000033161.10-Atp1a1',],
  stat1 = SCT["ENSMUSG00000026104.14-Stat1",]
)
d <- d[rowSums(d) > 1 & pbmc$seurat_clusters == 7,]
Heatmap(t(d))
plot(d$atp1a1,d$stat1)
cor.test(d$atp1a1,d$stat1, method = 'spearman',exact=FALSE)

df <- 
  data.frame(lb = pbmc@meta.data$seurat_clusters,
             lb0 = factor(pbmc@meta.data$clusterK12, 
                          levels = unique(pbmc@meta.data$clusterK12)),
             x = pbmc@reductions[["tsne"]]@cell.embeddings[,1],
             y = pbmc@reductions[["tsne"]]@cell.embeddings[,2]
  )
ggplot(df, aes(x = x, y = y, color = lb)) + geom_point()

clkp <- t(table(whole_obj$seurat_clusters,whole_obj$kp))
clkp <- apply(clkp, 2, function(x) x / sum(x))
pdf('cluster_kp.pdf')
Heatmap(clkp,col = circlize::colorRamp2(c(0,1), c("#6495ED","red")))
dev.off()
# ======================== ssgsea ================================
library(GSVA)
library(GSEABase)
library(tidyverse)
library(BiocParallel)
library(stats)

pbmc <- readRDS('HPSC_obj.Rds')
set <- GSEABase::getGmt('D:/GeneSet/mGSKB_Ensembl.gmt')
SCT <- as.matrix(pbmc@assays$SCT@data)
rn <- str_extract(rownames(SCT), pattern = 'ENSMUSG[0-9]+')
SCT <- SCT[!is.na(rn),]
rn <- rn[!is.na(rn)]
rownames(SCT) <- rn

nl <- strsplit(names(set), split = '_')
nl <- unlist(lapply(nl, function(x) x[1]))

table(nl)
BPPARAM = BiocParallel::SnowParam(5)
for(i in c('GO', 'KEGG', 'MSIGDB','REACTOME')){
  KEGG <- set[nl == i]
  KEGG <- filterGeneSets(KEGG, min.sz=5)
  KEGGres <- gsva(SCT, KEGG,method='ssgsea', 
                  parallel.sz = 0, min.sz = 5, 
                  BPPARAM = BPPARAM)
  
  fn <- paste0(i,'.Rds')
  saveRDS(KEGGres,file = fn)
}
# ============================== cellchat1 ============================
library(Seurat)
library(CellChat)
library(patchwork)
library(tidyverse)

pbmc <- readRDS('HPSC_obj.Rds')
data.input = pbmc@assays$SCT@data # normalized data matrix
meta = pbmc@meta.data # a dataframe with rownames containing cell mata data
data.input <- data.input[-grep(rownames(data.input), pattern = 'tdTomato', value = F),]
gn <- gsub(rownames(data.input),pattern = 'ENSMUSG[0-9]+?\\.[0-9]+?-',replacement = '')
data.input <- data.input[!duplicated(gn),]
rownames(data.input) <- gn[!duplicated(gn)]
meta$SCT_snn_res.0.5 <- 1 + as.numeric(as.character(meta$SCT_snn_res.0.5))
cellchat <- createCellChat(object = data.input, 
                           meta = meta, group.by = "SCT_snn_res.0.5")
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "SCT_snn_res.0.5") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
# showDatabaseCategory(CellChatDB)
# dplyr::glimpse(CellChatDB$interaction)
# CellChatDB.use <- subsetDB(CellChatDB, 
#                            search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use = CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)

cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
# saveRDS(cellchat, file = 'cellchat.Rds')
cellchat <- readRDS('cellchat.Rds')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Interaction weights/strength")


mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


pathways.show <- c("CXCL") 
# Hierarchy plot
# Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy")

# Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Chord diagram
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

netAnalysis_contribution(cellchat, signaling = pathways.show)

pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[1,] # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,4) # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)

# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")

# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }


# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), signaling = c("CCL","CXCL"), remove.isolate = FALSE)
# show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
# show all the interactions sending from Inflam.FIB
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(5:11), lab.cex = 0.5,legend.pos.y = 30)
#> Note: The first link end is drawn out of sector 'MIF'.
# show all the interactions received by Inflam.DC
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = 8, legend.pos.x = 15)
# show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_chord_gene(cellchat, sources.use = c(1,2,3,4), targets.use = c(5:11), signaling = c("CCL","CXCL"),legend.pos.x = 8)
#> Note: The second link end is drawn out of sector 'CXCR4 '.
#> Note: The first link end is drawn out of sector 'CXCL12 '.

plotGeneExpression(cellchat, signaling = "CXCL")
plotGeneExpression(cellchat, signaling = "CXCL", enriched.only = FALSE)



# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
#> Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
#> Signaling role analysis on the cell-cell communication network from user's input
gg1 + gg2

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

# =========================== cellchat2 ========================================
library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
# river plot
library(ggalluvial)
netAnalysis_river(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")

nPatterns = 4
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "incoming")
#> Please make sure you have load `library(ggalluvial)` when running this function
# dot plot
netAnalysis_dot(cellchat, pattern = "incoming")


# 相似信号
# reticulate::py_install(packages = 'umap-learn')
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "functional")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space

netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "functional", label.size = 3.5)



cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)
# ============================== monocle ============================
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle)

pbmc <- readRDS('HPSC_obj.Rds')
table(pbmc$seurat_clusters,pbmc$SCT_snn_res.0.5)
##提取表型信息--细胞信息(建议载入细胞的聚类或者细胞类型鉴定信息、实验条件等信息)
expr_matrix <- as(as.matrix(pbmc@assays$SCT@counts), 'sparseMatrix')
##提取表型信息到p_data(phenotype_data)里面 
p_data <- pbmc@meta.data 
p_data$celltype <- pbmc@meta.data$SCT_snn_res.0.5  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
##提取基因信息 如生物类型、gc含量等
f_data <- data.frame(gene_short_name = row.names(pbmc),row.names = row.names(pbmc))
##expr_matrix的行数与f_data的行数相同(gene number), expr_matrix的列数与p_data的行数相同(cell number)

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
#将p_data和f_data从data.frame转换AnnotatedDataFrame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds)))#此时有13714个基因
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10)) #过滤掉在小于10个细胞中表达的基因，还剩11095个基因。


##使用seurat选择的高变基因⚠️
seurat_express_genes <- VariableFeatures(pbmc)

##使用clusters差异表达基因
deg.cluster <- FindAllMarkers(pbmc)
deg.cluster_express_genes <- subset(deg.cluster, p_val_adj < 0.05)$gene
##使用monocle选择的高变基因⚠️
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & 
                       dispersion_empirical >= 1 * dispersion_fit)$gene_id
# saveRDS(cds, file = 'moncole_base.Rds')
# ---------------------------- 所有基因 -----------------------
cds <- readRDS('moncole_base.Rds')
#这一步输入的expressed_genes来自于步骤4。
#⚠️⚠️后续分析使用的是该方法
#也可输入seurat筛选出的高变基因：expressed_genes <- VariableFeatures(pbmc) 
colnames(p_data)
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~ SCT_snn_res.0.5",cores=4) 
#~后面是表示对谁做差异分析的变量，理论上可以为p_data的任意列名
head(diff)

##差异表达基因作为轨迹构建的基因,差异基因的选择标准是qval<0.01,decreasing=F表示按数值增加排序
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

##差异基因的结果文件保存
write.csv(deg,file="expressed_genes_train.monocle.DEG.csv")
## 轨迹构建基因可视化
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
#这一步是很重要的，在我们得到想要的基因列表后，我们需要使用setOrderingFilter将它嵌入cds对象，后续的一系列操作都要依赖于这个list。
#setOrderingFilter之后，这些基因被储存在cds@featureData@data[["use_for_ordering"]]，可以通过table(cds@featureData@data[["use_for_ordering"]])查看
plot_ordering_genes(cds)

#出的图黑色的点表示用来构建轨迹的差异基因，灰色表示背景基因。红色的线是根据第2步计算的基因表达大小和离散度分布的趋势(可以看到，找到的基因属于离散度比较高的基因)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
#⚠️使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
#cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点
head(cds@phenoData@data)
saveRDS(cds, file = 'expressed_genes_monocle.Rds')
# ---------------------------- seurat高变基因基因 -----------------------
cds <- readRDS('moncole_base.Rds')
diff <- differentialGeneTest(cds[seurat_express_genes,],fullModelFormulaStr="~ SCT_snn_res.0.5",cores=4) 
head(diff)
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

write.csv(deg,file="seurat_express_genes_train.monocle.DEG.csv")
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
head(cds@phenoData@data)
saveRDS(cds, file = 'seurat_express_genes_monocle.Rds')
# ---------------------------- moncole高变基因基因 -----------------------
cds <- readRDS('moncole_base.Rds')
diff <- differentialGeneTest(cds[disp.genes,],fullModelFormulaStr="~ SCT_snn_res.0.5",cores=4) 
head(diff)
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

write.csv(deg,file="moncole_genes_train.monocle.DEG.csv")
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
head(cds@phenoData@data)
saveRDS(cds, file = 'moncole_genes_monocle.Rds')

# ---------------------------- cluster差异基因 -----------------------
cds <- readRDS('moncole_base.Rds')
diff <- differentialGeneTest(cds[deg.cluster_express_genes,],fullModelFormulaStr="~ SCT_snn_res.0.5",cores=4) 
head(diff)
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

write.csv(deg,file="cluster_express_genes_train.monocle.DEG.csv")
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
head(cds@phenoData@data)
saveRDS(cds, file = 'cluster_express_genes_monocle.Rds')
# ---------------------------- stemness genes -----------------------
cds <- readRDS('moncole_base.Rds')
library(clusterProfiler)
library(biomaRt)
stemness <- read.gmt('stemness.gmt')

human <- readRDS('humandb.Rds')
mouse <- readRDS('mousedb.Rds')
h2m.driver <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
                     values = unique(stemness$gene), mart = human,
                     attributesL = c("ensembl_gene_id","mgi_symbol"),
                     martL = mouse,uniqueRows = T)

h2m.driver_f <- gsub(rownames(cds),pattern = '\\..*', replacement = '') %in% h2m.driver$Gene.stable.ID 

diff <- differentialGeneTest(cds[h2m.driver_f,],fullModelFormulaStr="~ SCT_snn_res.0.5",cores=4) 
head(diff)
deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

write.csv(deg,file="stemness_genes_train.monocle.DEG.csv")
ordergene <- rownames(deg) 
cds <- setOrderingFilter(cds, ordergene)  
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2, method = 'DDRTree')
cds <- orderCells(cds)
head(cds@phenoData@data)
saveRDS(cds, file = 'stemness_genes_monocle.Rds')
# ======================= monocle analysis ==================================
rm(list = ls())
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle)
library(AUCell)
library(biomaRt)
library(clusterProfiler)

stemness <- read.gmt('stemness.gmt')
# human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl")
# mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl")
# saveRDS(human, file = 'humandb.Rds')
# saveRDS(mouse, file = 'mousedb.Rds')
human <- readRDS('humandb.Rds')
mouse <- readRDS('mousedb.Rds')

stemnessL <- list()
for(i in unique(stemness$term)){
  .driver <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
                    values = stemness[stemness$term == i,2], mart = human,
                    attributesL = c("ensembl_gene_id","mgi_symbol"),
                    martL = mouse,uniqueRows = T)$Gene.stable.ID
  stemnessL[[i]] <- unique(.driver)
}

cds <- readRDS('stemness_genes_monocle.Rds')
deg <- read.csv("stemness_genes_train.monocle.DEG.csv")
ordergene <- deg$X

exp_data <- exprs(cds)
rn <- gsub(rownames(exp_data), pattern = '\\..*', replacement = '')
rownames(exp_data) <- rn
cells_rankings <- AUCell_buildRankings(exp_data)
cells_AUC <- list()
for(i in unique(stemness$term)){
  cells_AUC[[i]] <- getAUC(AUCell_calcAUC(stemnessL[[i]], cells_rankings, 
                                          aucMaxRank=nrow(cells_rankings)*0.1))
}
ns <- names(cells_AUC)
cells_AUC <- purrr::reduce(cells_AUC, rbind)
cells_AUC <- 10^(cells_AUC - 0.1)
rownames(cells_AUC) <- ns
cells_AUC <- t(cells_AUC)
cells_AUC <- as.data.frame(cells_AUC)
pData(cds) <- cbind(pData(cds), cells_AUC)

cds <- orderCells(cds, root_state = 6)
table(cds$State,cds$seurat_clusters)
cds$Pseudotime

plot_cell_trajectory(cds,color_by="Pseudotime ", size=1, show_backbone=TRUE) 
plot_cell_trajectory(cds, color_by = "State",size=1, show_backbone=TRUE)
plot_cell_trajectory(cds, color_by = "seurat_clusters")
plot_cell_trajectory(cds, color_by = "SCT_snn_res.0.5")
plot_cell_trajectory(cds, color_by = "w")
plot_cell_trajectory(cds, color_by = "State") + facet_wrap("~State", nrow = 1)

plot_complex_cell_trajectory(cds[, cds$SCT_snn_res.0.5 == 8],
                             color_by = "seurat_clusters",
                             root_state = 6)




pdf('moncole_tree.pdf')
plot_complex_cell_trajectory(cds, 
                             color_by = "Pseudotime",
                             root_state = 6)
plot_complex_cell_trajectory(cds, 
                             color_by = "w",
                             root_state = 6)
plot_complex_cell_trajectory(cds, 
                             color_by = "KORKOLA_EMBRYONAL_CARCINOMA_UP",
                             root_state = 6)
plot_complex_cell_trajectory(cds, 
                             color_by = "State",
                             root_state = 6)

twdd <- table(cds$State,cds$seurat_clusters)
twdd <- apply(twdd, 1, function(x) x / sum(x))
Heatmap(twdd,col = circlize::colorRamp2(c(0,1), c("#6495ED","red")))

twdd <- t(table(cds$State,cds$seurat_clusters))
twdd <- apply(twdd, 1, function(x) x / sum(x))
Heatmap(twdd,col = circlize::colorRamp2(c(0,1), c("#6495ED","red")))
dev.off()

library(tidyverse)
library(reshape2)
markers <- read.csv('HPSC_marker.csv')
markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> markers0
pbmc <- readRDS('HPSC_obj.Rds')
dim(markers0)
as.data.frame(
  t(as.matrix(pbmc@assays$SCT@scale.data[rownames(pbmc@assays$SCT@scale.data) %in% markers0$X,]))
) -> datagene
dim(datagene)
datagene$cluster <- pbmc$seurat_clusters 
datageneL <- melt(datagene)
datageneL %>% group_by(cluster,variable) %>% 
  summarise(mean_expression = mean(value)) -> 
  datageneL0
datageneL0 %>% group_by(variable) %>% 
  mutate(mean_expression0 = as.vector(scale(mean_expression))) -> 
  datageneL00

library(viridis)
datageneL00$variable <- gsub(datageneL00$variable,
                             pattern = 'ENSMUSG[.0-9]+-',
                             replacement = '')


markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> markers2

library(Seurat)
# pdf('feature.pdf', width = 10, height = 10)
# dt <- DoHeatmap(pbmc, features = markers0$gene, size = 5, angle = -50, hjust=0.8) + NoLegend()
# dev.off()

# pdf('featurexxxx.pdf', width = 10, height = 10)
# dt <- DoHeatmap(pbmc, features = markers0$gene, size = 5, angle = -50, hjust=0.8) + NoLegend()
# dev.off()

markers0 <- markers0[markers0$p_val_adj < 0.9,]
pdf('fearture.boubpl0.pdf',width   = 15)
DotPlot(pbmc, features = markers0$gene)+  theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  scale_x_discrete("")+scale_y_discrete("") + scale_color_viridis(
    option="C",end = 0.8)
dev.off()

dt <- DoHeatmap(pbmc, features = markers0$gene, size = 5, angle = -50, hjust=0.8) + NoLegend()
dt0 <- dt$data
dt0 <- na.omit(dt0)
dt0 %>% group_by(Identity,Feature) %>% 
  summarise(mean_expression = mean(Expression)) -> 
  datageneLdt
View(datageneLdt)
datageneLdt$Feature <- gsub(datageneLdt$Feature,
                            pattern = 'ENSMUSG[.0-9]+-',
                            replacement = '')
datageneLdt %>% group_by(Feature) %>% mutate(
  scalee = as.vector(scale(mean_expression))) -> datageneLdt
colnames(dt$data)
pdf('fearture.boubpl.pdf',height  = 12)
ggplot(datageneLdt,aes(x = Identity,
                       y = Feature,
                       color = scalee,
                       size = scalee)) +
  geom_point() + 
  scale_color_viridis(option="C",begin = 0.2) + theme_bw()
dev.off()




head(diff)
table(cds$State,cds$seurat_clusters)

pdata <- Biobase::pData(cds)
s.cells <- subset(pdata, State=="10") %>% rownames()

##选择前4个top基因并将其对象取出
keygenes <- head(ordergene,10)
cds_subset <- cds[keygenes,]
##可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset, color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset, color_by = "SCT_snn_res.0.5")
p3 <- plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")
plotc <- p1|p2|p3
plotc

#指定基因
s.genes <- c("ENSMUSG00000033161.10-Atp1a1")
p1 <- plot_genes_jitter(cds[s.genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(cds[s.genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(cds[s.genes,], color_by = "State")
p4 <- plot_complex_cell_trajectory(cds[s.genes,], x = 1, y = 2,
                                   color_by = "celltype")
plotc <- p1|p2|p3|p4
plotc

#这里是把排序基因（ordergene）提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds[ordergene,], cores = 4, 
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_genes <- top_n(Time_diff, n = 10, desc(qval)) %>% 
  pull(gene_short_name) %>% as.character()
p = plot_pseudotime_heatmap(cds[c(Time_genes,'ENSMUSG00000033161.10-Atp1a1'),], 
                            num_clusters=4, show_rownames=T, 
                            return_heatmap=T)
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes, c("gene_short_name", "pval", "qval")]
marker_genes <- row.names(subset(fData(cds),
                                 gene_short_name %in% c(
                                   "ENSMUSG00000000402.2-Egfl6",
                                   'ENSMUSG00000033161.10-Atp1a1')))

diff_test_res <- differentialGeneTest(cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)

plot_cell_trajectory(cds, color_by = "State")

# BEAM_res1 <- BEAM(cds[ordergene,], branch_point = 1, cores = 4)
# BEAM_res2 <- BEAM(cds[ordergene,], branch_point = 2, cores = 4)
# BEAM_res3 <- BEAM(cds[ordergene,], branch_point = 3, cores = 4)
# BEAM_res4 <- BEAM(cds[ordergene,], branch_point = 4, cores = 4)
# BEAM_res5 <- BEAM(cds[ordergene,], branch_point = 5, cores = 4)
# BEAM_res6 <- BEAM(cds[ordergene,], branch_point = 6, cores = 4)
# BEAM_res7 <- BEAM(cds[ordergene,], branch_point = 7, cores = 4)
# BEAM_res <- list(BEAM_res1 = BEAM_res1,
#                  BEAM_res2 = BEAM_res2,
#                  BEAM_res3 = BEAM_res3,
#                  BEAM_res4 = BEAM_res4,
#                  BEAM_res5 = BEAM_res5,
#                  BEAM_res6 = BEAM_res6,
#                  BEAM_res7 = BEAM_res7)
# saveRDS(BEAM_res, file = 'BEAM_res.Rds')

BEAM_res <- readRDS('BEAM_res.Rds')
BEAM_res6 <- BEAM_res[[4]]
BEAM_res6 <- BEAM_res6[order(BEAM_res6$qval),]
BEAM_res6 <- BEAM_res6[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res6,
                                                 qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 2,
                            use_gene_short_name = T,
                            show_rownames = T)#有632个gene，太多了

#选前100个基因可视化
grep(BEAM_res6$gene_short_name, pattern = 'atp1a1',ignore.case = T)
BEAM_res6[grep(BEAM_res6$gene_short_name, pattern = 'atp1a1',ignore.case = T),]
BEAM_genes6t <- top_n(BEAM_res6, n = 100, desc(qval)) %>% 
  pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes6t,],  
                                 branch_point = 2, 
                                 num_clusters = 3, 
                                 show_rownames = T, 
                                 return_heatmap = T)
#显著差异基因按热图结果排序并保存
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- BEAM_res6[hp.genes, c("gene_short_name", "pval", "qval")]
genes <- row.names(subset(fData(cds),
                          gene_short_name %in% c( "ENSMUSG00000000402.2-Egfl6",
                                                  'ENSMUSG00000033161.10-Atp1a1')))

plot_genes_branched_pseudotime(cds[genes,],
                               branch_point = 1,
                               color_by = "State",
                               ncol = 1)

# ============================== infcnv ============================
library(infercnv)
library(dplyr)
library(Seurat)
library(patchwork)
library(tidyverse)

dir.create('infcnv')
pbmc <- readRDS('HPSC_obj.Rds')
data0 <- data.frame(names(pbmc@active.ident), pbmc@active.ident, 
                    Tlv = str_extract(names(pbmc@active.ident),
                                      pattern = 'T[0-9]'))
cm <- as.matrix(table(data0$Tlv,data0$pbmc.active.ident))
class(cm)
t(cm) / colSums(cm)

data <- data.frame(names(pbmc@active.ident), pbmc@active.ident)
write.table(data,file = 'infcnv/sample.txt', row.names = F, col.names = F,quote = F,sep='\t')


raw_counts_matrixcounts <- as.matrix(pbmc@assays$SCT@counts)
write.csv(raw_counts_matrixcounts, file = 'infcnv/SCTcount.csv')

setwd('infcnv')
grep(rownames(raw_counts_matrixcounts),pattern = 'td',value = T)
raw_counts_matrixcounts <- raw_counts_matrixcounts[rownames(raw_counts_matrixcounts) != "tdTomato-tdTomato",]
rownames(raw_counts_matrixcounts) <- str_extract(rownames(raw_counts_matrixcounts), pattern = 'ENSMUSG[0-9]+')
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = raw_counts_matrixcounts, # 可以直接提供矩阵对象
                                    annotations_file="sample.txt",
                                    delim="\t",
                                    gene_order_file="ENSG_position.txt",
                                    ref_group_names = c('11'))

# perform infercnv operations to reveal cnv signal
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=1,  # use 1 for smart-seq, 0.1 for 10x-genomics
                             out_dir="output_dir_new",  # 输出文件夹
                             denoise=T, #去噪
                             HMM=T,
                             output_format = "pdf") # 是否基于HMM预测CNV
saveRDS(infercnv_obj, file = '../infercnv_new_obj.Rds')
# ======================================  AUCell  ===============================================
library(AUCell)
library(clusterProfiler)
library(ggplot2)
library(Seurat)
library(SeuratData)
library(GSEABase)
library(tidyverse)
library(GSVA)

##seurat RDS object
whole_obj <- readRDS('HPSC_obj.Rds')
SCT <- as.matrix(whole_obj@assays$SCT@data)
rn <- str_extract(rownames(SCT), pattern = 'ENSMUSG[0-9]+')
SCT <- SCT[!is.na(rn),]
rn <- rn[!is.na(rn)]
rownames(SCT) <- rn
cells_rankings <- AUCell_buildRankings(SCT)  # 关键一步
set <- GSEABase::getGmt('D:/GeneSet/mGSKB_Ensembl.gmt')


nl <- strsplit(names(set), split = '_')
nl <- unlist(lapply(nl, function(x) x[1]))

for(i in c('GO', 'KEGG', 'MSIGDB','REACTOME')){
  geneSets <- set[nl == i]
  geneSets <- filterGeneSets(geneSets, min.sz=5)
  
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
  fn <- paste0(i,'_AUCell.Rds')
  saveRDS(cells_AUC,file = fn)
}
# =============================== aucell cor ========================================
library(ggrepel)
REACTOME_AUCell <- readRDS('REACTOME_AUCell.Rds')
GO_AUCell <- readRDS('GO_AUCell.Rds')
KEGG_AUC <- readRDS('KEGG_AUCell.Rds')
MSIGDB_AUCell <- readRDS('MSIGDB_AUCell.Rds')

aucs <- as.matrix(getAUC(REACTOME_AUCell)[rownames(REACTOME_AUCell@assays@data$AUC),])
AUC <- REACTOME_AUCell@assays@data$AUC
AUCcp <- apply(AUC, 1, function(x) {
  tmpd <- data.frame(x = x / max(x),
                     Atp1a1 = SCT["ENSMUSG00000033161.10-Atp1a1",] / max(SCT["ENSMUSG00000033161.10-Atp1a1",]))
  tmpd <- tmpd[rowSums(tmpd) > 0.2,]
  tmpc <- cor.test(tmpd[,1],
                   tmpd[,2],
                   method = 'spearman',exact=FALSE)
  c(tmpc$p.value,tmpc$estimate)
  
})
AUCcp <- t(AUCcp)

SCTt <- SCT
SCTt <- SCTt[rowSums(SCTt) > 100,]
SCTtcp <- apply(SCTt, 1, function(x) {
  tmpd <- data.frame(x = x / max(x),
                     Atp1a1 = SCT["ENSMUSG00000033161.10-Atp1a1",] / max(SCT["ENSMUSG00000033161.10-Atp1a1",]))
  tmpd <- tmpd[rowSums(tmpd) > 0,]
  tmpc <- cor.test(tmpd[,1],
                   tmpd[,2],
                   method = 'spearman',exact=FALSE)
  c(tmpc$p.value,tmpc$estimate)
  
})
SCTtcp <- t(SCTtcp)

AUC_GO <- GO_AUCell@assays@data$AUC
AUC_GOcp <- apply(AUC_GO, 1, function(x) {
  tmpd <- data.frame(x = x / max(x),
                     Atp1a1 = SCT["ENSMUSG00000033161.10-Atp1a1",] / max(SCT["ENSMUSG00000033161.10-Atp1a1",]))
  tmpd <- tmpd[rowSums(tmpd) > 0.2,]
  tmpc <- cor.test(tmpd[,1],
                   tmpd[,2],
                   method = 'spearman',exact=FALSE)
  c(tmpc$p.value,tmpc$estimate)
  
})
AUC_GOcp <- t(AUC_GOcp)

AUC_KEGG <- KEGG_AUC@assays@data$AUC
AUC_KEGGcp <- apply(AUC_KEGG, 1, function(x) {
  tmpd <- data.frame(x = x / max(x),
                     Atp1a1 = SCT["ENSMUSG00000033161.10-Atp1a1",] / max(SCT["ENSMUSG00000033161.10-Atp1a1",]))
  tmpd <- tmpd[rowSums(tmpd) > 0.2,]
  tmpc <- cor.test(tmpd[,1],
                   tmpd[,2],
                   method = 'spearman',exact=FALSE)
  c(tmpc$p.value,tmpc$estimate)
  
})
AUC_KEGGcp <- t(AUC_KEGGcp)

AUC_MSIGDB <- MSIGDB_AUCell@assays@data$AUC
AUC_MSIGDBcp <- apply(AUC_MSIGDB, 1, function(x) {
  tmpd <- data.frame(x = x / max(x),
                     Atp1a1 = SCT["ENSMUSG00000033161.10-Atp1a1",] / max(SCT["ENSMUSG00000033161.10-Atp1a1",]))
  tmpd <- tmpd[rowSums(tmpd) > 0.2,]
  tmpc <- cor.test(tmpd[,1],
                   tmpd[,2],
                   method = 'spearman',exact=FALSE)
  c(tmpc$p.value,tmpc$estimate)
  
})
AUC_MSIGDBcp <- t(AUC_MSIGDBcp)

ttmp <- list(
  geneSet_ATP1A1 = AUCcp,
  gene_ATP1A1 = SCTtcp,
  AUC_GOcp = AUC_GOcp,
  AUC_MSIGDBcp = AUC_MSIGDBcp,
  AUC_KEGGcp = AUC_KEGGcp
)
# saveRDS(ttmp, file = 'Cor_set_gene.Rds')
# tmp <- readRDS('Cor_set_gene.Rds')
# AUCcp <- tmp$geneSet_ATP1A1
# SCTtcp <- tmp$gene_ATP1A1

grep(rownames(whole_obj@assays$SCT@data), pattern = 'runx', value = T,ignore.case = T)
AUC_GOcp[grep(rownames(AUC_GOcp),pattern = 'hypoxia',ignore.case = T),]
AUCcp[grep(rownames(AUCcp),pattern = 'stat',ignore.case = T),]
AUCcp[grep(rownames(AUCcp),pattern = 'nkx',ignore.case = T),]

whole_obj <- readRDS('HPSC_obj.Rds')

d <- 
  data.frame(
    tumor = whole_obj$tumor,
    w = whole_obj$w,
    seurat_clusters = whole_obj$seurat_clusters,
    atp1a1 = whole_obj@assays$SCT@data['ENSMUSG00000033161.10-Atp1a1',],
    stat1 = whole_obj@assays$SCT@data['ENSMUSG00000026104.14-Stat1',],
    Nkx2 = whole_obj@assays$SCT@data['ENSMUSG00000001496.15-Nkx2-1',],
    Runx2 = whole_obj@assays$SCT@data['ENSMUSG00000039153.16-Runx2',],
    set0 = AUC['REACTOME_MM_STAT1_PATHWAY_Ensembl',],
    set = AUC_GO['GO_BP_MM_CELLULAR_RESPONSE_TO_HYPOXIA_Ensembl',]
  )

d %>% group_by(seurat_clusters) %>% 
  summarise(mean(atp1a1),mean(Nkx2),mean(Runx2),mean(stat1),mean(set0)) -> tmp

ggplot(tmp, aes(x = `mean(atp1a1)`, y = `mean(set0)`)) + 
  geom_point() + geom_text_repel(mapping = aes(x = `mean(atp1a1)`, 
                                               y = `mean(set0)`, 
                                               label = seurat_clusters,
                                               color = 'red')) + 
  theme_bw()
tmp1 <- tmp[as.numeric(as.character(tmp$seurat_clusters)) %in% c(0,1,4,9),]
ggplot(tmp1, aes(x = `mean(atp1a1)`, y = `mean(set0)`)) + 
  geom_point() + geom_text_repel(mapping = aes(x = `mean(atp1a1)`, 
                                               y = `mean(set0)`, 
                                               label = seurat_clusters,
                                               color = 'red')) + 
  theme_bw()
cor.test(tmp2$`mean(atp1a1)`, tmp2$`mean(set0)`, method = 'spearman')
# =============================== aucell thread =======================================
library(AUCell)
library(tidyverse)

whole_obj <- readRDS('HPSC_obj.Rds')
REACTOME_AUCell <- readRDS('REACTOME_AUCell.Rds')
GO_AUCell <- readRDS('GO_AUCell.Rds')
KEGG_AUC <- readRDS('KEGG_AUCell.Rds')
MSIGDB_AUCell <- readRDS('MSIGDB_AUCell.Rds')


# cells_assignment <- AUCell_exploreThresholds(REACTOME_AUCell, 
#                                              plotHist=F, assign=TRUE) 
# # saveRDS(cells_assignment, file = 'REACTOME_AUCell_thre.Rds')
# GO_AUCell_assignment <- AUCell_exploreThresholds(GO_AUCell, 
#                                              plotHist=F, assign=TRUE) 
# saveRDS(GO_AUCell_assignment, file = 'GO_AUCell_thre.Rds')
# KEGG_AUC_assignment <- AUCell_exploreThresholds(KEGG_AUC, 
#                                                  plotHist=F, assign=TRUE) 
# saveRDS(KEGG_AUC_assignment, file = 'KEGG_AUC_thre.Rds')
# MSIGDB_AUCell_assignment <- AUCell_exploreThresholds(MSIGDB_AUCell, 
#                                                  plotHist=F, assign=TRUE) 
# saveRDS(MSIGDB_AUCell_assignment, file = 'MSIGDB_AUCell_thre.Rds')

# cells_assignmentv <- lapply(cells_assignment, function(x) length(x$assignment))
# fb100 <- unlist(cells_assignmentv) > 100
# grep(names(cells_assignment[fb100]),pattern = 'stat1',value = T,ignore.case = T)
# aucdata <- as.data.frame(t(getAUC(REACTOME_AUCell)))
# gathlist <- list()
# for(i in 1:2565){
#   tdt <- data.frame(geneset = aucdata[,i], group = whole_obj$seurat_clusters)
#   colnames(tdt) <- c('set', 'seurat_clusters')
#   gathlist[[colnames(aucdata)[i]]] <- aov(set ~ seurat_clusters, tdt)
# }
# pgathlist <- lapply(gathlist, function(x) summary(x)[[1]][["Pr(>F)"]][1])
# sum(pgathlist < 0.0001)
aucdata <- as.data.frame(t(getAUC(REACTOME_AUCell)))
aucdata$cluster <- whole_obj$seurat_clusters
ggplot(aucdata, aes(x = seurat_clusters, 
                    y = REACTOME_MM_STAT1_PATHWAY_Ensembl,
                    fill = seurat_clusters)) + 
  geom_violin() + theme_bw()
# ===================================================
library(reshape2)
cnv_table <- read.table("./infcnv/output_dir/infercnv.observations.txt", 
                        header=T)
# Score cells based on their CNV scores 
# Replicate the table 
cnv_score_table <- as.matrix(cnv_table)
cnv_score_mat <- as.matrix(cnv_table)
# Scoring
# CNV的5分类系统
cnv_score_table[cnv_score_mat >= 0 & cnv_score_mat < 0.3] <- "A" #complete loss. 2pts
cnv_score_table[cnv_score_mat >= 0.3 & cnv_score_mat < 0.7] <- "B" #loss of one copy. 1pts
cnv_score_table[cnv_score_mat >= 0.7 & cnv_score_mat < 1.3] <- "C" #Neutral. 0pts
cnv_score_table[cnv_score_mat >= 1.3 & cnv_score_mat <= 1.5] <- "D" #addition of one copy. 1pts
cnv_score_table[cnv_score_mat > 1.5 & cnv_score_mat <= 2] <- "E" #addition of two copies. 2pts
cnv_score_table[cnv_score_mat > 2] <- "F" #addition of more than two copies. 2pts

# cnv_score_table[cnv_score_mat >= 0 & cnv_score_mat < 2.5] <- "A" 
# cnv_score_table[cnv_score_mat >= 2.5 & cnv_score_mat < 3.5] <- "B" 
# cnv_score_table[cnv_score_mat >= 3.5 & cnv_score_mat < 4.5] <- "C" 
# cnv_score_table[cnv_score_mat >= 4.5 & cnv_score_mat < 5.5] <- "D" 
# Check
table(cnv_score_table[,2])
# Replace with score 
cnv_score_table_pts <- cnv_table

#  把5分类调整为 3分类系统
cnv_score_table_pts[cnv_score_table == "A"] <- 2
cnv_score_table_pts[cnv_score_table == "B"] <- 1
cnv_score_table_pts[cnv_score_table == "C"] <- 0
cnv_score_table_pts[cnv_score_table == "D"] <- 1
cnv_score_table_pts[cnv_score_table == "E"] <- 2
cnv_score_table_pts[cnv_score_table == "F"] <- 2

# cnv_score_table_pts[cnv_score_table == "A"] <- 1
# cnv_score_table_pts[cnv_score_table == "B"] <- 0
# cnv_score_table_pts[cnv_score_table == "C"] <- 1
# cnv_score_table_pts[cnv_score_table == "D"] <- 2
# Scores are stored in “cnv_score_table_pts”. Use colSums to add up scores for each cell and store as vector 
cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
colnames(cell_scores_CNV) <- "cnv_score"
head(cell_scores_CNV)
write.csv(x = cell_scores_CNV, file = "cnv_scores.csv")

cell_scores_CNV$rn <- rownames(cell_scores_CNV)
f <- rownames(whole_obj@meta.data) %in% rownames(cell_scores_CNV)
dcnv <- 
  data.frame(
    cluster = whole_obj@meta.data$seurat_clusters,
    genetype = whole_obj@meta.data$kp,
    time = whole_obj@meta.data$w,
    tumor = whole_obj@meta.data$tumor)
rownames(dcnv) <- rownames(whole_obj@meta.data)
dcnv <- dcnv[f,]
dcnv$rn <- rownames(dcnv)
dcnv %>% left_join(cell_scores_CNV) -> dcnv

L <- melt(cnv_score_mat)
pdf('cnvsupp.pdf')
ggplot(L, aes(x = value,fill = 'red')) + geom_density() + theme_bw()
ggplot(dcnv, aes(x = cluster, y = cnv_score, fill = cluster)) + 
  geom_boxplot() + theme_bw()
ggplot(dcnv, aes(x = genetype, y = cnv_score, fill = genetype)) + 
  geom_boxplot() + theme_bw()
ggplot(dcnv, aes(x = time, y = cnv_score, fill = time)) + 
  geom_boxplot() + theme_bw()
ggplot(dcnv, aes(x = tumor, y = cnv_score, fill = tumor)) + 
  geom_boxplot() + theme_bw()

ggplot(dcnv, aes(x = cnv_score, fill = cluster)) + 
  geom_density(alpha = 0.3) + theme_bw()
ggplot(dcnv, aes(x = cnv_score, fill = genetype)) + 
  geom_density(alpha = 0.3) + theme_bw()
ggplot(dcnv, aes(x = cnv_score, fill = time)) + 
  geom_density(alpha = 0.3) + theme_bw()
ggplot(dcnv, aes(x = cnv_score, fill = tumor)) + 
  geom_density(alpha = 0.3) + theme_bw()
dev.off()

dcnv %>% group_by(cluster) %>% summarise(cnv_score = mean(cnv_score)) -> cnv_score1
data.frame(Atp1a1 = SCT['ENSMUSG00000033161.10-Atp1a1',],
           group = pbmc$seurat_clusters) %>% group_by(group) %>% summarise(mean(Atp1a1)) -> atp1a12
atp1a12 <- atp1a12[c(1:11,13),]
cbind(cnv_score = cnv_score1,
      mean_atp1a1 = atp1a12) -> cnvatp1a1
pdf('atp1a1_cnv.pdf')
ggplot(cnvatp1a1,aes(x = cnv_score.cnv_score, y = `mean_atp1a1.mean(Atp1a1)`)) +
  geom_point() + geom_text_repel(mapping = aes(x = cnv_score.cnv_score, 
                                               y = `mean_atp1a1.mean(Atp1a1)`, 
                                               label = cnv_score.cluster,
                                               color = 'red')) + 
  theme_bw() + xlab('CNV_Score') + ylab('mean_Atp1a1')
dev.off()

# ===================================================
options(stringsAsFactors = F)
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Seurat)
library(infercnv)
library(miscTools)
whole_obj <- readRDS('HPSC_obj.Rds')
#  Import inferCNV dendrogram
infercnv.dend <- read.dendrogram(file = "infcnv/output_dir_new/infercnv.observations_dendrogram.txt")
# Cut tree 
infercnv.labels <- cutree(infercnv.dend, k = 9, order_clusters_as_data = FALSE)

nsd <- 
  data.frame(
    lable = infercnv.labels,
    rn = names(infercnv.labels))
dcnv <- 
  data.frame(
    rn = rownames(whole_obj@meta.data),
    cluster = whole_obj@meta.data$seurat_clusters,
    genetype = whole_obj@meta.data$kp,
    time = whole_obj@meta.data$w,
    tumor = whole_obj@meta.data$tumor)
nsd %>% left_join(dcnv) -> nsd
nsd$lable <- factor(nsd$lable)
table(nsd$lable, nsd$cluster)
table(nsd$lable, nsd$time)
table(nsd$lable, nsd$tumor)
table(nsd$lable, nsd$genetype)

pdf('cnv_cluster.pdf')
ggplot(nsd, aes(x = lable, fill = cluster)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = lable, fill = time)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = lable, fill = tumor)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = lable, fill = genetype)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()

ggplot(nsd, aes(x = genetype, fill = lable)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = tumor, fill = lable)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = cluster, fill = lable)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
dev.off()

# Color labels
the_bars <- as.data.frame(tableau_color_pal("Tableau 20")(20)[infercnv.labels])
colnames(the_bars) <- "inferCNV_tree"
the_bars$inferCNV_tree <- as.character(the_bars$inferCNV_tree)

pdf('cnv_tree.pdf',height = 10)
infercnv.dend %>% set("labels",rep("", 
                                   nobs(infercnv.dend)))  %>% 
  plot(main="inferCNV dendrogram") %>%
  colored_bars(colors = as.data.frame(the_bars), 
               dend = infercnv.dend, 
               sort_by_labels_order = F, 
               add = T, y_scale=100 , y_shift = 0)
dev.off()

cd <- data.frame(cnv = the_bars$inferCNV_tree,
                 lables = infercnv.labels)
cd$lables <- factor(cd$lables, levels = names(table(cd$lables)))
v <- names(table(cd$cnv))
names(v) <- apply(table(cd$cnv, cd$lables), 1, function(x) which(x != 0))
pdf('color_bar.pdf')
ggplot(cd, aes(x = lables, fill = lables)) + 
  geom_bar() + 
  scale_fill_manual(values = v) +
  theme_bw()
dev.off()


# ================= 分辨率补图 ===========================
# library(clustree)
# sce <- FindClusters(
#   object = pbmc,
#   resolution = c(seq(.3,1.6,.2))
# )
# clustree(sce@meta.data, prefix = "SCT_snn_res.")
# ==========================================================
library(Seurat)
library(ggplot2)
library(AUCell)
library(tidyverse)
library(ggrepel)

whole_obj <- readRDS('HPSC_obj.Rds')
REACTOME_AUCell <- readRDS('REACTOME_AUCell.Rds')
GO_AUCell <- readRDS('GO_AUCell.Rds')
KEGG_AUC <- readRDS('KEGG_AUCell.Rds')
MSIGDB_AUCell <- readRDS('MSIGDB_AUCell.Rds')

grep(rownames(whole_obj@assays$SCT@data),pattern = 'atp1a1',
     ignore.case = T, value = T)
grep(rownames(whole_obj@assays$SCT@data),pattern = 'stat',
     ignore.case = T, value = T)
grep(rownames(getAUC(MSIGDB_AUCell)),pattern = 'stat1',
     ignore.case = T, value = T)

pdf('MSIGDB_MM_BAKER_HEMATOPOESIS_STAT1_TARGETS_Ensembl__KEGG_MM_JAK-STAT_SIGNALING_PATHWAY_Ensembl.pdf')
##
(grep(rownames(getAUC(MSIGDB_AUCell)),pattern = 'MSIGDB_MM_BAKER_HEMATOPOESIS_STAT1_TARGETS_Ensembl',
      ignore.case = T, value = F) -> setn)
##
grep(rownames(whole_obj@assays$SCT@data),pattern = 'atp1a1',
     ignore.case = T, value = F) -> genen

d <- data.frame(ATP1A1 = whole_obj@assays$SCT@data[genen,],
                stat = getAUC(MSIGDB_AUCell)[setn,],
                cluster = whole_obj$seurat_clusters)

d0 <- d[d$cluster != 11,]
d %>% group_by(cluster) %>% summarise(mean_atp1a1 = mean(ATP1A1), 
                                      mean_stat = mean(stat)) -> tmp
ggplot(tmp,aes(x = mean_stat, y = mean_atp1a1)) + geom_point() + 
  geom_text_repel(mapping = aes(x = mean_stat, 
                                y = mean_atp1a1, 
                                label = cluster,
                                color = 'red')) + theme_bw()

ggplot(d,aes(x = cluster, y = stat, fill = cluster)) + geom_violin() + theme_bw()

(grep(rownames(getAUC(KEGG_AUC)),pattern = 'KEGG_MM_JAK-STAT_SIGNALING_PATHWAY_Ensembl',
      ignore.case = T, value = F) -> setn)
d <- data.frame(ATP1A1 = whole_obj@assays$SCT@data[genen,],
                stat = getAUC(MSIGDB_AUCell)[setn,],
                cluster = whole_obj$seurat_clusters)
d %>% group_by(cluster) %>% summarise(mean_atp1a1 = mean(ATP1A1), 
                                      mean_stat = mean(stat)) -> tmp
ggplot(tmp,aes(x = mean_stat, y = mean_atp1a1)) + geom_point() + 
  geom_text_repel(mapping = aes(x = mean_stat, 
                                y = mean_atp1a1, 
                                label = cluster,
                                color = 'red')) + theme_bw()

ggplot(d,aes(x = cluster, y = stat, fill = cluster)) + geom_violin() + theme_bw()
dev.off()

# ---------------第一列大于0，第二列大于0-----------
# d <- d[d[,1] > 0 & d[,2] > 0,]

# ---------------第一列大于0，第二列大于0-----------
# d <- d[rowSums(d) > 0 ,]

# library(RColorBrewer)
# colormap<-rev(brewer.pal(10,'Spectral'))
# pdf('helx.pdf')
# ggplot(d,aes(x = atp1a1, y = gene)) + 
#   theme_bw() + ylab('Stat1') + xlab('Atp1a1') + 
#   geom_hex(bins=30,na.rm=F)+ 
#   scale_fill_gradientn(colours=colormap,trans='log10')
# dev.off()

# pdf("NKA-PGC1.PDF")
# cor.test(
#   d[,"ATP1A1"],
#   d[,"Ppargc1a"])
# dev.off()
# sum(whole_obj@assays$SCT@data["ENSMUSG00000033161.10-Atp1a1",] == 0 &
#   whole_obj@assays$SCT@data["ENSMUSG00000029167.13-Ppargc1a",] == 0)
# 
# 
# pdf('violin.Ppargc1a.pdf',height = 10)
# VlnPlot(whole_obj, 
#         features = c(
#           "ENSMUSG00000033161.10-Atp1a1",
#           "ENSMUSG00000029167.13-Ppargc1a",
#           'ENSMUSG00000039153.16-Runx2',
#           'ENSMUSG00000001496.15-Nkx2-1',
#           'ENSMUSG00000071552.4-Tigit'),
#         ncol = 1)
# dev.off()












library(Seurat)
library("survival")
library(clusterProfiler)
library(biomaRt)
library(tidyverse)

cox_model <- readRDS('res.cox_model.Rds')
HPSC <- readRDS('../../肺癌数据/HPSC/scRNA/HPSC_obj.Rds')
View(HPSC)
# data.frame(HPSC$seurat_clusters,
#            names(HPSC$seurat_clusters)) -> colord
# write.csv(colord,'colord.csv')
sctdata <- HPSC@assays$SCT@data
rownames(sctdata)
data.frame(UMAP1 = HPSC@reductions$umap@cell.embeddings[,1],
           UMAP2 = HPSC@reductions$umap@cell.embeddings[,2],
           week = HPSC$w,
           cluster = HPSC$seurat_clusters) -> umapp

umapp$week <- factor(umapp$week,levels =  c('0w','2w','4w','12w','18w','20w','30w'))
ggplot(umapp,aes(x = UMAP1, y = UMAP2, color = week)) + 
  theme_bw() + geom_point() -> umapweek
# ggsave(umapweek,filename = 'umap_week.pdf')
DimPlot(HPSC,reduction = "umap",group.by = 'w')

tsmep <- data.frame(UMAP1 = HPSC@reductions$tsne@cell.embeddings[,1],
                    UMAP2 = HPSC@reductions$tsne@cell.embeddings[,2],
                    week = HPSC$w,
                    cluster = HPSC$seurat_clusters)

pdf('new_tsne_umap.pdf',width = 5,height = 5)
ggplot(tsmep,aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  theme_bw()+ 
  stat_density2d(contour_var = "ndensity",size = 0.01) + 
  geom_point(size = 0.5) + ylim(-45,45) + xlim(-55,55)

ggplot(umapp,aes(x = UMAP1, y = UMAP2, color = cluster)) + 
  theme_bw()+ 
  stat_density2d(contour_var = "ndensity",size = 0.01) + 
  geom_point(size = 0.5)  + ylim(-8,10) + xlim(-6,10)
dev.off()

mouse.gene <- gsub(rownames(sctdata), pattern = '.+?-',replacement = '')
# human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl")
# mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl")
# m2h.g <- getLDS(attributes = c("mgi_symbol"),filters = "mgi_symbol",
#                 values = mouse.gene, mart = mouse,
#                 attributesL = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position"),
#                 martL = human,uniqueRows = T)
# saveRDS(m2h.g, file = 'm2h.Rds')

m2h.g <- readRDS('m2h.Rds')
sctdata <- HPSC@assays$SCT@data
kmg <- m2h.g$HGNC.symbol[m2h.g$HGNC.symbol %in% names(cox_model$coefficients)]
names(kmg) <- m2h.g$MGI.symbol[m2h.g$HGNC.symbol %in% names(cox_model$coefficients)]
fg <- mouse.gene %in% names(kmg)
sctdata <- sctdata[fg,]
rownames(sctdata) <- kmg[gsub(rownames(sctdata), pattern = '.+?-',replacement = '')]
sctdata <- as.data.frame(t(as.matrix(sctdata)))
names(cox_model$coefficients)[!names(cox_model$coefficients) %in% colnames(sctdata)]
gn <- colnames(sctdata)
apply(sctdata, 2, sum)

sctdata <- apply(sctdata, 2, scale)
sctdata <- as.data.frame(sctdata)
apply(sctdata, 2, mean)
cox_model$coefficients
sctdata$OGFRP1 <- 0
sctdata$CHRNA5 <- 0
sctdata$FAM207A <- 0
sctdata$LSP1P4 <- 0
sctdata$C20orf197 <- 0
rownames(sctdata) <- colnames(HPSC@assays$SCT@data)
predictans <- predict(cox_model,sctdata)
# saveRDS(predictans, file = 'predictans_ans.Rds')
predictans <- readRDS('predictans_ans.Rds')
HPSC <- AddMetaData(HPSC,predictans,'cox_predictans')
# FeaturePlot(HPSC, features = 'cox_predictans')
VlnPlot(object = HPSC, 
        features = 'cox_predictans') + ylim(-3,3)
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 0]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 1]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 2]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 3]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 4]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 5]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 6]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 7]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 8]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 9]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 10]) -> tmp
tmp$p.value
wilcox.test(HPSC$cox_predictans[HPSC$seurat_clusters == 11],
            HPSC$cox_predictans[HPSC$seurat_clusters == 12]) -> tmp
tmp$p.value

# week genetype c
data.frame(week = HPSC$w,
           genetype = HPSC$kp,
           tumor = HPSC$tumor,
           cluster = HPSC$seurat_clusters,
           cox_predictans = HPSC$cox_predictans) -> dotheatmap
dotheatmap$week <- factor(dotheatmap$week, levels = c('0w','2w','4w','12w','18w','20w','30w'))

# pdf('mean.predict.pdf')
# dotheatmap %>% group_by(cluster,week) %>% summarise(mean = mean(cox_predictans),
#                                                     siz = length(cox_predictans)) -> 
#   tmp1
# ggplot(tmp1, aes(x = cluster, y =week, size = siz, color = mean)) + geom_point() +
#   scale_color_gradient2(low = "#6495ED",
#                         mid = 'grey',
#                         high = 'red') + theme_bw()
# 
# dotheatmap %>% group_by(cluster,genetype) %>% summarise(mean = mean(cox_predictans),
#                                                         siz = length(cox_predictans)) -> 
#   tmp2
# ggplot(tmp2, aes(x = cluster, y =genetype, size = siz, color = mean)) + geom_point() +
#   scale_color_gradient2(low = "#6495ED",
#                         mid = 'grey',
#                         high = 'red') + theme_bw()
# 
# dotheatmap %>% group_by(cluster,tumor) %>% summarise(mean = mean(cox_predictans),
#                                                      siz = length(cox_predictans)) -> 
#   tmp3
# ggplot(tmp3, aes(x = cluster, y =tumor, size = siz, color = mean)) + geom_point() +
#   scale_color_gradient2(low = "#6495ED",
#                         mid = 'grey',
#                         high = 'red') + theme_bw()
# dev.off()

# pdf('violin.big.pdf',height = 7)
# VlnPlot(HPSC, 
#         features = c(
#           "cox_predictans",
#           'ENSMUSG00000039153.16-Runx2',
#           'ENSMUSG00000001496.15-Nkx2-1',
#           'ENSMUSG00000071552.4-Tigit'),
#         ncol = 1)
# dev.off()

vl <- data.frame(cox_predictans = HPSC$cox_predictans,
                 genetype = HPSC$kp,
                 time = HPSC$w,
                 tumor = HPSC$tumor,
                 cluster = HPSC$SCT_snn_res.0.5)
# 
# ggplot(vl, aes(x = cluster, y = cox_predictans, fill = genetype)) + 
#   geom_violin() + theme_bw()
library(reshape2)
library(ComplexHeatmap)
library(scales)
vl %>% group_by(cluster, time) %>% summarise(mean_pre = mean(cox_predictans)) -> cti
# vl %>% group_by(cluster, tumor) %>% summarise(mean_pre = mean(cox_predictans)) -> ctu
vl %>% group_by(cluster, genetype) %>% summarise(mean_pre = mean(cox_predictans)) -> ctg

vl %>% group_by(cluster, time) %>% count() -> ctc
vl %>% group_by(cluster, genetype) %>% count() -> ctgc


reshape2::dcast(cti, time~cluster) -> tmp
rownames(tmp) <- tmp[,1]
tmp <- tmp[,2:dim(tmp)[2]]

dcast(ctc, time~cluster) -> tmp_ctc
tmp_ctc[is.na(tmp_ctc)] <- 0
rownames(tmp_ctc) <- tmp_ctc[,1]
tmp_ctc <- tmp_ctc[,2:dim(tmp_ctc)[2]]
tmp_ctc <- apply(tmp_ctc, 2, function(x) {x/sum(x)})
tmp[tmp_ctc <= 0.026] <- NA
tmp <- tmp[c("0w", "2w", "4w", "12w", "18w", "20w", "30w"),]
# pdf('time_cluster_cox.pdf',7,4)
# Heatmap(tmp, na_col = '#6495ED',
#         cluster_columns = F,
#         cluster_rows = F)
# dev.off()
# 
# dcast(ctu, tumor~cluster) -> tmp
# rownames(tmp) <- tmp[,1]
# tmp <- tmp[,2:dim(tmp)[2]]
# pdf('tumor_cluster_cox.pdf',7,4)
# Heatmap(tmp)
# dev.off()

dcast(ctg, genetype~cluster) -> tmp
rownames(tmp) <- tmp[,1]
tmp <- tmp[,2:dim(tmp)[2]]

dcast(ctgc, genetype~cluster) -> tmp_ctgc
tmp_ctgc[is.na(tmp_ctgc)] <- 0
rownames(tmp_ctgc) <- tmp_ctgc[,1]
tmp_ctgc <- tmp_ctgc[,2:dim(tmp_ctgc)[2]]
tmp_ctgc <- apply(tmp_ctgc, 2, function(x) {x/sum(x)})
tmp[tmp_ctgc <= 0.026] <- NA

# pdf('genetype_cluster_cox.pdf',7,4)
# Heatmap(tmp, na_col = '#6495ED',
#         cluster_columns = F,
#         cluster_rows = F)
# dev.off()


data.frame(cox_predictans = HPSC$cox_predictans,
           genetype = HPSC$kp,
           time = HPSC$w) -> vp
# pdf('cox_predictans_genetype_box.pdf',width = 3,height = 4)
# ggplot(vp, aes(x = genetype, 
#                y = cox_predictans, 
#                fill = genetype)) + geom_boxplot() + theme_bw()
# dev.off()

vp$time <- factor(vp$time, levels = c("0w", "2w", "4w", "12w", "18w", "20w", "30w"))
# pdf('cox_predictans_time_box.pdf',width = 4,height = 4)
# ggplot(vp, aes(x = time, 
#                y = cox_predictans, 
#                fill = time)) + geom_boxplot() + theme_bw()
# dev.off()


library(mclust)

emobj <- Mclust(HPSC$cox_predictans,G=6)
data.frame(class = as.character(emobj$classification),
           cox_predictans = HPSC$cox_predictans) -> dt
dt$group <- 'high'
dt$group[dt$cox_predictans < -0.1] <- 'low'
# pdf('scRNA_density.pdf',5,5)
# ggplot(dt,aes(cox_predictans,group = group, fill = group)) +
#   geom_density() + theme_bw()
# 
# ggplot(dt,aes(cox_predictans, fill = 'red')) +
#   geom_density() + theme_bw()
# dev.off()

data.frame(class = HPSC$SCT_snn_res.0.5,
           cox_predictans = HPSC$cox_predictans) -> dt
dt$group <- 'high'
dt$group[dt$cox_predictans < -0.1] <- 'low'
dt %>% group_by(class,group) %>% 
  count() -> m_cox


pdf('high_low.pdf',width = 6,height = 7)
pl <- c()
for(i in 0:12){
  tmp <- m_cox[m_cox$class == i,]
  pv <- binom.test(as.numeric(tmp[1,3]),sum(tmp[,3]))
  pl <- c(pl,pv$p.value)
  
}
pd <- data.frame(class = factor(0:12,levels = 0:12),
                 pval = signif(pl,3))
m_cox %>% left_join(pd) -> m_cox
ggplot(m_cox,aes(x = class,y = n,
                 group = group,
                 shape = group,
                 color = group)) +
  geom_point(size = 3) + 
  geom_text(aes(x = class,
                y=800,
                label = pval),
            color = 'black') + 
  theme_bw() -> p1
m_cox$n[m_cox$group == 'low']/(m_cox %>% 
                                 group_by(class) %>% 
                                 summarise(sum(n)) %>% 
                                 pull(`sum(n)`)) -> pct
data.frame(
  x = unique(m_cox$class),
  High = rep(1,times = 1),
  Low = pct
) -> dt0
p2 <- ggplot(dt0) + geom_bar(
  mapping = aes(x = x, y = High,
                fill = hue_pal()(2)[2]),
  stat = 'identity',
  color = 'black') + 
  geom_bar(
    mapping = aes(x = x, y = Low,
                  fill = hue_pal()(2)[1]),
    stat = 'identity',
    color = 'black') +
  theme_bw() + geom_hline(
    yintercept = 0.5,
    color = 'red',
    linetype = 2,
    size = 1)

p1 / p2
dev.off()




dim(dt)
dim(HPSC)
dt$UMAP_1 <- HPSC@reductions$umap@cell.embeddings[,1]
dt$UMAP_2 <- HPSC@reductions$umap@cell.embeddings[,2]

# pdf('umap.pdf', width = 7,height = 7)
# ggplot(dt, aes(x = UMAP_1, y = UMAP_2, color = group)) + 
#   geom_point(alpha = 1) + theme_bw() 
# dev.off()

top100_protein <- read.table('ExoCarta_top100_protein_details_5.txt',
                             quote = '',sep = '\t',header = T)
top100_protein <- top100_protein[,1]
top100_protein <- gsub(top100_protein, pattern = ' ',replacement = '')
tensor <- read.csv('张力annurev.csv')
tensor <- unique(tensor[,1])
fe <- readRDS('fedb.Rds')
ERstressL <- read.gmt('data/genesets.gmt')
library(AUCell)
sctdata <- HPSC@assays$SCT@data
rn <- str_extract(rownames(sctdata), pattern = 'ENSMUSG[0-9]+')
sctdata <- sctdata[!is.na(rn),]
rn <- rn[!is.na(rn)]
rownames(sctdata) <- rn
cells_rankings <- AUCell_buildRankings(sctdata)  # 关键一步

human <- useMart('ensembl',dataset = "hsapiens_gene_ensembl")
mouse <- useMart('ensembl',dataset = "mmusculus_gene_ensembl")
grep(listAttributes(mouse)[,1],pattern = 'symbol',value = T)

# h2m.driver <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
#                 values = fe$driver, mart = human,
#                 attributesL = c("ensembl_gene_id","mgi_symbol"),
#                 martL = mouse,uniqueRows = T)
# h2m.suppressor <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
#                      values = fe$suppressor, mart = human,
#                      attributesL = c("ensembl_gene_id","mgi_symbol"),
#                      martL = mouse,uniqueRows = T)
# 
# fe0 <- list(driver = h2m.driver$Gene.stable.ID,
#             suppressor = h2m.suppressor$Gene.stable.ID)
# cells_AUC_fe <- AUCell_calcAUC(fe0, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
# cells_assignment <- AUCell_exploreThresholds(cells_AUC_fe,
#                                              plotHist=T, assign=TRUE)
# HPSC <- AddMetaData(HPSC,getAUC(cells_AUC_fe)[1,],'driver')
# HPSC <- AddMetaData(HPSC,getAUC(cells_AUC_fe)[2,],'suppressor')
# VlnPlot(HPSC, features = 'cox_predictans')
# VlnPlot(HPSC, features = 'driver')
# VlnPlot(HPSC, features = 'suppressor')




tmpset <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
                 values = ERstressL$gene, mart = human,
                 attributesL = c("ensembl_gene_id","mgi_symbol"),
                 martL = mouse,uniqueRows = T)
colnames(tmpset) <- c('gene','select','mgene')
ERstressL %>% left_join(tmpset) -> tmpL
tmpL <- na.omit(tmpL[,c(1,3)])

erlist <- list()
for(i in unique(tmpL$term)){
  erlist[[i]] <- unique(tmpL[tmpL$term == i,2])
}
cells_AUC_er <- AUCell_calcAUC(erlist, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
# cells_assignment <- AUCell_exploreThresholds(cells_AUC_er,
#                                              plotHist=T, assign=TRUE)
cells_AUC_er <- t(as.matrix(getAUC(cells_AUC_er)))
HPSC@meta.data <- cbind(HPSC@meta.data, cells_AUC_er)

colnames(HPSC@meta.data)



# pdf('set_violin.pdf')
# for(i in 56:81){
#   print(VlnPlot(HPSC, features = colnames(HPSC@meta.data)[i]))
# }
# VlnPlot(HPSC, 
#         features = c(
#           'ENSMUSG00000039153.16-Runx2',
#           'ENSMUSG00000001496.15-Nkx2-1',
#           'ENSMUSG00000071552.4-Tigit'),
#         ncol = 1)
# dev.off()

h2m.top100_protein <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
                             values = top100_protein, mart = human,
                             attributesL = c("ensembl_gene_id","mgi_symbol"),
                             martL = mouse,uniqueRows = T)
h2m.tensor <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
                     values = tensor, mart = human,
                     attributesL = c("ensembl_gene_id","mgi_symbol"),
                     martL = mouse,uniqueRows = T)
fe0 <- list(tensor = h2m.tensor$Gene.stable.ID,
            h2m.top100_protein = h2m.top100_protein$Gene.stable.ID)
cells_AUC_fe <- AUCell_calcAUC(fe0, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
# cells_assignment <- AUCell_exploreThresholds(cells_AUC_fe,
#                                              plotHist=T, assign=TRUE)
HPSC <- AddMetaData(HPSC,getAUC(cells_AUC_fe)[1,],'tensor')
HPSC <- AddMetaData(HPSC,getAUC(cells_AUC_fe)[2,],'exosome')

hs <- read.gmt('D:/GeneSet/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_files_to_download_locally/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt')
hlist <- list()
for(i in unique(hs$term)){
  hlist[[i]] <- unique(hs[hs$term == i,2])
}

hlistm <- list()
for(i in names(hlist)){
  hlistm[[i]] <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
                        values = hlist[[i]], mart = human,
                        attributesL = c("ensembl_gene_id","mgi_symbol"),
                        martL = mouse,uniqueRows = T)$Gene.stable.ID
}
cells_AUC_H <- AUCell_calcAUC(hlistm, cells_rankings, aucMaxRank=nrow(cells_rankings)*0.1)
hsc <- as.data.frame(t(getAUC(cells_AUC_H)))
hsc$predict <- HPSC$cox_predictans
hsc$cluster <- HPSC$seurat_clusters
hscL <- split(hsc, hsc$cluster)
hscL0 <- lapply(hscL, function(x) {apply(x, 2, function(x) mean(as.numeric(x)))})
hscL0 <- purrr::reduce(hscL0,rbind)
rownames(hscL0) <- names(hscL)

hscL0 <- as.data.frame(hscL0)
hscL0$cluster <- as.character(hscL0$cluster)
library(reshape2)
hscL0L <- melt(hscL0,c('cluster', 'predict'))
hscL0L %>% group_by(variable) %>% mutate(scale(value)) -> tmp
tmp$cluster <- factor(tmp$cluster,levels = c(0,1,2,3,4,5,6,7,8,9,10,11,12))

# ggplot(tmp, aes(x=cluster, 
#                 y=variable,
#                 size = predict,
#                 color = `scale(value)`)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_gradient2(low = '#6495ED',mid = 'grey', high = 'red')

# hsc <- t(hsc)
# hsc <- as.data.frame(hsc)
# hsc$atp1a1 <- HPSC@assays$SCT@data["ENSMUSG00000033161.10-Atp1a1",]
# 
# hscd <- data.frame(atp1a1 = hsc$atp1a1, 
#                    HALLMARK_INTERFERON_ALPHA_RESPONSE = as.numeric(hsc$HALLMARK_INTERFERON_ALPHA_RESPONSE))
# hscd <- hscd[rowSums(hscd)>0,]
# plot(hscd$atp1a1, hscd$HALLMARK_INTERFERON_ALPHA_RESPONSE)

# pdf('Hsc0.pdf',height = 10)
# ggplot(tmp, aes(x=cluster, 
#                 y=variable,
#                 size = predict,
#                 color = `scale(value)`)) + 
#   geom_point() + 
#   theme_bw() + 
#   scale_color_gradient2(low = '#6495ED',mid = 'grey', high = 'red')
# dev.off()

ERD <- HPSC@meta.data[,57:79]
ERD <- as.data.frame(ERD)
ERD <- split(ERD, HPSC$seurat_clusters)
ERD <- lapply(ERD, function(x) apply(x,2,mean))
rn <- names(ERD)
ERD <- purrr::reduce(ERD, rbind)
rownames(ERD) <- rn
ERD_od <- apply(ERD, 2, function(x) rank(x))
rn %in% c(1,9,12)
od <- apply(ERD_od[c(1,9,12),],2,sum)
od <- sort(od,decreasing = F)
plot(hist(od))

pdf('long_ER_order.pdf', height = 50)
VlnPlot(HPSC, features = c(c('exosome'),c(
  'ENSMUSG00000039153.16-Runx2',
  'ENSMUSG00000001496.15-Nkx2-1',
  'ENSMUSG00000071552.4-Tigit'),
  names(od)),
  ncol = 1)
dev.off()

pn <- names(od[c(1:6, 19:23)])

pdf('tensor.pdf')
VlnPlot(HPSC, features = c('cox_predictans',
                           'top100_protein',
                           'tensor'),
        ncol = 1)

dev.off()

pdf('point.pdf')
# FeaturePlot(HPSC, features = c('top100_protein'),
#             label = T)
# FeaturePlot(HPSC, features = c('cox_predictans'),
#             label = T)
# 
# FeaturePlot(HPSC, features = c('top100_protein'),
#             label = T, reduction = "tsne")
# FeaturePlot(HPSC, features = c('cox_predictans'),
#             label = T, reduction = "tsne")
# 
# FeaturePlot(HPSC, features = c('tensor'),
#             label = T, reduction = "umap")
# FeaturePlot(HPSC, features = c('tensor'),
#             label = T, reduction = "tsne")

for(i in pn){
  print(
    FeaturePlot(HPSC, features = i,
                label = T, reduction = "umap"))
  print(
    FeaturePlot(HPSC, features = i,
                label = T, reduction = "tsne"))
}
# DimPlot(object = HPSC, reduction = "tsne", label = T)
# DimPlot(object = HPSC, reduction = "umap", label = T)
dev.off()

# ======================== CNV ===========================
pt <- 'D:/TCGA-LUAD/肺癌数据/HPSC/scRNA/'
library(reshape2)
library(ggrepel)

cnv_table <- read.table(paste0(pt,"/infcnv/output_dir/infercnv.observations.txt"), 
                        header=T)
# Score cells based on their CNV scores 
# Replicate the table 
cnv_score_table <- as.matrix(cnv_table)
cnv_score_mat <- as.matrix(cnv_table)
# Scoring
# CNV的5分类系统
cnv_score_table[cnv_score_mat >= 0 & cnv_score_mat < 0.3] <- "A" #complete loss. 2pts
cnv_score_table[cnv_score_mat >= 0.3 & cnv_score_mat < 0.7] <- "B" #loss of one copy. 1pts
cnv_score_table[cnv_score_mat >= 0.7 & cnv_score_mat < 1.3] <- "C" #Neutral. 0pts
cnv_score_table[cnv_score_mat >= 1.3 & cnv_score_mat <= 1.5] <- "D" #addition of one copy. 1pts
cnv_score_table[cnv_score_mat > 1.5 & cnv_score_mat <= 2] <- "E" #addition of two copies. 2pts
cnv_score_table[cnv_score_mat > 2] <- "F" #addition of more than two copies. 2pts

# cnv_score_table[cnv_score_mat >= 0 & cnv_score_mat < 2.5] <- "A" 
# cnv_score_table[cnv_score_mat >= 2.5 & cnv_score_mat < 3.5] <- "B" 
# cnv_score_table[cnv_score_mat >= 3.5 & cnv_score_mat < 4.5] <- "C" 
# cnv_score_table[cnv_score_mat >= 4.5 & cnv_score_mat < 5.5] <- "D" 
# Check
table(cnv_score_table[,2])
# Replace with score 
cnv_score_table_pts <- cnv_table

#  把5分类调整为 3分类系统
cnv_score_table_pts[cnv_score_table == "A"] <- 2
cnv_score_table_pts[cnv_score_table == "B"] <- 1
cnv_score_table_pts[cnv_score_table == "C"] <- 0
cnv_score_table_pts[cnv_score_table == "D"] <- 1
cnv_score_table_pts[cnv_score_table == "E"] <- 2
cnv_score_table_pts[cnv_score_table == "F"] <- 2

cell_scores_CNV <- as.data.frame(colSums(cnv_score_table_pts))
colnames(cell_scores_CNV) <- "cnv_score"
head(cell_scores_CNV)
# write.csv(x = cell_scores_CNV, file = "cnv_scores.csv")
whole_obj <- readRDS(paste0(pt,'/HPSC_obj.Rds'))
cell_scores_CNV$rn <- rownames(cell_scores_CNV)
f <- rownames(whole_obj@meta.data) %in% rownames(cell_scores_CNV)
dcnv <- 
  data.frame(
    cluster = whole_obj@meta.data$seurat_clusters,
    genetype = whole_obj@meta.data$kp,
    time = whole_obj@meta.data$w,
    tumor = whole_obj@meta.data$tumor)
rownames(dcnv) <- rownames(whole_obj@meta.data)
dcnv <- dcnv[f,]
dcnv$rn <- rownames(dcnv)
dcnv %>% left_join(cell_scores_CNV) -> dcnv

L <- melt(cnv_score_mat)
dcnv$time <- factor(dcnv$time, levels = c('0w','2w','4w','12w','18w','20w','30w'))
# pdf('cnvsupp.pdf')
# ggplot(L, aes(x = value,fill = 'red')) + geom_density() + theme_bw()
# ggplot(dcnv, aes(x = cluster, y = cnv_score, fill = cluster)) + 
#   geom_boxplot() + theme_bw()
# ggplot(dcnv, aes(x = genetype, y = cnv_score, fill = genetype)) + 
#   geom_boxplot() + theme_bw()
# ggplot(dcnv, aes(x = time, y = cnv_score, fill = time)) + 
#   geom_boxplot() + theme_bw()
# ggplot(dcnv, aes(x = tumor, y = cnv_score, fill = tumor)) + 
#   geom_boxplot() + theme_bw()
# 
# ggplot(dcnv, aes(x = cnv_score, fill = cluster)) + 
#   geom_density(alpha = 0.3) + theme_bw()
# ggplot(dcnv, aes(x = cnv_score, fill = genetype)) + 
#   geom_density(alpha = 0.3) + theme_bw()
# ggplot(dcnv, aes(x = cnv_score, fill = time)) + 
#   geom_density(alpha = 0.3) + theme_bw()
# ggplot(dcnv, aes(x = cnv_score, fill = tumor)) + 
#   geom_density(alpha = 0.3) + theme_bw()
# dev.off()

predictans <- readRDS('predictans_ans.Rds')
whole_obj <- AddMetaData(whole_obj,predictans,'cox_predictans')

dcnv %>% group_by(cluster) %>% 
  summarise(cnv_score = mean(cnv_score)) -> cnv_score1
data.frame(cox_predictans = whole_obj@meta.data$cox_predictans,
           group = whole_obj$seurat_clusters) %>% group_by(group) %>% 
  summarise(mean(cox_predictans)) -> atp1a12
atp1a12 <- atp1a12[c(1:11,13),]
cbind(cnv_score = cnv_score1,
      mean_atp1a1 = atp1a12) -> cnvatp1a1
cnvatp1a1$cnv_score.cluster <- factor(cnvatp1a1$cnv_score.cluster)
# pdf('mean_cox_predictans_cnv.pdf',width = 5,height = 4)
# ggplot(cnvatp1a1,aes(x = cnv_score.cnv_score,
#                      y = `mean_atp1a1.mean(cox_predictans)`,
#                      color = cnv_score.cluster)) +
#   geom_text_repel(mapping = aes(x = cnv_score.cnv_score,
#                                 y = `mean_atp1a1.mean(cox_predictans)`,
#                                 label = cnv_score.cluster),
#                                          color = 'black') +
#   geom_smooth(cnvatp1a1[cnvatp1a1$cnv_score.cluster != 6,],
#               mapping =aes(x = cnv_score.cnv_score,
#                            y = `mean_atp1a1.mean(cox_predictans)`),
#               method = 'lm', 
#               formula = y ~ x,
#               color = 'red',
#               se = F,
#               linetype = 2,
#               size = 2) +
#   geom_point(size = 6) + 
#   theme_bw() + xlab('CNV_Score') + ylab('mean_cox_predictans')
# dev.off()

cor(cnvatp1a1[cnvatp1a1$cnv_score.cluster != 6,]$cnv_score.cnv_score,
    cnvatp1a1[cnvatp1a1$cnv_score.cluster != 6,]$`mean_atp1a1.mean(cox_predictans)`,
    method = 'pearson')
cor.test(cnvatp1a1[cnvatp1a1$cnv_score.cluster != 6,]$cnv_score.cnv_score,
         cnvatp1a1[cnvatp1a1$cnv_score.cluster != 6,]$`mean_atp1a1.mean(cox_predictans)`,
         method = 'pearson')


FeaturePlot(HPSC,
            features = c(
              'ENSMUSG00000039153.16-Runx2',
              'ENSMUSG00000001496.15-Nkx2-1',
              'ENSMUSG00000071552.4-Tigit'))
# ===================== CNV2 =============================
options(stringsAsFactors = F)
library(phylogram)
library(gridExtra)
library(grid)
require(dendextend)
require(ggthemes)
library(tidyverse)
library(Seurat)
library(infercnv)
library(miscTools)
whole_obj <- readRDS(paste0(pt,'/HPSC_obj.Rds'))
#  Import inferCNV dendrogram
infercnv.dend <- read.dendrogram(file = paste0(pt,"infcnv/output_dir_new/infercnv.observations_dendrogram.txt"))
# Cut tree 
infercnv.labels <- cutree(infercnv.dend, k = 9, order_clusters_as_data = FALSE)

nsd <- 
  data.frame(
    lable = infercnv.labels,
    rn = names(infercnv.labels))
dcnv <- 
  data.frame(
    rn = rownames(whole_obj@meta.data),
    cluster = whole_obj@meta.data$seurat_clusters,
    genetype = whole_obj@meta.data$kp,
    time = whole_obj@meta.data$w,
    tumor = whole_obj@meta.data$tumor)
nsd %>% left_join(dcnv) -> nsd
nsd$lable <- factor(nsd$lable)
table(nsd$lable, nsd$cluster)
table(nsd$lable, nsd$time)
table(nsd$lable, nsd$tumor)
table(nsd$lable, nsd$genetype)

pdf('cnv_cluster.pdf')
ggplot(nsd, aes(x = lable, fill = cluster)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = lable, fill = time)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = lable, fill = tumor)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = lable, fill = genetype)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()

ggplot(nsd, aes(x = genetype, fill = lable)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = tumor, fill = lable)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
ggplot(nsd, aes(x = cluster, fill = lable)) + 
  geom_bar(position = 'fill', color = 'black') + theme_bw()
dev.off()

# Color labels
the_bars <- as.data.frame(tableau_color_pal("Tableau 20")(20)[infercnv.labels])
colnames(the_bars) <- "inferCNV_tree"
the_bars$inferCNV_tree <- as.character(the_bars$inferCNV_tree)

pdf('cnv_tree.pdf',height = 10)
infercnv.dend %>% set("labels",rep("", 
                                   nobs(infercnv.dend)))  %>% 
  plot(main="inferCNV dendrogram") %>%
  colored_bars(colors = as.data.frame(the_bars), 
               dend = infercnv.dend, 
               sort_by_labels_order = F, 
               add = T, y_scale=100 , y_shift = 0)
dev.off()

cd <- data.frame(cnv = the_bars$inferCNV_tree,
                 lables = infercnv.labels)
cd$lables <- factor(cd$lables, levels = names(table(cd$lables)))
v <- names(table(cd$cnv))
names(v) <- apply(table(cd$cnv, cd$lables), 1, function(x) which(x != 0))
pdf('color_bar.pdf')
ggplot(cd, aes(x = lables, fill = lables)) + 
  geom_bar() + 
  scale_fill_manual(values = v) +
  theme_bw()
dev.off()





# ===================== moncole2 ========
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle)
pt <- 'D:/TCGA-LUAD/肺癌数据/HPSC/scRNA/'
predictans <- readRDS('predictans_ans.Rds')

# # ---------------------------- 所有基因 -----------------------
# # 凑活太碎了
# deg <- read.csv(paste0(pt,"expressed_genes_train.monocle.DEG.csv"))
# cds <- readRDS(paste0(pt,'expressed_genes_monocle.Rds'))
# # ---------------------------- seurat高变基因基因 -----------------------
# # 还行
# deg <- read.csv(paste0(pt,"seurat_express_genes_train.monocle.DEG.csv"))
# cds <- readRDS(paste0(pt,'seurat_express_genes_monocle.Rds'))
# # ---------------------------- moncole高变基因基因 -----------------------
# # 还行
# deg <- read.csv(paste0(pt,"moncole_genes_train.monocle.DEG.csv"))
# cds <- readRDS(paste0(pt,'moncole_genes_monocle.Rds'))
# moncole_highval <- as.matrix(cds@reducedDimS)
# write.csv(moncole_highval,file = 'moncole_highval.csv')

# plot_complex_cell_trajectory(cds,root_state = 2) -> x
# x <- data.frame(dim1 = x$data$data_dim_1, dim2 = x$data$data_dim_2)
# write.csv(x, file = 'moncole_minspan_highval.csv')
# barcode0highval <- colnames(cds)
# write.csv(barcode0highval, file = 'barcode0highval.csv')
# # ---------------------------- cluster差异基因 -----------------------
# # 不好
# deg <- read.csv(paste0(pt,"cluster_express_genes_train.monocle.DEG.csv"))
# cds <- readRDS(paste0(pt,'cluster_express_genes_monocle.Rds'))
# # ---------------------------- stemness genes -----------------------
# # 最好
# deg <- read.csv(paste0(pt,"stemness_genes_train.monocle.DEG.csv"))
# cds <- readRDS(paste0(pt,'stemness_genes_monocle.Rds'))
# moncole_stemness <- as.matrix(cds@reducedDimS)
# write.csv(moncole_stemness,file = 'moncole_stemness.csv')
# plot_complex_cell_trajectory(cds,
#                              root_state = 6) -> x
# x <- data.frame(dim1 = x$data$data_dim_1, dim2 = x$data$data_dim_2)
# write.csv(x, file = 'moncole_minspan_stemness.csv')
# barcode0stemness <- colnames(cds)
# write.csv(barcode0stemness, file = 'barcode0stemness.csv')
# # ======================= monocle analysis ==================================
# library(dplyr)
# library(Seurat)
# library(patchwork)
# library(monocle)
# library(AUCell)
# library(biomaRt)
# library(clusterProfiler)
# 
# stemness <- read.gmt(paste0(pt,'stemness.gmt'))
# # human <- readRDS(paste0(pt,'humandb.Rds'))
# # mouse <- readRDS(paste0(pt,'mousedb.Rds'))
# # 
# # stemnessL <- list()
# # for(i in unique(stemness$term)){
# #   .driver <- getLDS(attributes = c("hgnc_symbol"),filters = "hgnc_symbol",
# #                     values = stemness[stemness$term == i,2], mart = human,
# #                     attributesL = c("ensembl_gene_id","mgi_symbol"),
# #                     martL = mouse,uniqueRows = T)$Gene.stable.ID
# #   stemnessL[[i]] <- unique(.driver)
# # }
# # saveRDS(stemnessL, file = 'stemness_mouse.Rds')
# stemnessL <- readRDS('stemness_mouse.Rds')
# # ----------------------------- 所有基因 ---------------------
# deg <- read.csv(paste0(pt,"expressed_genes_train.monocle.DEG.csv"))
# cds <- readRDS(paste0(pt,'expressed_genes_monocle.Rds'))
# ordergene <- deg$X
# 
# exp_data <- exprs(cds)
# rn <- gsub(rownames(exp_data), pattern = '\\..*', replacement = '')
# rownames(exp_data) <- rn
# cells_rankings <- AUCell_buildRankings(exp_data)
# cells_AUC <- list()
# for(i in unique(stemness$term)){
#   cells_AUC[[i]] <- getAUC(AUCell_calcAUC(stemnessL[[i]], cells_rankings,
#                                           aucMaxRank=nrow(cells_rankings)*0.1))
# }
# ns <- names(cells_AUC)
# cells_AUC <- purrr::reduce(cells_AUC, rbind)
# cells_AUC <- 10^(cells_AUC - 0.1)
# rownames(cells_AUC) <- ns
# cells_AUC <- t(cells_AUC)
# cells_AUC <- as.data.frame(cells_AUC)
# predictans <- predictans[abs(predictans) < 1]
# predictans0 <- 10^(scale(predictans[colnames(cds)]) - 0.1)
# predictans_high_low <- rep('High', times = length(predictans0))
# predictans_high_low[predictans0 < (10^(-0.1)-0.1)] <- 'Low'
# 
# pData(cds) <- cbind(pData(cds), cells_AUC,predictans0,predictans_high_low)
# 
# cds <- orderCells(cds, root_state = 3)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# pdf('total_moncole.pdf',7,7)
# cds <- orderCells(cds, root_state = 3)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# Heatmap(small_mat,
#         col = circlize::colorRamp2(c(0,2.5), c("#6495ED","red")),
#         layer_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.f", pindex(table(cds$State,cds$seurat_clusters), i, j)), x, y, gp = gpar(fontsize = 10))
#         }
#         )
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans0",
#                              root_state = 3)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans_high_low",
#                              root_state = 3)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "seurat_clusters",
#                              root_state = 3)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "BHATTACHARYA_EMBRYONIC_STEM_CELL",
#                              root_state = 3)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "Pseudotime",
#                              root_state = 3)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "State",
#                              root_state = 3)
# 
# 
# cds <- orderCells(cds, root_state = 5)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans0",
#                              root_state = 5)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans_high_low",
#                              root_state = 5)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "seurat_clusters",
#                              root_state = 5)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "BHATTACHARYA_EMBRYONIC_STEM_CELL",
#                              root_state = 5)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "Pseudotime",
#                              root_state = 5)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "State",
#                              root_state = 5)
# dev.off()
# # ----------------------------- seurat高变基因基因 ---------------------
# deg <- read.csv(paste0(pt,"seurat_express_genes_train.monocle.DEG.csv"))
# cds <- readRDS(paste0(pt,'seurat_express_genes_monocle.Rds'))
# ordergene <- deg$X
# 
# exp_data <- exprs(cds)
# rn <- gsub(rownames(exp_data), pattern = '\\..*', replacement = '')
# rownames(exp_data) <- rn
# cells_rankings <- AUCell_buildRankings(exp_data)
# cells_AUC <- list()
# for(i in unique(stemness$term)){
#   cells_AUC[[i]] <- getAUC(AUCell_calcAUC(stemnessL[[i]], cells_rankings,
#                                           aucMaxRank=nrow(cells_rankings)*0.1))
# }
# ns <- names(cells_AUC)
# cells_AUC <- purrr::reduce(cells_AUC, rbind)
# cells_AUC <- 10^(cells_AUC - 0.1)
# rownames(cells_AUC) <- ns
# cells_AUC <- t(cells_AUC)
# cells_AUC <- as.data.frame(cells_AUC)
# predictans <- predictans[abs(predictans) < 1]
# predictans0 <- 10^(scale(predictans[colnames(cds)]) - 0.1)
# predictans_high_low <- rep('High', times = length(predictans0))
# predictans_high_low[predictans0 < (10^(-0.1)-0.1)] <- 'Low'
# 
# pData(cds) <- cbind(pData(cds), cells_AUC,predictans0,predictans_high_low)
# 
# cds <- orderCells(cds)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# pdf('seurat_high_val_moncole.pdf',5,5)
# cds <- orderCells(cds, root_state = 1)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# Heatmap(small_mat,
#         col = circlize::colorRamp2(c(0,2.5), c("#6495ED","red")),
#         layer_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.f", pindex(table(cds$State,cds$seurat_clusters), i, j)), x, y, gp = gpar(fontsize = 10))
#         }
# )
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans0",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans_high_low",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "seurat_clusters",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "BHATTACHARYA_EMBRYONIC_STEM_CELL",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "Pseudotime",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "State",
#                              root_state = 1)
# 
# dev.off()
# ----------------------------- moncole高变基因基因 ---------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(AUCell)
library(ComplexHeatmap)
library(ggplot2)
library(scatterpie)
pt <- 'D:/TCGA-LUAD/肺癌数据/HPSC/scRNA/'

predictans <- readRDS('predictans_ans.Rds')
stemness <- read.gmt(paste0(pt,'stemness.gmt'))
stemnessL <- readRDS('stemness_mouse.Rds')
deg <- read.csv(paste0(pt,"moncole_genes_train.monocle.DEG.csv"))
cds <- readRDS(paste0('moncole_genes_monocle.Rds'))
cds <- orderCells(cds, root_state = 2)
cell_scores_CNV <- read.csv("cnv_scores.csv")
dir.create('./moncole_high_val_moncole')
# data.frame(predictans = cds$raw_predictans, State = cds$State) %>%
#   ggplot(aes(x = State, y = predictans, fill = State)) +
#   geom_violin(trim = F) +
#   geom_jitter(width = 0.3,size = 0.1) +
#   theme_bw() -> t
# ggsave(t,filename = './moncole_high_val_moncole/violin.pdf',height = 6,width = 7)
# 
data.frame(predictans = cds$raw_predictans, State = cds$State,
           risk_level = cds$predictans_high_low) %>%
  group_by(State,risk_level) %>%
  count -> td
plnb <- list()
for(i in 1:7){
  tmp <- td[td$State == i,]
  binom.test(as.numeric(tmp[1,3]),as.numeric(tmp[1,3]+tmp[2,3])) -> tmp0
  plnb[[i]] <- tmp0$p.value
}
bind <- data.frame(State = factor(1:7,levels = 1:7),p = signif(unlist(plnb),3))
td %>% left_join(bind) -> td0

ggplot(td0,aes(x = State, y = n,
               group = risk_level, color = risk_level,
               shape = risk_level)) +
  geom_point(size = 3,alpha=0.9) +
  theme_bw()+ geom_text(aes(x = State, 
                            y = 1000, 
                            label = p),
                        color = 'black',size = 4) -> t
t
td0
(td0$n[td0$risk_level == 'Low'])/
  (td0 %>% group_by(State) %>% summarise(sum(n)) %>% pull(`sum(n)`)) -> 
  pct
data <- data.frame(State = factor(1:7, levels = 1:7),
                   High = rep(1,times = 7),
                   Low = pct
)
library(scales)
ggplot()+ geom_bar(data, mapping = aes(x = State, 
                                       y = High),
                   stat = 'identity',
                   color = 'black',
                   fill = hue_pal()(2)[1]) + 
  geom_bar(data, mapping = aes(x = State, 
                               y = Low),
           color = 'black',
           stat = 'identity',
           fill = hue_pal()(2)[2]) + 
  geom_hline(
    yintercept = 0.5,
    color = 'red',
    linetype = 2,
    size = 1) +
  theme_bw() -> 
  t0
t/t0 -> t
ggsave(t,
       filename = './moncole_high_val_moncole/point_path.pdf',
       height = 6,width = 7)

# ordergene <- deg$X
# exp_data <- exprs(cds)
# rn <- gsub(rownames(exp_data), pattern = '\\..*', replacement = '')
# rownames(exp_data) <- rn
# cells_rankings <- AUCell_buildRankings(exp_data)
# cells_AUC <- list()
# for(i in unique(stemness$term)){
#   cells_AUC[[i]] <- getAUC(AUCell_calcAUC(stemnessL[[i]], cells_rankings,
#                                           aucMaxRank=nrow(cells_rankings)*0.1))
# }
# ns <- names(cells_AUC)
# cells_AUC <- purrr::reduce(cells_AUC, rbind)
# cells_AUC <- 10^(cells_AUC - 0.1)
# rownames(cells_AUC) <- ns
# cells_AUC <- t(cells_AUC)
# cells_AUC <- as.data.frame(cells_AUC)
# predictans <- predictans[abs(predictans) < 1]
# predictans0 <- 10^(scale(predictans[colnames(cds)]) - 0.1)
# predictans_high_low <- rep('High', times = length(predictans0))
# predictans_high_low[predictans0 < (10^(-0.1)-0.1)] <- 'Low'
# raw_predictans <- predictans[cds$sampleID]
# pData(cds) <- cbind(pData(cds), cells_AUC,predictans0,
#                     predictans_high_low,raw_predictans)
# 
# cds <- orderCells(cds)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# saveRDS(cds, file = 'moncole_genes_monocle.Rds')
# diff_test_res <- differentialGeneTest(cds,
#                                       fullModelFormulaStr = "~sm.ns(Pseudotime)")
# write.csv(diff_test_res, file = 'monocle_high_val_diff_test_res.csv')
diff_test_res <- read.csv('monocle_high_val_diff_test_res.csv')
# ===== ~~~~~~ 相关基因 ~~~~~~~~~~~~~~~~~ =====
diff_test_res <- diff_test_res[order(diff_test_res$qval,decreasing = F),]
memory.limit(62405)

sig_gene_names <- diff_test_res$X[diff_test_res$qval < 1e-20]
# length(sig_gene_names)
# trend_formula = "~sm.ns(Pseudotime, df=3)"
# cds_subset <- cds[sig_gene_names,]
# pseudocount <- 1
# newdata <- data.frame(Pseudotime = seq(
#   min(pData(cds_subset)$Pseudotime),
#   max(pData(cds_subset)$Pseudotime), length.out = 3891))
# m <- genSmoothCurves(cds_subset, cores = 4,
#                      trend_formula = trend_formula,
#                      relative_expr = T, new_data = newdata)
# m = m[!apply(m, 1, sum) == 0, ]
# m = vstExprs(cds_subset, expr_matrix = m)
# m = m[!apply(m, 1, sd) == 0, ]
# m = Matrix::t(scale(Matrix::t(m), center = TRUE))
# m = m[is.na(row.names(m)) == FALSE, ]
# m[is.nan(m)] = 0
# 
# scale_max = 3
# scale_min = -3
# m[m > scale_max] = scale_max
# m[m < scale_min] = scale_min
# heatmap_matrix <- m
# 
# rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix),
#                                  pattern = 'ENSMUSG[-.0-9]+',
#                                  replacement = '')
# name <- gsub(sig_gene_names[1:10], pattern = 'ENSMUSG[-.0-9]+',
#              replacement = '')

# pdf('./moncole_high_val_moncole/high_val_moncole_all.pdf',height = 7,width = 7)
# Heatmap(heatmap_matrix, cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         clustering_distance_rows = "pearson",
#         col = circlize::colorRamp2(c(-2,-1,0,1,2),
#                                    c("blue",
#                                      'Cyan',
#                                      "green",
#                                      'yellow',
#                                      "red")))+
#   rowAnnotation(link = anno_mark(at = which(rownames(
#     heatmap_matrix) %in% name),
#     labels = name,
#     labels_gp = gpar(fontsize = 15)))
# dev.off()

# dir.create('moncole_high_val_moncole/diff_test_res')
# for(i in 1:10){
#   fn <- paste0('moncole_high_val_moncole/diff_test_res/',sig_gene_names[i],'.pdf')
# 
#   pps <- data.frame(predictans = cds$raw_predictans,
#                     Pseudotime = cds$Pseudotime,
#                     gene = exprs(cds)[sig_gene_names[i],],
#                     stat = cds$State)
#   pps %>% na.omit() %>% group_by(stat) %>%
#     summarise(mean_predictans = mean(predictans),
#               mean_gene = mean(gene),
#               mean_time = mean(Pseudotime)) %>%
#     ggplot(aes(x = mean_predictans,
#                y = mean_time,
#                color = stat,
#                size = mean_gene)) +
#     geom_point() + theme_bw() +
#     geom_text_repel(mapping = aes(x = mean_predictans,
#                                   y = mean_time,
#                                   label = stat),
#                     size = 4,color = 'black') -> t
#   ggsave(t, filename =  fn)
# 
# }

high_val_total <- list()
for(i in 1:10){
  pps <- data.frame(predictans = cds$raw_predictans,
                    Pseudotime = cds$Pseudotime,
                    gene = exprs(cds)[sig_gene_names[i],],
                    stat = cds$State)
  pps %>% na.omit() %>% group_by(stat) %>%
    summarise(mean_predictans = mean(predictans),
              mean_gene = mean(gene),
              mean_time = mean(Pseudotime),
              gene_name = sig_gene_names[i],
              group = 'total',
              col_order = i) -> high_val_total[[i]]
}
high_val_total <- purrr::reduce(high_val_total, rbind)


# BEAM_res1 <- BEAM(cds, branch_point = 1, cores = 4)
# BEAM_res2 <- BEAM(cds, branch_point = 2, cores = 4)
# BEAM_res3 <- BEAM(cds, branch_point = 3, cores = 4)
# 
# list(BEAM_res1 = BEAM_res1,
#      BEAM_res2 = BEAM_res2,
#      BEAM_res3 = BEAM_res3) -> BEAM_moncole_high_val
# saveRDS(BEAM_moncole_high_val, file = 'BEAM_moncole_high_val.Rds')
BEAM_moncole_high_val <- readRDS('BEAM_moncole_high_val.Rds')
# ==================== ~~~~~~~~~~~~~~~ 分叉基因 ~~~~~~~~~ =====================

BEAM_res1t <- BEAM_moncole_high_val[[1]]
BEAM_res1t <- BEAM_res1t[order(BEAM_res1t$qval),]
BEAM_res1t <- BEAM_res1t[,c("gene_short_name", "pval", "qval")]
BEAM_res1t <- BEAM_res1t[order(BEAM_res1t$qval),]
rnf <- row.names(BEAM_res1t)[1:100]
dir.create('moncole_high_val_moncole/branched1')

# 分叉1小
branch_labels = c("Cell fate 1", "Cell fate 2")
cds_subset <- cds[row.names(subset(BEAM_res1t,
                                   qval < 1e-5)),]
new_cds <- buildBranchCellDataSet(cds_subset,
                                  branch_point = 1,
                                  progenitor_method = "duplicate")
new_cds@dispFitInfo <- cds_subset@dispFitInfo
progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)

col_gap_ind <- 101
newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
trend_formula <- "~sm.ns(Pseudotime, df=3) * Branch"
BranchAB_exprs <- genSmoothCurves(new_cds[,], cores = 4,
                                  trend_formula = trend_formula,
                                  relative_expr = T, new_data = rbind(newdataA,
                                                                      newdataB))
BranchA_exprs <- BranchAB_exprs[, 1:100]
BranchB_exprs <- BranchAB_exprs[, 101:200]
common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
                                                    setdiff(pData(new_cds)$State,
                                                            branch_states), ])
BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells,
                                               "Pseudotime"])))
BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells,
                                        "Pseudotime"]))
BranchB_num <- BranchA_num
BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1],
                        BranchB_exprs)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1,
                                       sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix),
                                 center = TRUE))
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
                                  FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0

scale_max = 3
scale_min = -3
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,
                                                          1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                             `Cell Type` = c(rep(branch_labels[1], BranchA_num),
                                             rep("Pre-branch", 2 * BranchP_num),
                                             rep(branch_labels[2],
                                                 BranchB_num)))
branch_colors = c("#979797", "#F05662", "#7990C8")
csplit <- c(rep('Cell fate 1', times = 100),
            rep('Cell fate 2', times = 100))
rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix),
                                 pattern = 'ENSMUSG[-.0-9]+',
                                 replacement = '')
name <- rownames(heatmap_matrix)[1:10]
library(ComplexHeatmap)
ha = HeatmapAnnotation(
  `Cell Type` = annotation_col$Cell.Type,
  col = list(`Cell Type` = c("Cell fate 1" = branch_colors[2],
                             "Cell fate 2" = branch_colors[3],
                             "Pre-branch" = branch_colors[1])),
  na_col = "grey")


pdf('moncole_high_val_moncole/branched1/branched1_top100.pdf',
    height = 6,width = 7)
Heatmap(heatmap_matrix, cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        column_split = csplit,
        top_annotation = ha,
        clustering_distance_rows = "pearson",
        col = circlize::colorRamp2(c(-2,-1,0,1,2),
                                   c("blue",'Cyan',
                                     "green",
                                     'yellow',
                                     "red")))+
  rowAnnotation(link = anno_mark(at = which(rownames(heatmap_matrix) %in% name),
                                 labels = name,
                                 labels_gp = gpar(fontsize = 10)))
dev.off()
library(ggrepel)
# 气泡图top10
# for(i in 1:10){
#   fn = paste0('moncole_high_val_moncole/branched1/', row.names(BEAM_res1t)[i],'.pdf')
#   pps <- data.frame(predictans = cds$raw_predictans,
#                     Pseudotime = cds$Pseudotime,
#                     gene = exprs(cds)[row.names(BEAM_res1t)[i],],
#                     stat = cds$State)
#   pps %>% na.omit() %>% group_by(stat) %>%
#     summarise(mean_predictans = mean(predictans),
#               mean_gene = mean(gene),
#               mean_time = mean(Pseudotime))  %>%
#     ggplot(aes(x = mean_predictans,
#                y = mean_time,
#                color = stat,
#                size = mean_gene)) +
#     geom_point() + theme_bw() +
#     geom_text_repel(mapping = aes(x = mean_predictans,
#                                   y = mean_time,
#                                   label = stat),
#                     size = 4,color = 'black') -> t
#   ggsave(t,filename = fn)
# }



high_val_branch1 <- list()
for(i in 1:10){
  pps <- data.frame(predictans = cds$raw_predictans,
                    Pseudotime = cds$Pseudotime,
                    gene = exprs(cds)[row.names(BEAM_res1t)[i],],
                    stat = cds$State)
  pps %>% na.omit() %>% group_by(stat) %>%
    summarise(mean_predictans = mean(predictans),
              mean_gene = mean(gene),
              mean_time = mean(Pseudotime),
              gene_name = row.names(BEAM_res1t)[i],
              group = 'branch1',
              col_order = i) -> high_val_branch1[[i]]
}
high_val_branch1 <- purrr::reduce(high_val_branch1, rbind)





BEAM_res2t <- BEAM_moncole_high_val[[2]]
BEAM_res2t <- BEAM_res2t[order(BEAM_res2t$qval),]
BEAM_res2t <- BEAM_res2t[,c("gene_short_name", "pval", "qval")]
BEAM_res2t <- BEAM_res2t[order(BEAM_res2t$qval),]
rnf <- row.names(BEAM_res2t)[1:100]
dir.create('moncole_high_val_moncole/branched2')

# 分叉2小
branch_labels = c("Cell fate 1", "Cell fate 2")
cds_subset <- cds[row.names(subset(BEAM_res2t,
                                   qval < 1e-5)),]
new_cds <- buildBranchCellDataSet(cds_subset,
                                  branch_point = 2,
                                  progenitor_method = "duplicate")
new_cds@dispFitInfo <- cds_subset@dispFitInfo
progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)

col_gap_ind <- 101
newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
trend_formula <- "~sm.ns(Pseudotime, df=3) * Branch"
BranchAB_exprs <- genSmoothCurves(new_cds[,], cores = 4,
                                  trend_formula = trend_formula,
                                  relative_expr = T, new_data = rbind(newdataA,
                                                                      newdataB))
BranchA_exprs <- BranchAB_exprs[, 1:100]
BranchB_exprs <- BranchAB_exprs[, 101:200]
common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
                                                    setdiff(pData(new_cds)$State,
                                                            branch_states), ])
BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells,
                                               "Pseudotime"])))
BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells,
                                        "Pseudotime"]))
BranchB_num <- BranchA_num
BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1],
                        BranchB_exprs)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1,
                                       sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix),
                                 center = TRUE))
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
                                  FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0

scale_max = 3
scale_min = -3
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,
                                                          1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                             `Cell Type` = c(rep(branch_labels[1], BranchA_num),
                                             rep("Pre-branch", 2 * BranchP_num),
                                             rep(branch_labels[2],
                                                 BranchB_num)))
branch_colors = c("#979797", "#F05662", "#7990C8")
csplit <- c(rep('Cell fate 1', times = 100),
            rep('Cell fate 2', times = 100))
rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix),
                                 pattern = 'ENSMUSG[-.0-9]+',
                                 replacement = '')
name <- rownames(heatmap_matrix)[1:10]
library(ComplexHeatmap)
ha = HeatmapAnnotation(
  `Cell Type` = annotation_col$Cell.Type,
  col = list(`Cell Type` = c("Cell fate 1" = branch_colors[2],
                             "Cell fate 2" = branch_colors[3],
                             "Pre-branch" = branch_colors[1])),
  na_col = "grey")
pdf('moncole_high_val_moncole/branched2/branched2_top100.pdf',height = 6,width = 7)
Heatmap(heatmap_matrix, cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        column_split = csplit,
        top_annotation = ha,
        clustering_distance_rows = "pearson",
        col = circlize::colorRamp2(c(-2,-1,0,1,2),
                                   c("blue",'Cyan',
                                     "green",
                                     'yellow',
                                     "red")))+
  rowAnnotation(link = anno_mark(at = which(rownames(heatmap_matrix) %in% name),
                                 labels = name,
                                 labels_gp = gpar(fontsize = 10)))
dev.off()

# 气泡图top10
# for(i in 1:10){
#   fn = paste0('moncole_high_val_moncole/branched2/', row.names(BEAM_res2t)[i],'.pdf')
#   pps <- data.frame(predictans = cds$raw_predictans,
#                     Pseudotime = cds$Pseudotime,
#                     gene = exprs(cds)[row.names(BEAM_res2t)[i],],
#                     stat = cds$State)
#   pps %>% na.omit() %>% group_by(stat) %>%
#     summarise(mean_predictans = mean(predictans),
#               mean_gene = mean(gene),
#               mean_time = mean(Pseudotime))  %>%
#     ggplot(aes(x = mean_predictans,
#                y = mean_time,
#                color = stat,
#                size = mean_gene)) +
#     geom_point() + theme_bw() +
#     geom_text_repel(mapping = aes(x = mean_predictans,
#                                   y = mean_time,
#                                   label = stat),
#                     size = 4,color = 'black') -> t
#   ggsave(t,filename = fn)
# }


high_val_branch2 <- list()
for(i in 1:10){
  pps <- data.frame(predictans = cds$raw_predictans,
                    Pseudotime = cds$Pseudotime,
                    gene = exprs(cds)[row.names(BEAM_res2t)[i],],
                    stat = cds$State)
  pps %>% na.omit() %>% group_by(stat) %>%
    summarise(mean_predictans = mean(predictans),
              mean_gene = mean(gene),
              mean_time = mean(Pseudotime),
              gene_name = row.names(BEAM_res2t)[i],
              group = 'branch2',
              col_order = i) -> high_val_branch2[[i]]
}
high_val_branch2 <- purrr::reduce(high_val_branch2, rbind)







BEAM_res3t <- BEAM_moncole_high_val[[3]]
BEAM_res3t <- BEAM_res3t[order(BEAM_res3t$qval),]
BEAM_res3t <- BEAM_res3t[,c("gene_short_name", "pval", "qval")]
BEAM_res3t <- BEAM_res3t[order(BEAM_res3t$qval),]
rnf <- row.names(BEAM_res3t)[1:100]
dir.create('moncole_high_val_moncole/branched3')

# 分叉3小
# branch_labels = c("Cell fate 1", "Cell fate 2")
# cds_subset <- cds[row.names(subset(BEAM_res3t,
#                                    qval < 1e-5)),]
# new_cds <- buildBranchCellDataSet(cds_subset,
#                                   branch_point = 3,
#                                   progenitor_method = "duplicate")
# new_cds@dispFitInfo <- cds_subset@dispFitInfo
# progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
# branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)
# 
# col_gap_ind <- 101
# newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
#                        Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
# newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
#                        Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
# trend_formula <- "~sm.ns(Pseudotime, df=3) * Branch"
# BranchAB_exprs <- genSmoothCurves(new_cds[,], cores = 4,
#                                   trend_formula = trend_formula,
#                                   relative_expr = T, new_data = rbind(newdataA,
#                                                                       newdataB))
# BranchA_exprs <- BranchAB_exprs[, 1:100]
# BranchB_exprs <- BranchAB_exprs[, 101:200]
# common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
#                                                     setdiff(pData(new_cds)$State,
#                                                             branch_states), ])
# BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells,
#                                                "Pseudotime"])))
# BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells,
#                                         "Pseudotime"]))
# BranchB_num <- BranchA_num
# BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
# BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
# heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1],
#                         BranchB_exprs)
# heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1,
#                                        sd) == 0, ]
# heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix),
#                                  center = TRUE))
# heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
#                                   FALSE, ]
# heatmap_matrix[is.nan(heatmap_matrix)] = 0
# 
# scale_max = 3
# scale_min = -3
# heatmap_matrix[heatmap_matrix > scale_max] = scale_max
# heatmap_matrix[heatmap_matrix < scale_min] = scale_min
# heatmap_matrix_ori <- heatmap_matrix
# heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,
#                                                           1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
# annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
#                              `Cell Type` = c(rep(branch_labels[1], BranchA_num),
#                                              rep("Pre-branch", 2 * BranchP_num),
#                                              rep(branch_labels[2],
#                                                  BranchB_num)))
# branch_colors = c("#979797", "#F05662", "#7990C8")
# csplit <- c(rep('Cell fate 1', times = 100),
#             rep('Cell fate 2', times = 100))
# rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix),
#                                  pattern = 'ENSMUSG[-.0-9]+',
#                                  replacement = '')
# name <- rownames(heatmap_matrix)[1:10]
# library(ComplexHeatmap)
# ha = HeatmapAnnotation(
#   `Cell Type` = annotation_col$Cell.Type,
#   col = list(`Cell Type` = c("Cell fate 1" = branch_colors[2],
#                              "Cell fate 2" = branch_colors[3],
#                              "Pre-branch" = branch_colors[1])),
#   na_col = "grey")
# 
# 
# pdf('moncole_high_val_moncole/branched3/branched3_top100.pdf',height = 6,width = 7)
# Heatmap(heatmap_matrix, cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         column_split = csplit,
#         top_annotation = ha,
#         clustering_distance_rows = "pearson",
#         col = circlize::colorRamp2(c(-2,-1,0,1,2),
#                                    c("blue",'Cyan',
#                                      "green",
#                                      'yellow',
#                                      "red")))+
#   rowAnnotation(link = anno_mark(at = which(rownames(heatmap_matrix) %in% name),
#                                  labels = name,
#                                  labels_gp = gpar(fontsize = 10)))
# dev.off()
# 
# 气泡图top10
# for(i in 1:10){
#   fn = paste0('moncole_high_val_moncole/branched3/', row.names(BEAM_res3t)[i],'.pdf')
#   pps <- data.frame(predictans = cds$raw_predictans,
#                     Pseudotime = cds$Pseudotime,
#                     gene = exprs(cds)[row.names(BEAM_res3t)[i],],
#                     stat = cds$State)
#   pps %>% na.omit() %>% group_by(stat) %>%
#     summarise(mean_predictans = mean(predictans),
#               mean_gene = mean(gene),
#               mean_time = mean(Pseudotime))  %>%
#     ggplot(aes(x = mean_predictans,
#                y = mean_time,
#                color = stat,
#                size = mean_gene)) +
#     geom_point() + theme_bw() +
#     geom_text_repel(mapping = aes(x = mean_predictans,
#                                   y = mean_time,
#                                   label = stat),
#                     size = 4,color = 'black') -> t
#   ggsave(t,filename = fn)
# }




high_val_branch3 <- list()
for(i in 1:10){
  pps <- data.frame(predictans = cds$raw_predictans,
                    Pseudotime = cds$Pseudotime,
                    gene = exprs(cds)[row.names(BEAM_res3t)[i],],
                    stat = cds$State)
  pps %>% na.omit() %>% group_by(stat) %>%
    summarise(mean_predictans = mean(predictans),
              mean_gene = mean(gene),
              mean_time = mean(Pseudotime),
              gene_name = row.names(BEAM_res3t)[i],
              group = 'branch3',
              col_order = i) -> high_val_branch3[[i]]
}
high_val_branch3 <- purrr::reduce(high_val_branch3, rbind)

high_val <- rbind(high_val_total,high_val_branch1,high_val_branch2,high_val_branch3)
library(ggrepel)
high_val$gene_name <- gsub(high_val$gene_name, pattern = 'ENSMUSG[-.0-9]+',replacement = '')
ggplot(high_val, aes(x = mean_predictans,
                     y = mean_time,
                     color = stat,
                     size = mean_gene, group = group)) +
  geom_point() + theme_bw() +
  geom_text_repel(mapping = aes(x = mean_predictans,
                                y = mean_time,
                                label = stat),
                  size = 4,color = 'black') + 
  facet_wrap( ~ group + gene_name, 
              ncol = 5) -> 
  big_val

ggsave(big_val, filename = 'high_val_bouble.pdf',
       height = 20,width = 12,limitsize = FALSE)

# pdf('moncole_high_val_moncole.pdf',5,5)
# cds <- orderCells(cds, root_state = 2)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# Heatmap(small_mat,
#         col = circlize::colorRamp2(c(0,2.5), c("#6495ED","red")),
#         layer_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.f", pindex(table(cds$State,cds$seurat_clusters), i, j)), x, y, gp = gpar(fontsize = 10))
#         }
# )
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans0",
#                              root_state = 2)
# plot_complex_cell_trajectory(cds,
#                              color_by = "predictans_high_low",
#                              root_state = 2)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "seurat_clusters",
#                              root_state = 2)
# plot_complex_cell_trajectory(cds,
#                              color_by = "BHATTACHARYA_EMBRYONIC_STEM_CELL",
#                              root_state = 2)
# plot_complex_cell_trajectory(cds,
#                              color_by = "Pseudotime",
#                              root_state = 2)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "State",
#                              root_state = 2)
# 
# dev.off()
cds$w <- factor(cds$w,levels = c('0w','2w','4w','12w','18w','20w','30w'))
plot_complex_cell_trajectory(cds,
                             color_by = "w",
                             root_state = 2) -> moncolhighval
ggsave(moncolhighval,filename = 'moncolhighval.pdf',height = 7,width = 7)

plot_complex_cell_trajectory(cds,
                             color_by = "State",
                             root_state = 2)
pData(cds) -> cln
rownames(cell_scores_CNV) <- cell_scores_CNV$X
pData(cds) <- cbind(pData(cds),cell_scores_CNV[cln$sampleID,])
plot_complex_cell_trajectory(cds,
                             color_by = "cnv_score",
                             root_state = 2) + 
  scale_color_gradient2(mid = 'blue', high = 'red')



data.frame(barcode = cds$sampleID,
           predictans_high_low = cds$predictans_high_low, 
           State_high_val = cds$State,
           time_high_val = cds$raw_predictans,
           Pseudotime = cds$Pseudotime) -> moncole_color
write.csv(moncole_color, file = 'moncole_color.csv')
# # ----------------------------- cluster差异基因 ---------------------
# deg <- read.csv(paste0(pt,"cluster_express_genes_train.monocle.DEG.csv"))
# cds <- readRDS(paste0(pt,'cluster_express_genes_monocle.Rds'))
# ordergene <- deg$X
# 
# exp_data <- exprs(cds)
# rn <- gsub(rownames(exp_data), pattern = '\\..*', replacement = '')
# rownames(exp_data) <- rn
# cells_rankings <- AUCell_buildRankings(exp_data)
# cells_AUC <- list()
# for(i in unique(stemness$term)){
#   cells_AUC[[i]] <- getAUC(AUCell_calcAUC(stemnessL[[i]], cells_rankings,
#                                           aucMaxRank=nrow(cells_rankings)*0.1))
# }
# ns <- names(cells_AUC)
# cells_AUC <- purrr::reduce(cells_AUC, rbind)
# cells_AUC <- 10^(cells_AUC - 0.1)
# rownames(cells_AUC) <- ns
# cells_AUC <- t(cells_AUC)
# cells_AUC <- as.data.frame(cells_AUC)
# predictans <- predictans[abs(predictans) < 1]
# predictans0 <- 10^(scale(predictans[colnames(cds)]) - 0.1)
# predictans_high_low <- rep('High', times = length(predictans0))
# predictans_high_low[predictans0 < (10^(-0.1)-0.1)] <- 'Low'
# 
# pData(cds) <- cbind(pData(cds), cells_AUC,predictans0,predictans_high_low)
# 
# cds <- orderCells(cds)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# table(cds$State,cds$seurat_clusters)
# pdf('seurate_cluster_moncole.pdf',5,5)
# cds <- orderCells(cds, root_state = 1)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# Heatmap(small_mat,
#         col = circlize::colorRamp2(c(0,2.5), c("#6495ED","red")),
#         layer_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.f", pindex(table(cds$State,cds$seurat_clusters), i, j)), x, y, gp = gpar(fontsize = 10))
#         }
# )
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans0",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans_high_low",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "seurat_clusters",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "BHATTACHARYA_EMBRYONIC_STEM_CELL",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "Pseudotime",
#                              root_state = 1)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "State",
#                              root_state = 1)
# 
# dev.off()
# ----------------------------- stemness genes ---------------------
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle)
library(clusterProfiler)
library(AUCell)
library(ggrepel)

pt <- 'D:/TCGA-LUAD/肺癌数据/HPSC/scRNA/'

predictans <- readRDS('predictans_ans.Rds')
stemness <- read.gmt(paste0(pt,'stemness.gmt'))
stemnessL <- readRDS('stemness_mouse.Rds')
deg <- read.csv(paste0(pt,"stemness_genes_train.monocle.DEG.csv"))
cds <- readRDS(paste0('stemness_genes_monocle.Rds'))
cds <- orderCells(cds, root_state = 6)

# plot(cds@reducedDimS[1,],cds@reducedDimS[2,])
# plot(cds@reducedDimW[,1],cds@reducedDimW[,2])
# plot(cds@reducedDimK[1,],cds@reducedDimK[2,])
# dim(cds@cellPairwiseDistances)


dir.create('./moncole_fig_stemness')
# data.frame(predictans = cds$raw_predictans, State = cds$State) %>%
#   ggplot(aes(x = State, y = predictans, fill = State)) +
#   geom_violin(trim = F) +
#   geom_jitter(width = 0.3,size = 0.1) +
#   theme_bw() -> t
# ggsave(t,filename = './moncole_fig_stemness/violin.pdf',height = 6,width = 7)
# 
data.frame(predictans = cds$raw_predictans, State = cds$State,
           risk_level = cds$predictans_high_low) %>%
  group_by(State,risk_level) %>%
  count -> td
plnb <- list()
for(i in 1:7){
  tmp <- td[td$State == i,]
  binom.test(as.numeric(tmp[1,3]),as.numeric(tmp[1,3]+tmp[2,3])) -> tmp0
  plnb[[i]] <- tmp0$p.value
}
bind <- data.frame(State = factor(1:7,levels = 1:7),p = signif(unlist(plnb),3))
td %>% left_join(bind) -> td0

ggplot(td0,aes(x = State, y = n,
               group = risk_level, color = risk_level,
               shape = risk_level)) +
  geom_point(size = 3,alpha=0.9) +
  theme_bw()+ geom_text(aes(x = State, 
                            y = 1000, 
                            label = p),
                        color = 'black',size = 4) -> t
t
td0
(td0$n[td0$risk_level == 'Low'])/
  (td0 %>% group_by(State) %>% summarise(sum(n)) %>% pull(`sum(n)`)) -> 
  pct
data <- data.frame(State = factor(1:7, levels = 1:7),
                   High = rep(1,times = 7),
                   Low = pct
)
library(scales)
ggplot()+ geom_bar(data, mapping = aes(x = State, 
                                       y = High),
                   stat = 'identity',
                   color = 'black',
                   fill = hue_pal()(2)[1]) + 
  geom_bar(data, mapping = aes(x = State, 
                               y = Low),
           color = 'black',
           stat = 'identity',
           fill = hue_pal()(2)[2]) + 
  geom_hline(
    yintercept = 0.5,
    color = 'red',
    linetype = 2,
    size = 1) +
  theme_bw() -> 
  t0
t/t0 -> t

# ggsave(t,
#   filename = './moncole_fig_stemness/point_path.pdf',
#   height = 6,width = 7)

# ordergene <- deg$X
# exp_data <- exprs(cds)
# rn <- gsub(rownames(exp_data), pattern = '\\..*', replacement = '')
# rownames(exp_data) <- rn
# cells_rankings <- AUCell_buildRankings(exp_data)
# cells_AUC <- list()
# for(i in unique(stemness$term)){
#   cells_AUC[[i]] <- getAUC(AUCell_calcAUC(stemnessL[[i]], cells_rankings,
#                                           aucMaxRank=nrow(cells_rankings)*0.1))
# }
# ns <- names(cells_AUC)
# cells_AUC <- purrr::reduce(cells_AUC, rbind)
# cells_AUC <- 10^(cells_AUC - 0.1)
# rownames(cells_AUC) <- ns
# cells_AUC <- t(cells_AUC)
# cells_AUC <- as.data.frame(cells_AUC)
# predictans <- predictans[abs(predictans) < 1]
# predictans0 <- 10^(scale(predictans[colnames(cds)]) - 0.1)
# predictans_high_low <- rep('High', times = length(predictans0))
# predictans_high_low[predictans0 < (10^(-0.1)-0.1)] <- 'Low'
# raw_predictans <- predictans[cds$sampleID]
# pData(cds) <- cbind(pData(cds), cells_AUC,predictans0,
#                     predictans_high_low,raw_predictans)
# 
# cds <- orderCells(cds)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# saveRDS(cds, file = 'stemness_genes_monocle.Rds')

# diff_test_res <- differentialGeneTest(cds,
#                                       fullModelFormulaStr = "~sm.ns(Pseudotime)")
# write.csv(diff_test_res, file = 'stemness_diff_test_res.csv')
diff_test_res <- read.csv('stemness_diff_test_res.csv')

# ===== ~~~~~~相关基因 =====
# diff_test_res <- diff_test_res[order(diff_test_res$qval,decreasing = F),]
# # 大热图
# memory.limit(62405)
# sig_gene_names <- diff_test_res$X[diff_test_res$qval < 1e-45]
# length(sig_gene_names)
# trend_formula = "~sm.ns(Pseudotime, df=3)"
# cds_subset <- cds[sig_gene_names,]
# pseudocount <- 1
# newdata <- data.frame(Pseudotime = seq(
#   min(pData(cds_subset)$Pseudotime),
#   max(pData(cds_subset)$Pseudotime), length.out = 3891))
# m <- genSmoothCurves(cds_subset, cores = 4,
#                      trend_formula = trend_formula,
#                      relative_expr = T, new_data = newdata)
# m = m[!apply(m, 1, sum) == 0, ]
# m = vstExprs(cds_subset, expr_matrix = m)
# m = m[!apply(m, 1, sd) == 0, ]
# m = Matrix::t(scale(Matrix::t(m), center = TRUE))
# m = m[is.na(row.names(m)) == FALSE, ]
# m[is.nan(m)] = 0
# 
# scale_max = 3
# scale_min = -3
# m[m > scale_max] = scale_max
# m[m < scale_min] = scale_min
# heatmap_matrix <- m
# 
# rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix),
#                                  pattern = 'ENSMUSG[-.0-9]+',
#                                  replacement = '')
# name <- gsub(sig_gene_names[1:10], pattern = 'ENSMUSG[-.0-9]+',
#              replacement = '')
# 
# pdf('./moncole_fig_stemness/stemness_moncole_all.pdf',height = 7,width = 7)
# Heatmap(heatmap_matrix, cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         clustering_distance_rows = "pearson",
#         col = circlize::colorRamp2(c(-2,-1,0,1,2),
#                                    c("blue",
#                                      'Cyan',
#                                      "green",
#                                      'yellow',
#                                      "red")))+
#   rowAnnotation(link = anno_mark(at = which(rownames(
#     heatmap_matrix) %in% name),
#     labels = name,
#     labels_gp = gpar(fontsize = 15)))
# dev.off()

# dir.create('moncole_fig_stemness/diff_test_res')
# for(i in 1:10){
#   fn <- paste0('moncole_fig_stemness/diff_test_res/',sig_gene_names[i],'.pdf')
# 
#     pps <- data.frame(predictans = cds$raw_predictans,
#                       Pseudotime = cds$Pseudotime,
#                       gene = exprs(cds)[sig_gene_names[i],],
#                       stat = cds$State)
#     pps %>% na.omit() %>% group_by(stat) %>%
#       summarise(mean_predictans = mean(predictans),
#                 mean_gene = mean(gene),
#                 mean_time = mean(Pseudotime)) %>%
#       ggplot(aes(x = mean_predictans,
#                      y = mean_time,
#                      color = stat,
#                      size = mean_gene)) +
#       geom_point() + theme_bw() +
#       geom_text_repel(mapping = aes(x = mean_predictans,
#                                     y = mean_time,
#                                     label = stat),
#                       size = 4,color = 'black') -> t
#     ggsave(t, filename =  fn)
# 
# }



stemness_total <- list()
for(i in 1:10){
  pps <- data.frame(predictans = cds$raw_predictans,
                    Pseudotime = cds$Pseudotime,
                    gene = exprs(cds)[sig_gene_names[i],],
                    stat = cds$State)
  pps %>% na.omit() %>% group_by(stat) %>%
    summarise(mean_predictans = mean(predictans),
              mean_gene = mean(gene),
              mean_time = mean(Pseudotime),
              gene_name = sig_gene_names[i],
              group = 'total',
              col_order = i) -> stemness_total[[i]]
}
stemness_total <- purrr::reduce(stemness_total, rbind)




# BEAM_res1 <- BEAM(cds, branch_point = 1, cores = 4)
# BEAM_res2 <- BEAM(cds, branch_point = 2, cores = 4)
# BEAM_res3 <- BEAM(cds, branch_point = 3, cores = 4)
# 
# list(BEAM_res1 = BEAM_res1,
#      BEAM_res2 = BEAM_res2,
#      BEAM_res3 = BEAM_res3) -> BEAM_stemness
# saveRDS(BEAM_stemness, file = 'BEAM_stemness.Rds')

BEAM_stemness <- readRDS('BEAM_stemness.Rds')
# ==================== ~~~~~~~~~~~~~~~ 分叉基因 ~~~~~~~~~ =====================

BEAM_res1t <- BEAM_stemness[[1]]
BEAM_res1t <- BEAM_res1t[order(BEAM_res1t$qval),]
BEAM_res1t <- BEAM_res1t[,c("gene_short_name", "pval", "qval")]
BEAM_res1t <- BEAM_res1t[order(BEAM_res1t$qval),]
rnf <- row.names(BEAM_res1t)[1:100]
dir.create('moncole_fig_stemness/branched1')

# 分叉1小
branch_labels = c("Cell fate 1", "Cell fate 2")
cds_subset <- cds[row.names(subset(BEAM_res1t,
                                   qval < 1e-5)),]
new_cds <- buildBranchCellDataSet(cds_subset,
                                  branch_point = 1,
                                  progenitor_method = "duplicate")
new_cds@dispFitInfo <- cds_subset@dispFitInfo
progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)

col_gap_ind <- 101
newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
trend_formula <- "~sm.ns(Pseudotime, df=3) * Branch"
BranchAB_exprs <- genSmoothCurves(new_cds[,], cores = 4,
                                  trend_formula = trend_formula,
                                  relative_expr = T, new_data = rbind(newdataA,
                                                                      newdataB))
BranchA_exprs <- BranchAB_exprs[, 1:100]
BranchB_exprs <- BranchAB_exprs[, 101:200]
common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
                                                    setdiff(pData(new_cds)$State,
                                                            branch_states), ])
BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells,
                                               "Pseudotime"])))
BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells,
                                        "Pseudotime"]))
BranchB_num <- BranchA_num
BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1],
                        BranchB_exprs)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1,
                                       sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix),
                                 center = TRUE))
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
                                  FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0

scale_max = 3
scale_min = -3
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,
                                                          1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                             `Cell Type` = c(rep(branch_labels[1], BranchA_num),
                                             rep("Pre-branch", 2 * BranchP_num),
                                             rep(branch_labels[2],
                                                 BranchB_num)))
branch_colors = c("#979797", "#F05662", "#7990C8")
csplit <- c(rep('Cell fate 1', times = 100),
            rep('Cell fate 2', times = 100))
rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix),
                                 pattern = 'ENSMUSG[-.0-9]+',
                                 replacement = '')
name <- rownames(heatmap_matrix)[1:10]
library(ComplexHeatmap)
ha = HeatmapAnnotation(
  `Cell Type` = annotation_col$Cell.Type,
  col = list(`Cell Type` = c("Cell fate 1" = branch_colors[2],
                             "Cell fate 2" = branch_colors[3],
                             "Pre-branch" = branch_colors[1])),
  na_col = "grey")

pdf('moncole_fig_stemness/branched1/branched1_top100.pdf',height = 8,width = 7)
Heatmap(heatmap_matrix, cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        column_split = csplit,
        top_annotation = ha,
        clustering_distance_rows = "pearson",
        col = circlize::colorRamp2(c(-2,-1,0,1,2),
                                   c("blue",'Cyan',
                                     "green",
                                     'yellow',
                                     "red")))+
  rowAnnotation(link = anno_mark(at = which(rownames(heatmap_matrix) %in% name),
                                 labels = name,
                                 labels_gp = gpar(fontsize = 10)))
dev.off()

# 气泡图top10
# for(i in 1:10){
#   fn = paste0('moncole_fig_stemness/branched1/', row.names(BEAM_res1t)[i],'.pdf')
#   pps <- data.frame(predictans = cds$raw_predictans,
#                     Pseudotime = cds$Pseudotime,
#                     gene = exprs(cds)[row.names(BEAM_res1t)[i],],
#                     stat = cds$State)
#   pps %>% na.omit() %>% group_by(stat) %>%
#     summarise(mean_predictans = mean(predictans),
#               mean_gene = mean(gene),
#               mean_time = mean(Pseudotime))  %>%
#     ggplot(aes(x = mean_predictans,
#                    y = mean_time,
#                    color = stat,
#                    size = mean_gene)) +
#     geom_point() + theme_bw() +
#     geom_text_repel(mapping = aes(x = mean_predictans,
#                                   y = mean_time,
#                                   label = stat),
#                     size = 4,color = 'black') -> t
#   ggsave(t,filename = fn)
# }


cds$w <- factor(cds$w,levels = c('0w','2w','4w','12w','18w','20w','30w'))
plot_complex_cell_trajectory(cds,
                             color_by = "w",
                             root_state = 6) -> stemnessweek
ggsave(stemnessweek,filename = 'stemnessweek.pdf',height = 7,width = 7)


high_val_branch1 <- list()
for(i in 1:10){
  pps <- data.frame(predictans = cds$raw_predictans,
                    Pseudotime = cds$Pseudotime,
                    gene = exprs(cds)[row.names(BEAM_res1t)[i],],
                    stat = cds$State)
  pps %>% na.omit() %>% group_by(stat) %>%
    summarise(mean_predictans = mean(predictans),
              mean_gene = mean(gene),
              mean_time = mean(Pseudotime),
              gene_name = row.names(BEAM_res1t)[i],
              group = 'branch1',
              col_order = i) -> high_val_branch1[[i]]
}
high_val_branch1 <- purrr::reduce(high_val_branch1, rbind)




BEAM_res2t <- BEAM_stemness[[2]]
BEAM_res2t <- BEAM_res2t[order(BEAM_res2t$qval),]
BEAM_res2t <- BEAM_res2t[,c("gene_short_name", "pval", "qval")]
BEAM_res2t <- BEAM_res2t[order(BEAM_res2t$qval),]
rnf <- row.names(BEAM_res2t)[1:100]
dir.create('moncole_fig_stemness/branched2')

# 分叉2小
branch_labels = c("Cell fate 1", "Cell fate 2")
cds_subset <- cds[row.names(subset(BEAM_res2t,
                                   qval < 1e-5)),]
new_cds <- buildBranchCellDataSet(cds_subset,
                                  branch_point = 2,
                                  progenitor_method = "duplicate")
new_cds@dispFitInfo <- cds_subset@dispFitInfo
progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)

col_gap_ind <- 101
newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
trend_formula <- "~sm.ns(Pseudotime, df=3) * Branch"
BranchAB_exprs <- genSmoothCurves(new_cds[,], cores = 4,
                                  trend_formula = trend_formula,
                                  relative_expr = T, new_data = rbind(newdataA,
                                                                      newdataB))
BranchA_exprs <- BranchAB_exprs[, 1:100]
BranchB_exprs <- BranchAB_exprs[, 101:200]
common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
                                                    setdiff(pData(new_cds)$State,
                                                            branch_states), ])
BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells,
                                               "Pseudotime"])))
BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells,
                                        "Pseudotime"]))
BranchB_num <- BranchA_num
BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1],
                        BranchB_exprs)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1,
                                       sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix),
                                 center = TRUE))
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
                                  FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0

scale_max = 3
scale_min = -3
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,
                                                          1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                             `Cell Type` = c(rep(branch_labels[1], BranchA_num),
                                             rep("Pre-branch", 2 * BranchP_num),
                                             rep(branch_labels[2],
                                                 BranchB_num)))
branch_colors = c("#979797", "#F05662", "#7990C8")
csplit <- c(rep('Cell fate 1', times = 100),
            rep('Cell fate 2', times = 100))
rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix),
                                 pattern = 'ENSMUSG[-.0-9]+',
                                 replacement = '')
name <- rownames(heatmap_matrix)[1:10]
library(ComplexHeatmap)
ha = HeatmapAnnotation(
  `Cell Type` = annotation_col$Cell.Type,
  col = list(`Cell Type` = c("Cell fate 1" = branch_colors[2],
                             "Cell fate 2" = branch_colors[3],
                             "Pre-branch" = branch_colors[1])),
  na_col = "grey")

pdf('moncole_fig_stemness/branched2/branched2_top100.pdf',height = 8,width = 7)
Heatmap(heatmap_matrix, cluster_columns = F,
        show_row_names = F,
        show_column_names = F,
        column_split = csplit,
        top_annotation = ha,
        clustering_distance_rows = "pearson",
        col = circlize::colorRamp2(c(-2,-1,0,1,2),
                                   c("blue",'Cyan',
                                     "green",
                                     'yellow',
                                     "red")))+
  rowAnnotation(link = anno_mark(at = which(rownames(heatmap_matrix) %in% name),
                                 labels = name,
                                 labels_gp = gpar(fontsize = 10)))
dev.off()
# 
# 气泡图top10
# for(i in 1:10){
#   fn = paste0('moncole_fig_stemness/branched2/', row.names(BEAM_res2t)[i],'.pdf')
#   pps <- data.frame(predictans = cds$raw_predictans,
#                     Pseudotime = cds$Pseudotime,
#                     gene = exprs(cds)[row.names(BEAM_res2t)[i],],
#                     stat = cds$State)
#   pps %>% na.omit() %>% group_by(stat) %>%
#     summarise(mean_predictans = mean(predictans),
#               mean_gene = mean(gene),
#               mean_time = mean(Pseudotime))  %>%
#     ggplot(aes(x = mean_predictans,
#                y = mean_time,
#                color = stat,
#                size = mean_gene)) +
#     geom_point() + theme_bw() +
#     geom_text_repel(mapping = aes(x = mean_predictans,
#                                   y = mean_time,
#                                   label = stat),
#                     size = 4,color = 'black') -> t
#   ggsave(t,filename = fn)
# }

high_val_branch2 <- list()
for(i in 1:10){
  pps <- data.frame(predictans = cds$raw_predictans,
                    Pseudotime = cds$Pseudotime,
                    gene = exprs(cds)[row.names(BEAM_res2t)[i],],
                    stat = cds$State)
  pps %>% na.omit() %>% group_by(stat) %>%
    summarise(mean_predictans = mean(predictans),
              mean_gene = mean(gene),
              mean_time = mean(Pseudotime),
              gene_name = row.names(BEAM_res2t)[i],
              group = 'branch2',
              col_order = i) -> high_val_branch2[[i]]
}
high_val_branch2 <- purrr::reduce(high_val_branch2, rbind)



BEAM_res3t <- BEAM_stemness[[3]]
BEAM_res3t <- BEAM_res3t[order(BEAM_res3t$qval),]
BEAM_res3t <- BEAM_res3t[,c("gene_short_name", "pval", "qval")]
BEAM_res3t <- BEAM_res3t[order(BEAM_res3t$qval),]
rnf <- row.names(BEAM_res3t)[1:100]
dir.create('moncole_fig_stemness/branched3')

# 分叉3小
branch_labels = c("Cell fate 1", "Cell fate 2")
cds_subset <- cds[row.names(subset(BEAM_res3t,
                                   qval < 1e-5)),]
new_cds <- buildBranchCellDataSet(cds_subset,
                                  branch_point = 3,
                                  progenitor_method = "duplicate")
new_cds@dispFitInfo <- cds_subset@dispFitInfo
progenitor_state <- subset(pData(cds_subset), Pseudotime == 0)[, "State"]
branch_states <- setdiff(pData(cds_subset)$State, progenitor_state)

col_gap_ind <- 101
newdataA <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[1]))
newdataB <- data.frame(Pseudotime = seq(0, 100, length.out = 100),
                       Branch = as.factor(unique(as.character(pData(new_cds)$Branch))[2]))
trend_formula <- "~sm.ns(Pseudotime, df=3) * Branch"
BranchAB_exprs <- genSmoothCurves(new_cds[,], cores = 4,
                                  trend_formula = trend_formula,
                                  relative_expr = T, new_data = rbind(newdataA,
                                                                      newdataB))
BranchA_exprs <- BranchAB_exprs[, 1:100]
BranchB_exprs <- BranchAB_exprs[, 101:200]
common_ancestor_cells <- row.names(pData(new_cds)[pData(new_cds)$State ==
                                                    setdiff(pData(new_cds)$State,
                                                            branch_states), ])
BranchP_num <- (100 - floor(max(pData(new_cds)[common_ancestor_cells,
                                               "Pseudotime"])))
BranchA_num <- floor(max(pData(new_cds)[common_ancestor_cells,
                                        "Pseudotime"]))
BranchB_num <- BranchA_num
BranchA_exprs <- vstExprs(new_cds, expr_matrix = BranchA_exprs)
BranchB_exprs <- vstExprs(new_cds, expr_matrix = BranchB_exprs)
heatmap_matrix <- cbind(BranchA_exprs[, (col_gap_ind - 1):1],
                        BranchB_exprs)
heatmap_matrix = heatmap_matrix[!apply(heatmap_matrix, 1,
                                       sd) == 0, ]
heatmap_matrix = Matrix::t(scale(Matrix::t(heatmap_matrix),
                                 center = TRUE))
heatmap_matrix = heatmap_matrix[is.na(row.names(heatmap_matrix)) ==
                                  FALSE, ]
heatmap_matrix[is.nan(heatmap_matrix)] = 0

scale_max = 3
scale_min = -3
heatmap_matrix[heatmap_matrix > scale_max] = scale_max
heatmap_matrix[heatmap_matrix < scale_min] = scale_min
heatmap_matrix_ori <- heatmap_matrix
heatmap_matrix <- heatmap_matrix[is.finite(heatmap_matrix[,
                                                          1]) & is.finite(heatmap_matrix[, col_gap_ind]), ]
annotation_col <- data.frame(row.names = c(1:ncol(heatmap_matrix)),
                             `Cell Type` = c(rep(branch_labels[1], BranchA_num),
                                             rep("Pre-branch", 2 * BranchP_num),
                                             rep(branch_labels[2],
                                                 BranchB_num)))
branch_colors = c("#979797", "#F05662", "#7990C8")
csplit <- c(rep('Cell fate 1', times = 100),
            rep('Cell fate 2', times = 100))
rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix),
                                 pattern = 'ENSMUSG[-.0-9]+',
                                 replacement = '')
name <- rownames(heatmap_matrix)[1:10]
library(ComplexHeatmap)
ha = HeatmapAnnotation(
  `Cell Type` = annotation_col$Cell.Type,
  col = list(`Cell Type` = c("Cell fate 1" = branch_colors[2],
                             "Cell fate 2" = branch_colors[3],
                             "Pre-branch" = branch_colors[1])),
  na_col = "grey")


# pdf('moncole_fig_stemness/branched3/branched3_top100.pdf',height = 8,width = 7)
# Heatmap(heatmap_matrix, cluster_columns = F,
#         show_row_names = F,
#         show_column_names = F,
#         column_split = csplit,
#         top_annotation = ha,
#         clustering_distance_rows = "pearson",
#         col = circlize::colorRamp2(c(-2,-1,0,1,2),
#                                    c("blue",'Cyan',
#                                      "green",
#                                      'yellow',
#                                      "red")))+
#   rowAnnotation(link = anno_mark(at = which(rownames(heatmap_matrix) %in% name),
#                                  labels = name,
#                                  labels_gp = gpar(fontsize = 10)))
# dev.off()
# 
# 气泡图top10
# for(i in 1:10){
#   fn = paste0('moncole_fig_stemness/branched3/', row.names(BEAM_res3t)[i],'.pdf')
#   pps <- data.frame(predictans = cds$raw_predictans,
#                     Pseudotime = cds$Pseudotime,
#                     gene = exprs(cds)[row.names(BEAM_res3t)[i],],
#                     stat = cds$State)
#   pps %>% na.omit() %>% group_by(stat) %>%
#     summarise(mean_predictans = mean(predictans),
#               mean_gene = mean(gene),
#               mean_time = mean(Pseudotime))  %>%
#     ggplot(aes(x = mean_predictans,
#                y = mean_time,
#                color = stat,
#                size = mean_gene)) +
#     geom_point() + theme_bw() +
#     geom_text_repel(mapping = aes(x = mean_predictans,
#                                   y = mean_time,
#                                   label = stat),
#                     size = 4,color = 'black') -> t
#   ggsave(t,filename = fn)
# }



high_val_branch3 <- list()
for(i in 1:10){
  pps <- data.frame(predictans = cds$raw_predictans,
                    Pseudotime = cds$Pseudotime,
                    gene = exprs(cds)[row.names(BEAM_res3t)[i],],
                    stat = cds$State)
  pps %>% na.omit() %>% group_by(stat) %>%
    summarise(mean_predictans = mean(predictans),
              mean_gene = mean(gene),
              mean_time = mean(Pseudotime),
              gene_name = row.names(BEAM_res3t)[i],
              group = 'branch3',
              col_order = i) -> high_val_branch3[[i]]
}
high_val_branch3 <- purrr::reduce(high_val_branch3, rbind)



stemness_gather <- rbind(stemness_total,high_val_branch1,high_val_branch2,high_val_branch3)
library(ggrepel)
stemness_gather$gene_name <- gsub(stemness_gather$gene_name, pattern = 'ENSMUSG[-.0-9]+',replacement = '')
ggplot(stemness_gather, aes(x = mean_predictans,
                            y = mean_time,
                            color = stat,
                            size = mean_gene, group = group)) +
  geom_point() + theme_bw() +
  geom_text_repel(mapping = aes(x = mean_predictans,
                                y = mean_time,
                                label = stat),
                  size = 4,color = 'black') + 
  facet_wrap( ~ group + gene_name, 
              ncol = 5) -> 
  big_val

ggsave(big_val, filename = 'stemness_gather_bouble.pdf',
       height = 20,width = 12,limitsize = FALSE)
# pdf('stemness_moncole.pdf',5,5)
# cds <- orderCells(cds, root_state = 6)
# small_mat <- log10(1 + table(cds$State,cds$seurat_clusters))
# Heatmap(small_mat,
#         col = circlize::colorRamp2(c(0,2.5), c("#6495ED","red")),
#         layer_fun = function(j, i, x, y, width, height, fill) {
#           grid.text(sprintf("%.f", pindex(table(cds$State,cds$seurat_clusters), i, j)), x, y, gp = gpar(fontsize = 10))
#         }
# )
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans0",
#                              root_state = 6)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "predictans_high_low",
#                              root_state = 6)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "seurat_clusters",
#                              root_state = 6)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "BHATTACHARYA_EMBRYONIC_STEM_CELL",
#                              root_state = 6)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "Pseudotime",
#                              root_state = 6)
# plot_complex_cell_trajectory(cds, 
#                              color_by = "State",
#                              root_state = 6)
# 
# dev.off()
data.frame(barcode = cds$sampleID,
           predictans_high_low = cds$predictans_high_low, 
           State_high_val = cds$State,
           time_high_val = cds$raw_predictans,
           Pseudotime = cds$Pseudotime) -> moncole_stemness
write.csv(moncole_stemness, file = 'moncole_stemness_color.csv')
# ============================== moncole2 拟风险 ===================================
library(monocle)
pt <- 'D:/TCGA-LUAD/肺癌数据/HPSC/scRNA/'

predictans <- readRDS('predictans_ans.Rds')
HPSC <- readRDS(paste0(pt,'HPSC_obj.Rds'))
expr_matrix <- as(as.matrix(HPSC@assays$SCT@data), 'sparseMatrix')
p_data <- HPSC@meta.data 
p_data$celltype <- HPSC@meta.data$SCT_snn_res.0.5  ##整合每个细胞的细胞鉴定信息到p_data里面。如果已经添加则不必重复添加
p_data$predictans <- predictans[HPSC$sampleID]
f_data <- data.frame(gene_short_name = row.names(HPSC),
                     row.names = row.names(HPSC))

pd <- new('AnnotatedDataFrame', data = p_data) 
fd <- new('AnnotatedDataFrame', data = f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr = 0.1) #这一操作会在fData(cds)中添加一列num_cells_expressed

pData(cds)$Risk_Level <- rep('High',times = dim(pData(cds))[1])
pData(cds)$Risk_Level[pData(cds)$predictans < -0.1] <- 'Low'

cor_risk <- 
  apply(expr_matrix,1,function(x) {
    cor.test(x,cds$predictans) -> t
    c(p = t$p.value, R = t$estimate)
  })
cor_risk <- t(cor_risk)
cor_risk <- cor_risk[order(cor_risk[,1]),]
sig_gene_names <- rownames(cor_risk[cor_risk[,1] < 1e-25,])
write.csv(cor_risk,file = 'cor_risk_相关热图的基因.csv')
length(sig_gene_names)
library(ComplexHeatmap)

trend_formula = "~sm.ns(predictans, df=3)"
cds_subset <- cds[sig_gene_names,]
cds_subset <- cds_subset[,order(cds_subset$predictans)]
pData(cds_subset) <- pData(cds_subset)[
  order(cds_subset$predictans),]

pseudocount <- 1
newdata <- data.frame(predictans = seq(
  min(pData(cds_subset)$predictans), 
  max(pData(cds_subset)$predictans), length.out = 3891))
m <- genSmoothCurves(cds_subset, cores = 4, 
                     trend_formula = trend_formula,
                     relative_expr = T, new_data = newdata)

m = m[!apply(m, 1, sum) == 0,]
m = vstExprs(cds_subset, expr_matrix = m)
m = m[!apply(m, 1, sd) == 0, ]
m = Matrix::t(scale(Matrix::t(m), center = TRUE))
m = m[is.na(row.names(m)) == FALSE, ]
m[is.nan(m)] = 0
scale_max = 3
scale_min = -3
m[m > scale_max] = scale_max
m[m < scale_min] = scale_min
heatmap_matrix <- m
rownames(heatmap_matrix) <- gsub(rownames(heatmap_matrix), 
                                 pattern = 'ENSMUSG[-.0-9]+', 
                                 replacement = '')
name <- gsub(sig_gene_names[1:10], pattern = 'ENSMUSG[-.0-9]+', 
             replacement = '')
View(sig_gene_names)
ha <- HeatmapAnnotation(
  Risk_Level = cds_subset$Risk_Level,
  col = list(
    Risk_Level = c(
      'Low' = 'blue',
      'High' = 'red')
  )
)
pdf('./diff_risk_all.pdf',height = 7,width = 7)
Heatmap(heatmap_matrix, 
        cluster_columns = F, 
        show_row_names = F,
        show_column_names = F,
        top_annotation = ha,
        clustering_distance_rows = "pearson",
        col = circlize::colorRamp2(c(-2,-1,0,1,2), 
                                   c("blue",
                                     'Cyan',
                                     "green",
                                     'yellow',
                                     "red")))+ 
  rowAnnotation(link = anno_mark(at = which(rownames(
    heatmap_matrix) %in% name), 
    labels = name, 
    labels_gp = gpar(fontsize = 15)))
dev.off()


# ================================= HPSC_marker ====================================
markers <- read.csv('HPSC_marker.csv')
markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) -> markers0
markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) -> markers2

pdf('feature.pdf', width = 10, height = 10)
DoHeatmap(pbmc, features = markers0$gene, size = 5, angle = -50, hjust=0.8) + NoLegend()
dev.off()
markers0$gene

marker_mx <- as.matrix(HPSC@assays$SCT@data[markers0$gene,])
marker_mx <- as.data.frame(t(marker_mx))
marker_mx$cluster <- HPSC$SCT_snn_res.0.5
marker_mx <- melt(marker_mx)
marker_mx %>% group_by(cluster, variable) %>% 
  summarise(mean_score = mean(value)) -> fm
library(RColorBrewer)
pdf('boble.pdf',6.5,10)
ggplot(fm,aes(x = cluster, y = variable, size = mean_score, color = mean_score)) + 
  geom_point() + theme_bw() + scale_color_gradientn(colours = rev(c("#FFD92F","#FEE391",brewer.pal(11, "Spectral")[7:11])))
dev.off()
# ========================== 上游 ==============================
library(tidyverse)
rbind(read.csv('sample.csv'),read.csv('sample (1).csv')) -> tmp
tmp0 <- tmp[,c(1,2)]
tmp0$Title <- str_extract(tmp0$Title, pattern = '(?<=: ).+')
write.csv(tmp0,file = 'sampletmp0.csv')
tmp$sample <- str_extract(tmp$Title, pattern = '(?<=: ).+')
read.table('SraRunTable.txt',sep = ',',header = T) -> tb

HPSC <- readRDS('../../肺癌数据/HPSC/scRNA/HPSC_obj.Rds')

tmp$Accession[
  tmp$sample %in% colnames(HPSC)[colnames(HPSC) %in% tmp$sample]
] -> GSM

tb$Run[tb$GEO_Accession..exp. %in% GSM] -> SRR
tb$GEO_Accession..exp.[tb$GEO_Accession..exp. %in% GSM] -> SGMid
tb$timepoint[tb$GEO_Accession..exp. %in% GSM] -> timepoint
data.frame(SRR = SRR,
           GSM = SGMid,
           timepoint = timepoint) -> download
download$gt <- str_extract(download$timepoint, pattern = '.+?(?=_)')
download$time <- str_extract(download$timepoint, pattern = '(?<=_).+?(?=w_)')

gn <- names(table(download$GSM))
singled <- download[download$GSM %in% gn[table(download$GSM) == 1],]
doubled <- download[download$GSM %in% gn[table(download$GSM) == 2],]


write.table(singled, file = 'singledtable.txt',quote = F,row.names = F,col.names = F)
write.table(doubled, file = 'doubledtable.txt',quote = F,row.names = F,col.names = F)























