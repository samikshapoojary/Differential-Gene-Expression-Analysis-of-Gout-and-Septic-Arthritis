#installing packages
install.packages("ggplot2")
install.packages("ggfortify")

library(ggplot2)
library(ggfortify)

#opening all the files
anno = read.table("Annotations.csv", header=TRUE,row.names=1,sep="\t")
de_g_vs_hc = read.table("DE_GOUT_vs_HC.csv", header=TRUE,row.names=1,sep="\t")
de_sa_vs_hc = read.table("DE_SA_vs_HC.csv", header=TRUE,row.names=1,sep="\t")
exp = read.table("Expression_Table.csv", header=TRUE,row.names=1,sep="\t")
samp_info = read.table("Sample_Information.csv", header=TRUE,row.names=1,sep="\t")

#count of genes whose expression levels were measured
nrow(exp) #29814-1

#separating by sample group
samp_hc =subset (samp_info, SAMPLE_GROUP == "HC")
samp_g =subset (samp_info, SAMPLE_GROUP == "GOUT")
samp_sa = subset (samp_info, SAMPLE_GROUP == "SEPSIS")

exp_hc =subset (exp[,1:9])
exp_g =subset (exp[,10:18])
exp_sa = subset (exp[,19:27])

#QUESTION 1

#ratio of males to female
summary(samp_info)

samp_males = subset (samp_info, SEX == "M") #13 samples
hc_males = subset (samp_info, SEX == "M" & SAMPLE_GROUP == "HC" )#5 samples
g_males = subset (samp_info, SEX == "M" & SAMPLE_GROUP == "GOUT" ) #4 samples
sa_males = subset (samp_info, SEX == "M" & SAMPLE_GROUP == "SEPSIS" )#4 samples

samp_females = subset (samp_info, SEX == "F") #14 samples
hc_females = subset (samp_info, SEX == "F" & SAMPLE_GROUP == "HC" )#4 samples
g_females = subset (samp_info, SEX == "F" & SAMPLE_GROUP == "GOUT" )#5 samples
sa_females = subset (samp_info, SEX == "F" & SAMPLE_GROUP == "SEPSIS" )#5 samples

#comparison of neutrophil count
summary(samp_males) #mean neutrophils = 8.538
summary(samp_females) #mean neutrophils = 8.336

p_neutro_g_vs_hc= t.test(samp_hc$NEUTROPHILS,samp_g$NEUTROPHILS) 
p_neutro_sa_vs_hc = t.test(samp_hc$NEUTROPHILS,samp_sa$NEUTROPHILS)
p_neutro_g_vs_hc= p_neutro_g_vs_hc$p.value
p_neutro_sa_vs_hc= p_neutro_sa_vs_hc$p.value

p_neutro_g_vs_hc # 0.1108907
p_neutro_sa_vs_hc #1.219768e-05 = 0.0000122


#QUESTION 2

#finding significant genes in gout and SA respectively 
gout_sig = subset(de_g_vs_hc, p.adj < 0.05) 
nrow(gout_sig) #69-1
sa_sig = subset(de_sa_vs_hc, p.adj <0.05)
nrow(sa_sig) #13046-1

#getting expression of all significant genes of gout and SA
gout_exp = merge(gout_sig,exp_g, by.x=0, by.y=0)
gene_id_g = data.frame(gout_exp[,1])
gene_exp_g = data.frame(gout_exp[,5:13])
all_g_sig_exp = cbind(gene_id_g, gene_exp_g)
row.names(all_g_sig_exp) = all_g_sig_exp[,1]
all_g_sig_exp = all_g_sig_exp[,-1]

sa_exp =merge(sa_sig,exp_sa, by.x=0, by.y=0)
gene_id_sa = data.frame(sa_exp[,1])
gene_exp_sa = data.frame(sa_exp[,5:13])
all_sa_sig_exp = cbind(gene_id_sa, gene_exp_sa)
row.names(all_sa_sig_exp) = all_sa_sig_exp[,1]
all_sa_sig_exp = all_sa_sig_exp[,-1]

#pca on the two data sets to check for\clustering
pca_g= prcomp(t(all_g_sig_exp))
pca_coordinates_g = data.frame(pca_g$x)
pca_coordinates_g$sample = "Gout"

pca_sa= prcomp(t(all_sa_sig_exp))
pca_coordinates_sa = data.frame(pca_sa$x)
pca_coordinates_sa$sample = "SA"

both_pca_coord = rbind(pca_coordinates_g, pca_coordinates_sa)
pca_plot = ggplot(both_pca_coord, aes(x = PC1, y = PC2, color = sample)) + geom_point() 
pca_plot

#finding common significant genes between the two
g_sa_sig = merge(gout_sig,sa_sig, by.x=0, by.y=0)
names(g_sa_sig) = c("GENE_ID","Log2fold_G", "p_G","p.adj_G","Log2fold_SA","p_SA","p.adj_SA")  
row.names(g_sa_sig) = g_sa_sig[,1]
nrow (g_sa_sig) #47-1

#the most significantly differential genes between HC and Gout/SA
most_sig_g = rownames(gout_sig)[which.min(gout_sig$p)]
most_sig_g #ENSG00000179023 - as it has lowest p-value even though same p.adj value as others
most_sig_sa = rownames(sa_sig)[which.min(sa_sig$p.adj)]
most_sig_sa #ENSG00000198074 - lowest p and p.adj value

#plotting the most significantly differential genes to check their distribution and similarity
g_sig_exp = data.frame(t(exp_g['ENSG00000179023',]))
sa_sig_exp = data.frame(t(exp_sa['ENSG00000198074',]))

g1 = ggplot(g_sig_exp, aes(x=log10(ENSG00000179023))) + geom_histogram(colour="black", fill="pink", linewidth=0.5, alpha=1.0, bins =10) + labs(x="Gout Samples", y="Expression levels of ENSG00000179023", title="Expression levels of ENSG00000179023 in Gout Patients")
sa1 = ggplot(sa_sig_exp, aes(x=log10(ENSG00000198074))) + geom_histogram(colour="black", fill="lavender", linewidth=0.5, alpha=1.0, bins = 10) + labs(x="SA Samples", y="Expression levels of ENSG00000198074", title="Expression levels of ENSG00000198074 in SA Patients")
g1
sa1 #both are roughly normally distributed

#making violin boxplots for the two genes 
compare = data.frame(Expression = c(g_sig_exp[,1],sa_sig_exp[,1]), Samples = c("Gout", "SA"))
ggp1 = ggplot(compare, aes(x = Samples, y = Expression, fill = Samples)) +
  geom_violin() + geom_boxplot(width=0.1, alpha=0.2) +
  stat_summary(fun=mean, colour="black") +
  labs(title = " Gene Expression", x = "Samples", y = "Expression Level")
ggp1
#violin boxplots for the two are different - genes not similar


#QUESTION 3

#to check if gene expression is affected by sex and neutrophil count (as age and monocyte not given)
gout_clinical = merge(g_sig_exp,samp_info, by.x=0, by.y=0)
sa_clinical = merge(sa_sig_exp,samp_info, by.x=0, by.y=0)

#H0 : sex does not affect gene exp of ENSG00000179023
cor(as.numeric(factor(gout_clinical$SEX)),gout_clinical$ENSG00000179023)
# -0.4969473 - imperfect negative correlation
model1 = lm(gout_clinical$ENSG00000179023 ~ as.numeric(factor(gout_clinical$SEX)))
anova(model1)
#Pr(>F) = 0.1735 - as p value more than 0.05, accept Ho
summary(model1)
#r_sq value is 0.247, ie 24.7% variance explained by sex

#H0: neutrophil count does not affect gene exp of ENSG00000179023
cor(gout_clinical$NEUTROPHILS,gout_clinical$ENSG00000179023)
#0.1735595- imperfect positive correlation
model2 = lm(gout_clinical$ENSG00000179023 ~ gout_clinical$NEUTROPHILS)
anova(model2)
#Pr(>F) = 0.6552 - as p value more than 0.05, accept Ho
summary(model2)
#r_sq value is 0.03012, ie 3.012% variance explained by neutrophils

#H0: sex does not affect gene exp of ENSG00000198074
cor(as.numeric(factor(sa_clinical$SEX)),sa_clinical$ENSG00000198074)
#-0.04762512- imperfect negative correlation
model3 = lm(sa_clinical$ENSG00000198074 ~ as.numeric(factor(sa_clinical$SEX)))
anova(model3)
#Pr(>F) = 0.9032 - as p value more than 0.05, accept Ho
summary(model3)
#r_sq value is 0.002268, ie 0.23% variance explained by sex

#H0: neutrophil count does not affect gene exp of ENSG00000198074
cor(sa_clinical$NEUTROPHILS,sa_clinical$ENSG00000198074)
#0.2572772 - imperfect positive correlation
model4 = lm(sa_clinical$ENSG00000198074 ~ sa_clinical$NEUTROPHILS)
anova(model4)
#Pr(>F) = 0.5039 - as p value more than 0.05, accept Ho
#as p value less than 0.05, reject Ho
summary(model4)
#r_sq value is 0.06619, ie 6.619% variance explained by neutrophil


#QUESTION 4

#finding significantly different genes between gout and SA 
de_g_vs_sa = data.frame(Gene_ID = rownames(exp), log2Fold = NA, p_g_vs_sa = NA, p.adj = NA)

for (row in 1:nrow(de_g_vs_sa)) 
{
  gout = as.numeric(exp_g[row, ])
  sa = as.numeric(exp_sa[row, ])
  mean_g = mean(gout)
  mean_sa = mean(sa)
  log2fold = log2(mean_g) - log2(mean_sa) 
  p = t.test(gout,sa) 
  p = p$p.value
  de_g_vs_sa[row,"log2Fold"] = log2fold
  de_g_vs_sa[row,"p_g_vs_sa"] = p
  de_g_vs_sa[row,"p.adj"] = p.adjust(p)
  
}

sig_genes = subset(de_g_vs_sa, p.adj < 0.05) 
nrow(sig_genes) #13233-1
most_sig = sig_genes[which.min(sig_genes$p.adj), 'Gene_ID']
most_sig #ENSG00000146242

# Create a volcano plot
volcano_sig = ggplot(de_g_vs_sa, aes(x = log2Fold, y = -log10(p.adj))) + geom_point(aes(color = p.adj < 0.05), alpha = 0.5) +
  labs(title = "Expression of Significantly different genes between Gout and SA", x = "Log2Fold Change", y = "-Log10 of Adjusted P-value") 
volcano_sig

#plotting expression of the most significantly differential gene between Gout and SA
gsa_sig_gexp = data.frame(t(exp_g['ENSG00000146242',]))
gsa_sig_saexp = data.frame(t(exp_sa['ENSG00000146242',]))

g2 = ggplot(gsa_sig_gexp, aes(x=log10(ENSG00000146242))) + geom_histogram(colour="black", fill="pink", linewidth=0.5, alpha=1.0, bins =10) + labs(x="Gout Samples", y="Expression levels of ENSG00000146242", title="Expression levels of ENSG00000146242 in Gout Patients")
sa2 = ggplot(gsa_sig_saexp, aes(x=log10(ENSG00000146242))) + geom_histogram(colour="black", fill="lavender", linewidth=0.5, alpha=1.0, bins = 10) + labs(x="SA Samples", y="Expression levels of ENSG00000146242", title="Expression levels of ENSG00000146242 in SA Patients")
g2 #skewed to right
sa2 #bimodal
#plots are different

#Wilcoxon test to determine significant difference in gene expression between the two groups
w = wilcox.test(gsa_sig_gexp[,1], gsa_sig_saexp[,1])
w
#p value < 0.05 so null hypothesis rejected

#making violin boxplots for same gene 
compare2 = data.frame(Expression = c(gsa_sig_gexp[,1],gsa_sig_saexp[,1]), Samples = c("Gout", "SA"))
ggp2 = ggplot(compare2, aes(x = Samples, y = Expression, fill = Samples)) +
  geom_violin() + geom_boxplot(width=0.1, alpha=0.2) +
  stat_summary(fun=mean, colour="black") +
  labs(title = " Gene Expression of ENSG00000146242", x = "Samples", y = "Expression Level")
ggp2 #mean and median of gout are significantly lower than SA


