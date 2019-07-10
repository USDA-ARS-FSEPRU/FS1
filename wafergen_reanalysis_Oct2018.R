#nr re-analysis of wafergen data so it is all combined and to include a limit of detection
#note that day 7 and day 14 spreadsheets were formatted differently
#also, there are shared pigs between day 7 and 14 therefore have to keep them separate
#my initial plasmidome metadata sheet was missing two control pigs(due to bad sequence in lane 1) so had to redo the data with all pigs


#install_github('Jtrachsel/funfuns')
library(funfuns)
library(tidyverse)
library(ggplot2)

# setwd("/Users/nicolericker/Documents/Pig experiments/FY2016/Wafergen samples")
day7_fecal<-read.csv("wafergen_CTs.csv", header=TRUE, stringsAsFactors=FALSE)
#define a rep column and populate it with rep1 and rep2
day7_fecal$rep<- sub('.*(Rep[1-2])', '\\1', day7_fecal$gene)
head(day7_fecal$rep)

#now remove the rep from the original gene name
day7_fecal$gene<-sub('(.*)_Rep[1-2]', '\\1', day7_fecal$gene)
head(day7_fecal$gene)

#then order the data by gene to match reps
day7_fecal <- day7_fecal[order(day7_fecal$gene),]
head(day7_fecal$gene)
dim(day7_fecal)  #96 56
#convert all cycles above 28 to NA and remove rows with NA all the way across
day7_fecal[, 2:55][day7_fecal[,2:55] > 28] <-NA
day7_fecal <- day7_fecal[which(!(rowSums(is.na(day7_fecal)) == max(rowSums(is.na(day7_fecal))))),]
dim(day7_fecal)   #72  56

#next put this in long format to do the QC on missing replicates and high variability
#the starts_with('X') is b/c all the sample numbers had an X appended to the beginning
#then Jules linked the sample with the gene and summarized them
day7_fecal.gather <- gather(day7_fecal, key=sample_ID, value=CT, starts_with('X')) %>%
  mutate(sample_gene=paste(sample_ID, gene, sep='_')) %>%
  group_by(sample_gene) %>% summarise(meanCT=mean(CT),
                                      diff_btw_reps=sort(CT)[2]-sort(CT)[1],
                                      gene=gene[1],
                                      sample_ID=sample_ID[1],
                                      rep1=CT[1],
                                      rep2=CT[2],
                                      numbad=sum(is.na(CT)),
                                      onegoodCT=ifelse(sum(is.na(CT))==1,CT[which(!is.na(CT))], NA))
write.table(day7_fecal.gather, file="day7_fecal_gather.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
day14<-read.csv("Wafergen-March2018.csv", header=TRUE, stringsAsFactors=FALSE)
#note this has day14 fecals and day7 plasmidome data

#parse out just the assay, sample and Ct values then mutate to include the reps
day14 <- day14[,c(3,4,6)]
library(stringr)  
day14$Rep <- str_sub(day14$Assay, -4)
day14$Rep=parse_number(day14$Rep)
#now we need to remove the Rep information from the gene name
day14 <- separate(day14, Assay, into=c("Gene", "Primer"), extra = "drop", sep="_", fill="warn", remove=TRUE)
day14 <- unite(day14, Genes, Gene, Primer, sep="-", remove=TRUE)
head(day14)
#it would appear that the ermA genes had "-" instead of "_" for their Rep
unique(day14$Genes)
day14$Genes <- sub("-Rep[1-2]", "", day14$Genes)

#collect only the day7 plasmidome samples
plasmidome<- day14 %>% filter(str_detect(Sample, "o$"))
unique(plasmidome$Sample)

#and separate out the day14 fecal samples
day14_fecal<- day14 %>% filter(!str_detect(Sample, "o$"))
unique(day14_fecal$Sample)
unique(day7_fecal.gather$sample_ID)
#Analyze the fecals first.  First order by gene to match reps
#then order the data by gene to match reps
day14_fecal <- day14_fecal[order(day14_fecal$Genes),]
head(day14_fecal$Genes)
dim(day14_fecal) #2208 4

#need to spread out the samples and populate with the Ct values
day14_fecal_spread <- spread(day14_fecal, Sample, Ct)
dim(day14_fecal_spread) #96 25

#convert all cycles above 28 to NA and remove all rows with NA throughout
day14_fecal_spread[, 2:25][day14_fecal_spread[,2:25] > 28] <-NA
day14_fecal_spread <- day14_fecal_spread[which(!(rowSums(is.na(day14_fecal_spread)) == max(rowSums(is.na(day14_fecal_spread))))),]
dim(day14_fecal_spread) #57 25

#from looking at the samples, 16Sold looks more consistent between reps therefore remove 16Snew
day14_fecal_one16S <- day14_fecal_spread[-c(1,2), ]
dim(day14_fecal_one16S)  #55 25

day14_fecal.gather <- gather(day14_fecal_one16S, key=sample_ID, value=CT, starts_with('P')) %>%
  mutate(sample_gene=paste(sample_ID, Genes, sep='_')) %>%
  group_by(sample_gene) %>% summarise(meanCT=mean(CT),
                                      diff_btw_reps=sort(CT)[2]-sort(CT)[1],
                                      gene=Genes[1],
                                      sample_ID=sample_ID[1],
                                      rep1=CT[1],
                                      rep2=CT[2],
                                      numgood=sum(!is.na(CT)),
                                      onegoodCT=ifelse(numgood==1,CT[which(!is.na(CT))], NA))
write.table(day14_fecal.gather, file="day14_fecal_gather.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)
#now for the plasmidome samples (which are from day 7 also and DNase treated but NOT genomiphi'd)
#First order by gene to match reps
#then order the data by gene to match reps
plasmidome <- plasmidome[order(plasmidome$Genes),]
head(plasmidome$Genes)
dim(plasmidome)  #1920 4
#need to spread out the samples and populate with the Ct values
pldme_spread <- spread(plasmidome, Sample, Ct)
dim(pldme_spread) #96 22
#convert all cycles above 28 to NA and remove all rows with NA for all samples
pldme_spread[, 3:22][pldme_spread[,3:22] > 28] <-NA
pldme_spread <- pldme_spread[which(!(rowSums(is.na(pldme_spread)) == max(rowSums(is.na(pldme_spread))))),]
dim(pldme_spread) #70 22  #re-ran this and it is 58 22?
#from looking at the samples, 16Sold looks more consistent between reps therefore remove 16Snew
pldme_one16S <- pldme_spread[-c(1,2), ]
dim(pldme_one16S)  #rerun is 56 22, not sure what it was the first time
#summarise
pldme.gather <- gather(pldme_one16S, key=sample_ID, value=CT, ends_with('o')) %>%
  mutate(sample_gene=paste(sample_ID, Genes, sep='_')) %>%
  group_by(sample_gene) %>% summarise(meanCT=mean(CT),
                                      diff_btw_reps=sort(CT)[2]-sort(CT)[1],
                                      gene=Genes[1],
                                      sample_ID=sample_ID[1],
                                      rep1=CT[1],
                                      rep2=CT[2],
                                      numgood=sum(!is.na(CT)),
                                      onegoodCT=ifelse(numgood==1,CT[which(!is.na(CT))], NA))
write.table(pldme.gather, file="pldme.gather.tsv", sep="\t", quote=FALSE, col.names=TRUE, row.names=FALSE)

#substitute in the one good CT if there isn't a second one
#remember that ifelse has three details 1) the condition, 2) if yes, 3) if no
day7_fecal.gather$meanCTOG <- ifelse(is.na(day7_fecal.gather$meanCT), day7_fecal.gather$onegoodCT, day7_fecal.gather$meanCT)
day14_fecal.gather$meanCTOG <- ifelse(is.na(day14_fecal.gather$meanCT), day14_fecal.gather$onegoodCT, day14_fecal.gather$meanCT)
pldme.gather$meanCTOG <- ifelse(is.na(pldme.gather$meanCT), pldme.gather$onegoodCT, pldme.gather$meanCT)
#removing the data with 2 bad reps (keeping the good and the one-good-rep)
day7_fecal.goods <- filter(day7_fecal.gather, day7_fecal.gather$numbad!=2)
#note that day7 had numbad instead of numgood!
day14_fecal.goods <- filter(day14_fecal.gather, day14_fecal.gather$numgood !=0)
pldme.goods <- filter(pldme.gather, pldme.gather$numgood !=0)

#fixing the names: note that sub needs pattern, replacement, x
day7_fecal.goods$sample_ID <- sub('X', '', day7_fecal.goods$sample_ID)

Entero <- day7_fecal.goods[grep('E', day7_fecal.goods$sample_ID),]
day7_fecal.goods <- day7_fecal.goods[-grep('E', day7_fecal.goods$sample_ID),]

colnames(day7_fecal.goods)[5]<- 'Pignum' 

#now for day 14 and plasmidome
colnames(day14_fecal.goods)[5]<-"Pignum"
day14_fecal.goods$Pignum <- sub('P', '', day14_fecal.goods$Pignum)
colnames(pldme.goods)[5] <- "Pignum"
pldme.goods$Pignum <- sub('o', '', pldme.goods$Pignum)

#now to select just the genes we want
day7_fecal.goods.spread <- day7_fecal.goods %>% select(gene, Pignum, meanCTOG) %>% spread(gene, meanCTOG)
day14_fecal.goods.spread <- day14_fecal.goods %>% select(gene, Pignum, meanCTOG) %>% spread(gene, meanCTOG)
pldme.goods.spread <- pldme.goods %>% select(gene, Pignum, meanCTOG) %>% spread(gene, meanCTOG)

#bringing in the metadata - this has to be done separately for each run
#bringing in the metadata
unique(day7_fecal.goods.spread$Pignum)
FS1mettt<-read.csv('FS1.qprc.meta.csv', stringsAsFactors=FALSE) 
colnames(FS1mettt)<-c("Pignum", "Treatment")
day7_fecal.metaQPCR<- merge(day7_fecal.goods.spread, FS1mettt, by = 'Pignum')
write.table(day7_fecal.metaQPCR, file="day7_fecals_good_with_metadata.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

FS1mettt$Pignum
pldme.goods.spread$Pignum<-sub("P","",pldme.goods.spread$Pignum)
pldme.goods.spread$Pignum<-as.character(pldme.goods.spread$Pignum)
pldme.metaQPCR<-merge(pldme.goods.spread, FS1mettt, by = "Pignum")
pldme.metaQPCR$Pignum
write.table(pldme.metaQPCR, file="day7_plasmidome_good_with_metadata.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)
dim(pldme.metaQPCR)  #20 38 - different than when first analyzed?  20 is the right number....
#when I analyzed the first time I used the metadata sheet without pig 22 and 25 (those two were problematic in the plasmidome sequencing) 
#however looking at the pldme.metaQPCR file it is clear that they have amplified consistent with the others

#analyzed a third time and got 20 32????
# View(pldme.metaQPCR)

d14_metadata<-read.csv('day14_metadata.csv', stringsAsFactors=FALSE) 
day14_fecal.metaQPCR<- merge(day14_fecal.goods.spread, d14_metadata, by = 'Pignum')
day14_fecal.metaQPCR$Pignum<-as.character(day14_fecal.metaQPCR$Pignum)
# View(day14_fecal.metaQPCR)
write.table(day14_fecal.metaQPCR, file="day14_fecals_good_with_metadata.txt", quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

#still keeping them all separate so that there's no confusion - may combine later

#need to assign the controls for each
day7_fecal_controls <- grep('Control', day7_fecal.metaQPCR$Treatment)
#something went wrong since I have NA instead of InFeed for two of my three figures!
range(day7_fecal_controls) #note this is row numbers, not pignumbers (rows 1-7 and 22-29 are control)

day14_fecal_controls <- grep('Control', day14_fecal.metaQPCR$Treatment)
range(day14_fecal_controls)

pldme_controls <- grep('Control', pldme.metaQPCR$Treatment)
range(pldme_controls)

#Jules' ddCT function can't deal with NA therefore set all NA's to 28 since that was our limit of detection that we chose
day7_fecal.metaQPCR[(is.na(day7_fecal.metaQPCR))]<- 28
day14_fecal.metaQPCR[(is.na(day14_fecal.metaQPCR))]<- 28
pldme.metaQPCR[(is.na(pldme.metaQPCR))]<- 28

#This function has help documentation available
dim(day7_fecal.metaQPCR)  #44 26
day7_fecal_ddct<-qpcr_DDCT(data=day7_fecal.metaQPCR,
                exp_genes_indicies = c(4:25),
                housekeeping_gene_index = 3,
                control_row_indicies = day7_fecal_controls,
                metadata_column_indicies = c(1,26))

dim(day14_fecal.metaQPCR)  #23 31
day14_fecal_ddct<-qpcr_DDCT(data=day14_fecal.metaQPCR,
                           exp_genes_indicies = c(3:30),
                           housekeeping_gene_index = 2,
                           control_row_indicies = day14_fecal_controls,
                           metadata_column_indicies = c(1,31))

dim(pldme.metaQPCR)  #20 38   #20 32 when re-run?
pldme_ddct<-qpcr_DDCT(data=pldme.metaQPCR,
                           exp_genes_indicies = c(3:31),
                           housekeeping_gene_index = 2,
                           control_row_indicies = pldme_controls,
                           metadata_column_indicies = c(1,32))
#now to gather each dataframe
day7_fecal.ddct.gather <- day7_fecal_ddct %>% select(Pignum, Treatment, everything()) %>% gather(key = gene, value = response, -(Pignum:Treatment))

day7_fecal.ddct.gather$Treatment <-factor(day7_fecal.ddct.gather$Treatment, levels=c('Control', 'Inject', 'In-Feed'))

day14_fecal.ddct.gather <- day14_fecal_ddct %>% select(Pignum, Treatment, everything()) %>% gather(key = gene, value = response, -(Pignum:Treatment))

day14_fecal.ddct.gather$Treatment <-factor(day14_fecal.ddct.gather$Treatment, levels=c('Control', 'Inject', 'InFeed'))

pldme.ddct.gather <- pldme_ddct %>% select(Pignum, Treatment, everything()) %>% gather(key = gene, value = response, -(Pignum:Treatment))

pldme.ddct.gather$Treatment <-factor(pldme.ddct.gather$Treatment, levels=c('Control', 'Inject', 'In-Feed'))

day7_fecal.ddct.gather %>% ggplot(aes(x=Treatment, y=response, group=Treatment, fill=Treatment)) +
  geom_boxplot() + facet_wrap(~gene, scales = 'free') + ylab('ddCt') + 
  geom_jitter(shape = 21, width = .15) +
  ggtitle('Day 7 Fecal DeltaDelta CT relative to control mean')

day14_fecal.ddct.gather %>% ggplot(aes(x=Treatment, y=response, group=Treatment, fill=Treatment)) +
  geom_boxplot() + facet_wrap(~gene, scales = 'free') + ylab('ddCt') + 
  geom_jitter(shape = 21, width = .15) +
  ggtitle('Day 14 Fecal DeltaDelta CT relative to control mean')

pldme.ddct.gather %>% ggplot(aes(x=Treatment, y=response, group=Treatment, fill=Treatment)) +
  geom_boxplot() + facet_wrap(~gene, scales = 'free') + ylab('ddCt') + 
  geom_jitter(shape = 21, width = .15) +
  ggtitle('Day 7 Plasmidome DeltaDelta CT relative to control mean')

#can I do a group_by and get mean and sd for each gene by treatment.. no b/c can't figure it out!
day7_fecal.metaQPCR %>% group_by(Treatment) %>% summary()
#write.table(day7_fecal.metaQPCR, file="day7_fecal.metaQPCR", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
#now to determine significant ones and plot those

#Jules' helper function for outputing pairwise.wilcox.test pvalues for each gene
get_pairs <- function(df){
  pp <- pairwise.wilcox.test(df$response, df$Treatment)
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <-paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps
}

colnames(day7_fecal.ddct.gather)
day7_fecal.PW_wilc_per_gene <- day7_fecal.ddct.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data, get_pairs)) %>%
  select(gene, pps) %>% unnest()
colnames(day7_fecal.PW_wilc_per_gene)
range(day7_fecal.PW_wilc_per_gene$`Control_vs_In-Feed`)
colnames(day7_fecal.PW_wilc_per_gene)[3]<-"Control_vs_Feed"
day7_fecal.sigs_feed <- day7_fecal.PW_wilc_per_gene %>% filter(Control_vs_Feed < 0.05)
day7_fecal.sig_inject <- day7_fecal.PW_wilc_per_gene %>% filter(Control_vs_Inject < 0.05)

D7_wilc <- c(day7_fecal.sigs_feed$gene, day7_fecal.sig_inject$gene)

day14_fecal.PW_wilc_per_gene <- day14_fecal.ddct.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data, get_pairs)) %>%
  select(gene, pps) %>% unnest()
colnames(day14_fecal.PW_wilc_per_gene)
range(day14_fecal.PW_wilc_per_gene$`Control_vs_InFeed`)
colnames(day14_fecal.PW_wilc_per_gene)[3]<-"Control_vs_Feed"
day14_fecal.sigs_feed <- day14_fecal.PW_wilc_per_gene %>% filter(Control_vs_Feed < 0.05)
day14_fecal.sig_inject <- day14_fecal.PW_wilc_per_gene %>% filter(Control_vs_Inject < 0.05)
day14_marginal_feed<-day14_fecal.PW_wilc_per_gene %>% filter(Control_vs_Feed < 0.1)

D14_wilc <- c(day14_fecal.sigs_feed$gene, day14_fecal.sig_inject$gene)

# pldme.PW_wilc_per_gene <- pldme.ddct.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data, get_pairs)) %>%
#   select(gene, pps) %>% unnest()
# colnames(pldme.PW_wilc_per_gene)
# range(pldme.PW_wilc_per_gene$`Control_vs_In-Feed`)
# colnames(pldme.PW_wilc_per_gene)[3]<-"Control_vs_Feed"
# pldme.marg_sigs_feed <- pldme.PW_wilc_per_gene %>% filter(Control_vs_Feed < 0.1)
# pldme.marg_sig_inject <- pldme.PW_wilc_per_gene %>% filter(Control_vs_Inject < 0.1)
# #note that aph and intI3 are no longer significant with the addition of the two control pigs that were missed the first time
# write.table(pldme_ddct, file="Plasmidome_ddct.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# write.table(day7_fecal_ddct, file="Day7_fecal_ddct.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
# #not redoing the other figures since Jules' is going to add the LOD lines on them
# 
#modifying Jules' function to do a t-test instead of wilcox
get_paired_t <- function(df){
  pp <- pairwise.t.test(df$response, df$Treatment, p.adjust.method="fdr")
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <-paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps
}
bonferroni_day7_fecal.PW_ttest_per_gene <- day7_fecal.ddct.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data, get_paired_t)) %>%
  select(gene, pps) %>% unnest()
bonferroni_day14_fecal.PW_ttest_per_gene <- day14_fecal.ddct.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data, get_paired_t)) %>%
  select(gene, pps) %>% unnest()
bonferroni_pldme.PW_ttest_per_gene <- pldme.ddct.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data, get_paired_t)) %>%
  select(gene, pps) %>% unnest()

colnames(bonferroni_day7_fecal.PW_ttest_per_gene)[3]<-"Control_vs_Feed"
sig.day7.fecal.PW.ttest <- subset(bonferroni_day7_fecal.PW_ttest_per_gene, Control_vs_Inject < 0.05 | Control_vs_Feed < 0.05)

colnames(bonferroni_day14_fecal.PW_ttest_per_gene)[3]<-"Control_vs_Feed"
sig.day14.fecal.PW.ttest <- subset(bonferroni_day14_fecal.PW_ttest_per_gene, Control_vs_Inject < 0.05 | Control_vs_Feed < 0.05)

# colnames(bonferroni_pldme.PW_ttest_per_gene)[3]<-"Control_vs_Feed"
# sig.pldme.PW.ttest <- subset(bonferroni_pldme.PW_ttest_per_gene, Control_vs_Inject < 0.05 | Control_vs_Feed < 0.05)

sig_day7<-sig.day7.fecal.PW.ttest$gene
day7_fecal.ddct.gather %>% subset(gene %in% sig_day7) %>%
  ggplot(aes(x=Treatment, y=response, group=Treatment, fill=Treatment)) +
  geom_boxplot(outlier.color = NA) + facet_wrap(~gene, scales = 'free') + ylab('log2 fold change') + 
  geom_jitter(shape = 21, width = .15) + theme_bw() #+
  # ggtitle('Day 7 Significant Genes delta delta CT relative to control mean')
#shape = 21 must be open circles

sig_day14<-sig.day14.fecal.PW.ttest$gene
day14_fecal.ddct.gather %>% subset(gene %in% sig_day14) %>%
  ggplot(aes(x=Treatment, y=response, group=Treatment, fill=Treatment)) +
  geom_boxplot(outlier.color = NA) + facet_wrap(~gene, scales = 'free') + ylab('log2 fold change') + 
  geom_jitter(shape = 21, width = .15) + theme_bw()


#################

day7_fecal.ddct.gather$Treatment <- as.character(day7_fecal.ddct.gather$Treatment)
day14_fecal.ddct.gather$Treatment <- as.character(day14_fecal.ddct.gather$Treatment)

day7_fecal.ddct.gather$Treatment[day7_fecal.ddct.gather$Treatment == 'In-Feed'] <- 'Feed'
day14_fecal.ddct.gather$Treatment[day14_fecal.ddct.gather$Treatment == 'InFeed'] <- 'Feed'

day7_fecal.ddct.gather$Treatment[day7_fecal.ddct.gather$Treatment == 'Control'] <- 'NM'
day14_fecal.ddct.gather$Treatment[day14_fecal.ddct.gather$Treatment == 'Control'] <- 'NM'

day7_fecal.ddct.gather$Treatment <- factor(day7_fecal.ddct.gather$Treatment, levels = c('NM', 'Inject', 'Feed'))
day14_fecal.ddct.gather$Treatment <- factor(day14_fecal.ddct.gather$Treatment, levels = c('NM', 'Inject', 'Feed'))

day7_fecal.ddct.gather %>% subset(gene %in% D7_wilc) %>%
  ggplot(aes(x=Treatment, y=response, group=Treatment, fill=Treatment)) +
  geom_boxplot(outlier.color = NA) + facet_wrap(~gene, scales = 'free') + ylab('log2 fold change') + 
  geom_jitter(shape = 21, width = .15) + theme_bw() #+
# ggtitle('Day 7 Significant Genes delta delta CT relative to control mean')
#shape = 21 must be open circles

day14_fecal.ddct.gather %>% subset(gene %in% D14_wilc) %>%
  ggplot(aes(x=Treatment, y=response, group=Treatment, fill=Treatment)) +
  geom_boxplot(outlier.color = NA) + facet_wrap(~gene, scales = 'free') + ylab('log2 fold change') + 
  geom_jitter(shape = 21, width = .15) + theme_bw()















###################




#is the data normally distributed?
shapiro.test(day7_fecal.ddct.gather$response)
aph <- day7_fecal.ddct.gather %>% subset(gene=="aph2'-Id_104")
shapiro.test(aph$response)
aph_pldme <- pldme.ddct.gather %>% subset(gene=="aph2-Id-104")
shapiro.test(aph_pldme$response)
bla<-pldme.ddct.gather %>% subset(gene=="blaROB-41")
shapiro.test(bla$response)
hist(aph_pldme$response)

qqnorm(day7_fecal.ddct.gather$response)
qqnorm(day14_fecal.ddct.gather$response)
qqnorm(pldme.ddct.gather$response)
hist(day7_fecal.ddct.gather$response)
hist(day14_fecal.ddct.gather$response)
hist(pldme.ddct.gather$response)
 #histograms for day 14 and pldme are each skewed in different directions.  Day 14 has higher +1 range compared to -1.  Pldme is the opposite.
#day 7 fecal has a normal distribution, with much higher numbers in teh -1 to +1 range and then dropping off.

#Jules suggested trying the fdr correction instead of holm for the data to see if it makes a difference (still using wilcox)
get_pairs_wilcox_fdr <- function(df){
  pp <- pairwise.wilcox.test(df$response, df$Treatment, p.adjust.method="fdr")
  ppdf <- as.data.frame(pp$p.value)
  ps <- data.frame(matrix(c(pp$p.value), nrow = 1))
  names(ps) <-paste(c(rep(names(ppdf), each = nrow(ppdf))), "_vs_", rep(rownames(ppdf), ncol(ppdf)), sep = "")
  ps
}

colnames(day7_fecal.ddct.gather)
day7_fecal.PW_wilc_fdr <- day7_fecal.ddct.gather %>% group_by(gene) %>% nest() %>% mutate(pps = map(data, get_pairs_wilcox_fdr)) %>%
  select(gene, pps) %>% unnest()
pldme.PW_wilc_fdr <- pldme.ddct.gather %>% group_by(gene) %>% nest() %>% mutate(pps=map(data, get_pairs_wilcox_fdr)) %>%
  select(gene, pps) %>% unnest()

#need to make a table with all the meanCT values
# View(day7_fecal.goods)
day7_meanCTOG<- day7_fecal.goods %>% group_by(gene) %>% summarize(meanCT_gene=mean(meanCTOG), sd_gene=sd(meanCTOG)) %>% select(gene, meanCT_gene, sd_gene) %>% unique()
day7_meanCTOG<-day7_meanCTOG[-1, ]
#added standard deviation but haven't included it in final tables
colnames(day7_meanCTOG)<-c("Gene","Day 7 MeanCT")
day7_meanCTOG$Gene<-gsub("'","",day7_meanCTOG$Gene)
day7_meanCTOG$Gene<-gsub("-","_",day7_meanCTOG$Gene)

day14_meanCTOG<- day14_fecal.goods %>% group_by(gene) %>% summarize(meanCT_gene=mean(meanCTOG)) %>% select(gene, meanCT_gene) %>% unique()
head(day14_meanCTOG)
colnames(day14_meanCTOG)<-c("Gene","Day 14 MeanCT")
day14_meanCTOG$Gene<-gsub("-", "_", day14_meanCTOG$Gene)

pldme_meanCTOG<- pldme.goods %>% group_by(gene) %>% summarize(meanCT_gene=mean(meanCTOG), sd_gene=sd(meanCTOG)) %>% select(gene, meanCT_gene, sd_gene) %>% unique()
head(pldme_meanCTOG)
#added in sd later but didn't redo the tables to include it
colnames(pldme_meanCTOG)<-c("Gene","Plasmidome MeanCT")
pldme_meanCTOG$Gene<-gsub("-", "_", pldme_meanCTOG$Gene)

fecal_meanCTOG<-full_join(day7_meanCTOG, day14_meanCTOG, by = "Gene")
all_meanCTOG<-full_join(fecal_meanCTOG, pldme_meanCTOG, by = "Gene")
write.table(all_meanCTOG, file="FS1_meanCTOG.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)
#note that if we alter the N/A it messes up the decimals
#if we want to use this we can add in standard deviation and adjust decimals?
