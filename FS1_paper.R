


###### ABX conc. data  #########



setwd('~/FS1/reanal/')

library(tidyverse)
library(reshape2)

# functions used:
# grep()
# gsub()
# colnames()
# na.exclude()
# paste()


tis <- read.csv('HansResults_tissues.csv', stringsAsFactors = FALSE)  # reads in data, already cleaned a little
meta <- read.csv('hans_meta.csv', stringsAsFactors = FALSE)           # reads in some metadata

inj <- meta$Pig.ID.Numbers[grep('IM', meta$Control.treatment)]        # creates a vector of pig numbers that have 'IM' in their treatment column
oral <- meta$Pig.ID.Numbers[grep('In-feed', meta$Control.treatment)]    # creates a vector of pig numbers that have 'In-feed' in their treatment column
control <- meta$Pig.ID.Numbers[grep('Control', meta$Control.treatment)] # creates a vector of pig numbers that have 'Control' in their treatment column


tis$Time.Point <- as.numeric(gsub('D', '', tis$Time.Point))             # replaces the 'D' in the time column with '' (nothing)

tis.melt <- melt(tis, id.vars = c(1,2), measure.vars = c(3:6))          # converts to long dataframe format for easy plotting
tis.melt                                                                # just checking on the new long dataframe
colnames(tis.melt) <- c('pig', 'day', 'tissue', 'concentration')        # changing the column names to something I like
tis.melt$concentration <- gsub('NF', 0, tis.melt$concentration)         # replaces 'NF' with 0
tis.melt$concentration <- gsub('<', '', tis.melt$concentration)         # replaces '<' with '' (nothing)


tis.melt$concentration <- as.numeric(tis.melt$concentration)            # forces the concentration column to be numeric 
tis.melt <- na.exclude(tis.melt)                                        # removes NAs which were introduced in the previous line

tis.melt$treatment <- NA                                                # creates a new column 'treatment'
tis.melt$treatment[tis.melt$pig %in% control] <- 'control'              # assigns 'control' to pigs whose numbers are in the control vector we generated on line 11
tis.melt$treatment[tis.melt$pig %in% inj] <- 'inject'                   # assigns 'inject' to pigs whose numbers are in the control vector we generated on line 9
tis.melt$treatment[tis.melt$pig %in% oral] <- 'oral'                    # assigns 'oral' to pigs whose numbers are in the control vector we generated on line 10


tis.melt$dayXtreat <- paste(tis.melt$day, tis.melt$treatment, sep = '_') # creates a dayXtreatment column, it's just the day and treatment columns pasted together




tis.melt$day <- factor(tis.melt$day)                                     # makes the day column a factor (this means it is categorical data, not continuous data)


# this section uses tidyverse type pipes, filters, and ggplot to produce boxplots


tis.melt %>% filter(tissue == 'Serum') %>%                                 
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=treatment)) +
  geom_boxplot() + ylab('concentration ng/mL') + ggtitle('Oxytet concentrations: Serum')

tis.melt %>% filter(tissue == 'feces') %>%
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=treatment)) +
  geom_boxplot() + ylab('concentration ng/mL') + ggtitle('Oxytet concentrations: Feces')

tis.melt %>% filter(tissue == 'Nasal') %>%
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=treatment)) +
  geom_boxplot() + ylab('concentration ng/mL') + ggtitle('Oxytet concentrations: Nasal wash')


tis.melt %>% filter(tissue == 'Ileum') %>%
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=treatment)) +
  geom_boxplot() + ylab('concentration ng/mL') + ggtitle('Ileum')

tis.melt <- tis.melt[!(tis.melt$concentration >129000 & tis.melt$tissue == 'Ileum'),]  # removed one crazy high outlier

tis.melt %>% filter(tissue == 'Ileum') %>%
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=treatment)) +
  geom_boxplot() + ylab('concentration ng/mL') + ggtitle('Oxytet concentrations: Ileum')


tis.melt$group <- paste(tis.melt$pig, tis.melt$day, tis.melt$tissue, sep = '_')

feces.abx <- tis.melt %>% filter(tissue == 'feces')

#  Duplicated Pig numbers are an issue here. This is some fallout from when pig 50 was taken
# instead of pig 5 for a D4 necropsy.

#tis.melt %>% spread(tissue, concentration)


#tis.spread <- spread(tis.melt, key = tissue, value = concentration)

tis.melt[c(381,514),]
tis.melt[c(20,249),]
tis.melt[c(660,718),]
tis.melt[c(661,719),]

# this should fix the problem

tis.melt[c(381),]$pig <- 5
tis.melt[c(20),]$pig <- 5
tis.melt[c(718),]$pig <- 50
tis.melt[c(719),]$pig <- 50

filter(tis.melt, pig %in% c(5,50) & tissue %in% c('Nasal', 'Serum', 'feces') & day %in% c(4, 11, 14))

ileum.abx <- tis.melt %>% filter(tissue == 'Ileum') 
feces.abx <- tis.melt %>% filter(tissue == 'feces') 
Nasal.abx <- tis.melt %>% filter(tissue == 'Nasal') 
serum.abx <- tis.melt %>% filter(tissue == 'Serum') 

###################### 16S data ###############





setwd('~/FS1/V4/')

### this script requires that you have mothur in your system PATH  ###
library(vegan)
library(dplyr)
library(ggplot2)

library(Hmisc)

library(geomnet)
library(funfuns)

gather_nodes <- function(x, typ=NA){
  x <- as.data.frame(colnames(x))
  colnames(x)[1] <- 'node'
  x$type <- typ
  return(x)
  
}

########## functions ###########


rcorr_to_ggnet <- function(rcorr.list, pcut=0.05, spearcut=0.6){
  pval <- as.data.frame(rcorr.list$P)
  sim <- as.data.frame(rcorr.list$r)
  
  sim$to <- rownames(sim)
  pval$to <- rownames(pval)
  
  pval <- melt(pval, id.vars = 'to')
  sim <- melt(sim, id.vars = 'to')
  
  colnames(pval) <- c('to', 'from', 'pval')
  colnames(sim) <- c('to', 'from', 'spear')
  
  pval <- pval[,c(2,1,3)]
  sim <- sim[,c(2,1,3)]
  
  sigcor <- merge(pval,sim )
  sigcor <- na.exclude(sigcor)
  sigcor <- sigcor[sigcor$pval < pcut,]
  sigcor <- sigcor[sigcor$spear > spearcut,]
  
  return(sigcor)
  
}


extract_mothur_tax <- function(filename){
  
  tax <- read.table(filename, header = TRUE, stringsAsFactors = FALSE)
  tna <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
  tax[tna] <- NA
  rowcount <- 1
  
  for (i in strsplit(tax$Taxonomy, split = ';')){
    levelcount <- 1
    for (x in i){
      tax[rowcount,3+levelcount] <- x
      levelcount <- levelcount + 1
      #print(x)
    }
    rowcount <- rowcount+1
  }
  return(tax)
}

otu_tax_labels <- function(mothur_tax){
  tax <- mothur_tax
  swip <- tax$OTU
  names(swip) <- paste(tax$Genus, tax$OTU, sep = ':')
  names(swip)[grepl("NA", names(swip))] <- paste(tax$Family[grepl("NA", names(swip))], tax$OTU[grepl("NA", names(swip))], sep = ':')
  swap <- names(swip)
  names(swap) <- swip
  return(swap)
}


read.mothur.ovas <- function(filename){
  mothurfile <- scan(file = filename, what = "character", sep = "\t")
  name <- filename
  groups <- c()
  fvalues <- c()
  pvalues <- c()
  groupcount<- 1
  fvaluecount <- 17
  pvaluecount <- 18
  for (x in mothurfile){
    groups <- c(groups, mothurfile[groupcount])
    groupcount <- groupcount +18
    fvalues <- c(fvalues, mothurfile[fvaluecount])
    fvaluecount <- fvaluecount +18
    pvalues <- c(pvalues, mothurfile[pvaluecount])
    pvaluecount <- pvaluecount +18
  } 
  results.table <- data.frame(Groups=groups, Fvalue=fvalues, Pvalue=pvalues, stringsAsFactors = FALSE)  # results table, but has empty cells because of lengths of vectors?
  results.table <- results.table[1:(length(mothurfile)/18),]   # cleans up empty cells
  results.table$Pvalue <- gsub("p-value: ", "", results.table$Pvalue)
  results.table$Pvalue <- gsub("\\*", "", results.table$Pvalue)  # may need to remove
  results.table$Pvalue <- gsub("<", "", results.table$Pvalue)    # may need to remove
  results.table$Pvalue <- as.numeric(results.table$Pvalue)
  results.table$Fvalue <- as.numeric(results.table$Fvalue)
  
  
  results.table
}


pairwise.adonis <- function(x,factors, sim.function = 'vegdist', sim.method = 'bray', p.adjust.m ='holm')
{
  library(vegan)
  
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])], permutations = 999);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  sig = c(rep('',length(p.adjusted)))
  sig[p.adjusted <= 0.05] <-'.'
  sig[p.adjusted <= 0.01] <-'*'
  sig[p.adjusted <= 0.001] <-'**'
  sig[p.adjusted <= 0.0001] <-'***'
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted,sig)
  print("Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1")
  return(pairw.res)
  
} 


veganCovEllipse <- function (cov, center = c(0,0), scale = 1, npoints = 100){
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}



############ end functions ############



####### Metadata extraction from sample names #######

shared <- read.table('FS1.V4.shared', header = TRUE)   # this is the OTU table output by mothur

meta <- shared[,c(1,2,3)]
meta$day <- as.numeric(gsub('.*D([:digit:]*)','\\1', meta$Group))
meta$pig_number <- as.numeric(sub('X.P(.*)N.*','\\1', meta$Group))
meta$pen <- sub('X.P.*N(..).*D.*', '\\1', meta$Group)
meta$room <- sub('(.).', '\\1', meta$pen )
meta$experiment <- sub('(X.)P[:digit:]*.*', '\\1', meta$Group)
meta$tissue <- sub('X.P.*N..(.*)D.*', '\\1', meta$Group)
meta$tissue[meta$tissue == 'F'] <- 'Feces'
meta$tissue[meta$tissue == 'C'] <- 'Colon'
meta$tissue[meta$tissue == 'I'] <- 'Ileum'

meta <- meta[which(meta$experiment == 'X1'),]

OTC_nums <- c(15:21,36:42, 59:66)
control_nums <- c(1:7,22:28, 43:50, 67:73)
inj_nums <- c(8:14, 29:35, 51:58)
RPS_nums <- c(74:80)

controls <- meta$pig_number %in% control_nums
orals <- meta$pig_number %in% OTC_nums
injections <- meta$pig_number %in% inj_nums


meta$treatment[controls] <- 'control'
meta$treatment[orals] <- 'oral'
meta$treatment[injections] <- 'inject'

meta$Group # checking on the sample names....

# this is a metadata table we generated during the experiment

master <- read.csv("FS1master.csv", header = TRUE, stringsAsFactors = TRUE) 

colnames(master)[1] <- 'pig_number'

mastermeta <- merge(meta, master, by = 'pig_number')
colnames(mastermeta)[12] <- 'birth_room'
colnames(mastermeta)[15] <- 'initial_weight'
mastermeta <- mastermeta[,-c(2,4,16,18,19,21)]

day <- paste('D', mastermeta$day, sep = '')
mastermeta$set <- paste(day, mastermeta$tissue, mastermeta$treatment, sep = '_')

######### Vegan Style #############


rownames(mastermeta) <- mastermeta$Group

#mastermeta <- mastermeta[,-2]
mastermeta
rownames(shared) <- shared[,2]
shared <- shared[,-c(1,2,3)]
shared <- shared[rownames(shared) %in% rownames(mastermeta),]

tax[grep('Mitochondria', tax$Genus),]
shared <- shared[,-grep('Mitochondria', tax$Genus)]


shared <- shared[rowSums(shared) > 2000,] # removes samples with less than 2000 reads
shared <- shared[,colSums(shared) > 10] # removes OTUs occuring less than 10 times globally
mastermeta <- mastermeta[rownames(mastermeta) %in% rownames(shared),]

# 314 samples
# 2616 OTUs why does this not match the phyloseq object later on?

all_groups <- mastermeta %>% group_by(set) %>% summarise(no_rows = length(set)) # the number of samples falling into each set
mastermeta <- mastermeta[mastermeta$set != 'D7_Ileum_inject',] # Only 2 samples in the D7_Ileum_inject group, removed both.
write.csv(all_groups, 'pigs_per_group.csv')
shared <- shared[rownames(shared) %in% rownames(mastermeta),]

min(rowSums(shared))

share.rare <- rrarefy(shared, min(rowSums(shared)))  # rarrefied to 2202 sequences for alpha and beta diversity calcs

rownames(share.rare) == rownames(mastermeta)  # these two dataframes are in different orders, we need to fix this

share.rare <- share.rare[match(rownames(mastermeta), rownames(share.rare)),]
rownames(share.rare) == rownames(mastermeta)  # this confirms that they are in the same order.

share.rare2 <- share.rare # for boxplots later on
write.table(mastermeta, 'FS1.mastermeta.txt', quote = FALSE, sep = '\t')

#read.table('FS1.mastermeta.txt', header = TRUE)
########## alpha diversity?  ##########

rowSums(share.rare)

mastermeta$shannon <- diversity(share.rare)


mastermeta %>% filter(tissue=='Feces') %>%  
  ggplot(aes(x=day, y=shannon, group=set, fill=treatment)) + geom_boxplot() + ggtitle('Alpha diversity: Feces') + ylab("Shannon Diversity Index")

mastermeta %>% filter(tissue=='Colon') %>%  
  ggplot(aes(x=day, y=shannon, group=set, fill=treatment)) + geom_boxplot() + ggtitle('Alpha diversity: Colon')+ ylab("Shannon Diversity Index")

mastermeta %>% filter(tissue=='Ileum') %>% 
  ggplot(aes(x=day, y=shannon, group=set, fill=treatment)) + geom_boxplot() + ggtitle('Alpha diversity: Ileum')+ ylab("Shannon Diversity Index")


### Only doing feces here ###

mastermeta.feces <- mastermeta[mastermeta$tissue == 'Feces',]

mastermeta.feces <- mastermeta.feces[mastermeta.feces$experiment == 'X1',]
mastermeta.feces$set <- factor(mastermeta.feces$set)
share.rare <- share.rare[rownames(share.rare) %in% rownames(mastermeta.feces),]
rownames(share.rare) == rownames(mastermeta.feces)

#group_numbers <- mastermeta.feces %>% group_by(set) %>% summarise(no_rows = length(set))

#colnames(group_numbers) <- c('set', 'number_samples')
#write.table(group_numbers, 'group_nums.txt', sep = '\t', quote = FALSE, row.names = FALSE)
############ ################

PWadon <- pairwise.adonis(share.rare, mastermeta.feces$set)



fecesVSfeces <- grep(".*_Feces_.* vs .*_Feces_.*", PWadon$pairs)
# colonVScolon <- grep(".*_Colon_.* vs .*_Colon_.*", PWadon$pairs)
# ileumVSileum <- grep(".*_Ileum_.* vs .*_Ileum_.*", PWadon$pairs)
# 
# goodones <- c(colonVScolon, ileumVSileum)
# PWadon.tissue <- PWadon[goodones,]
# 
# D4t <- grep("D4_.*_.* vs D4_.*_.*", PWadon.tissue$pairs)
# D7t <- grep("D7_.*_.* vs D7_.*_.*", PWadon.tissue$pairs)
# D14t <- grep("D14_.*_.* vs D14_.*_.*", PWadon.tissue$pairs)
# 
# PWadon.tissue <- PWadon.tissue[c(D4t, D7t, D14t),]

D0 <- grep("D0_Feces_.* vs D0_Feces_.*", PWadon$pairs)
D4 <- grep("D4_Feces_.* vs D4_Feces_.*", PWadon$pairs)
D7 <- grep("D7_Feces_.* vs D7_Feces_.*", PWadon$pairs)
D11 <- grep("D11_Feces_.* vs D11_Feces_.*", PWadon$pairs)
D14 <- grep("D14_Feces_.* vs D14_Feces_.*", PWadon$pairs)

PWadon.good <- PWadon[c(D0, D4, D7, D11, D14),]

write.table(PWadon.good, file = 'FS1_all_PERMANOVA.txt', quote = FALSE, row.names = FALSE, sep = '\t')



######  
time.inj <- grep(".*_Feces_inject vs .*_Feces_inject", PWadon$pairs)
time.oral <- grep(".*_Feces_oral vs .*_Feces_oral", PWadon$pairs)
time.control <- grep(".*_Feces_control vs .*_Feces_control", PWadon$pairs)

PWadon.time <- PWadon[c(time.inj, time.control, time.oral),]
write.table(PWadon.time, file = 'PWadonis.time', quote = FALSE, row.names = FALSE, sep = '\t')


################  Just plotting adonis Fval to control vs day  ##################

PWadon.tocont <- PWadon.good[grep('control vs ', PWadon.good$pairs),]
PWadon.tocont$day <- gsub('D([0-9]+)_.*', '\\1', PWadon.tocont$pairs)
PWadon.tocont$treatment <- gsub('D[0-9]+_Feces_control vs D[0-9]+_Feces_([A-Za-z]+)', '\\1', PWadon.tocont$pairs)
PWadon.tocont$group <- gsub('D[0-9]+_Feces_control vs (D[0-9]+_Feces_[A-Za-z]+)', '\\1', PWadon.tocont$pairs)
#PWadon.tocont <- PWadon.tocont[order(as.numeric(PWadon.tocont$day)),]
PWadon.tocont$day <- factor(PWadon.tocont$day, levels = c('0', '4', '7','11', '14'))


p2 <- ggplot(PWadon.tocont, aes(x=day, y=F.Model, color=treatment)) +
  geom_path(aes(group=treatment), size=1.3) + 
  scale_fill_manual(values = c('#00BA38', '#619CFF')) + 
  geom_vline(xintercept = 3, color='purple', size = 1, alpha=.75) +
  geom_label(aes(label=p.value, fill=treatment), color='black', show.legend = FALSE) +
  scale_color_manual(values = c('#00BA38', '#619CFF'))
p2 + ggtitle('Route of antibiotic administration affects the magnitude of fecal microbiota community disturbance') + 
  ylab('PERMANOVA F vs control\n(difference relative to control)')


##### feces ordination ######

mastermeta.feces$Group <- rownames(mastermeta.feces)


rownames(mastermeta.feces) == rownames(share.rare)

NMDS <- NMDS_ellipse(metadata = mastermeta.feces, OTU_table = share.rare, grouping_set = 'set')



feces.metanmds <- NMDS[[1]]
df_ell <- NMDS[[2]]

df_ell$day <- gsub('(D[0-9]+)_([A-Za-z]+)_([A-Za-z]+)', '\\1', df_ell$group)
df_ell$tissue <- gsub('(D[0-9]+)_([A-Za-z]+)_([A-Za-z]+)', '\\2', df_ell$group)
df_ell$treatment <- gsub('(D[0-9]+)_([A-Za-z]+)_([A-Za-z]+)', '\\3', df_ell$group)

df_ell$day <- as.numeric(gsub('D', '', df_ell$day))
feces.metanmds$day
feces.metanmds$set



# filter(feces.metanmds, day ==0) %>% ggplot(aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color=set), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
#   geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set), alpha=.5) +
#   geom_path(data=filter(df_ell, day == 0), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
#   ggtitle('Day 0', subtitle = 'NMDS projection of Bray-Curtis community structure similarities') + geom_text(aes(label=pig_number))
# 
# filter(feces.metanmds, day ==4) %>% ggplot(aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color=set), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
#   geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set), alpha=.5) + 
#   geom_path(data=filter(df_ell, day == 4), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
#   ggtitle('Day 4', subtitle = 'NMDS projection of Bray-Curtis community structure similarities') + geom_text(aes(label=pig_number))
# 
# filter(feces.metanmds, day ==7) %>% ggplot(aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color=set), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
#   geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set), alpha=.5) + 
#   geom_path(data=filter(df_ell, day == 7), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
#   ggtitle('Day 7', subtitle = 'NMDS projection of Bray-Curtis community structure similarities') + geom_text(aes(label=pig_number))
# 
# filter(feces.metanmds, day ==11) %>% ggplot(aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color=set), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
#   geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set), alpha=.5) + 
#   geom_path(data=filter(df_ell, day == 11), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
#   ggtitle('Day 11', subtitle = 'NMDS projection of Bray-Curtis community structure similarities') + geom_text(aes(label=pig_number))
# 
# 
# filter(feces.metanmds, day ==14) %>% ggplot(aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color=set), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
#   geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set), alpha=.5) + 
#   geom_path(data=filter(df_ell, day == 14), aes(x=NMDS1, y=NMDS2, color=group), size=1.25) +
#   ggtitle('Day 14', subtitle = 'NMDS visualization of Bray-Curtis community structure similarities') + geom_text(aes(label=pig_number))
# 
# 
# 
# ####paths #####
# 
# feces.metanmds$day
# 
# feces.metanmds <- feces.metanmds[order(feces.metanmds$day),]
# #feces.metanmds$pigXday <- paste(feces.metanmds$pig_number, feces.metanmds$day, sep = '')
# filter(feces.metanmds, Necropsy.day == 14 & day != 0) %>% 
#   ggplot(aes(x=MDS1, y=MDS2)) + geom_path(aes(group=pig_number, color=treatment)) + 
#   geom_text(aes(label=day)) + facet_grid(~treatment)
# 
# df_ell2 <- df_ell
# df_ell2$set <- df_ell2$group
# 
# ggplot(feces.metanmds, aes(x=MDS1, y=MDS2)) + 
#   geom_point(aes(color=treatment), size=2) + #geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=group), size=1.25)+
#   geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) + 
#   geom_path(data=df_ell2, aes(x=NMDS1, y=NMDS2, color=treatment), size=1.25) +
#   facet_wrap(~day, scales = 'free', ncol = 5, nrow=1) + ggtitle('Community similarities over time', subtitle = 'NMDS projection of Bray-Curtis similarities') + 
#   labs(caption='ordination stress: 0.185')
# 

##### Trying to mimick phil's grey stuff ####
day0 <- filter(feces.metanmds, day == 0)
day4 <- filter(feces.metanmds, day == 4)
day7 <- filter(feces.metanmds, day == 7)
day11 <- filter(feces.metanmds, day == 11)
day14 <- filter(feces.metanmds, day == 14)

df_ell0 <- filter(df_ell, day == 0)
df_ell4 <- filter(df_ell, day == 4)
df_ell7 <- filter(df_ell, day == 7)
df_ell11 <- filter(df_ell, day == 11)
df_ell14 <- filter(df_ell, day == 14)

feces.metanmds2 <- feces.metanmds




ggplot(feces.metanmds, aes(x=MDS1, y=MDS2)) +  annotate(x=feces.metanmds2$MDS1, y=feces.metanmds2$MDS2, color='grey57', geom = 'point')+
  geom_path(data = df_ell, aes(x=NMDS1, y=NMDS2, color=treatment), size=1.25) + 
  geom_point(aes(color = treatment), size = 2) + 
  geom_segment(aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + facet_wrap(~day, nrow = 1, ncol = 5) + 
  ggtitle('NMDS visualization of Bray-Curtis similarities', subtitle = 'ordination stress = 0.179')



D0 <- ggplot(feces.metanmds, aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') +
  geom_path(data = df_ell0, aes(x=NMDS1, y=NMDS2, color=treatment), size=1.25) + 
  geom_point(data=day0, aes(color = treatment), size = 2) + 
  geom_segment(data = day0, aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) + 
  annotate(geom='text', x=-.6, y=.25, label='Day 0', size = 8) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + 
  ggtitle('NMDS visualization of Bray-Curtis similarities', subtitle = 'ordination stress = 0.179')

D0

D4 <- ggplot(feces.metanmds, aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') +
  geom_path(data = df_ell4, aes(x=NMDS1, y=NMDS2, color=treatment), size=1.25) + 
  geom_point(data=day4, aes(color = treatment), size = 2) + 
  geom_segment(data = day4, aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) + 
  annotate(geom='text', x=-.6, y=.25, label='Day 4', size = 8) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank())+
  ggtitle('NMDS visualization of Bray-Curtis similarities', subtitle = 'ordination stress = 0.179')

D4

D7 <- ggplot(feces.metanmds, aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') +
  geom_path(data = df_ell7, aes(x=NMDS1, y=NMDS2, color=treatment), size=1.25) + 
  geom_point(data=day7, aes(color = treatment), size = 2) + 
  geom_segment(data = day7, aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) +
  annotate(geom='text', x=-.6, y=.25, label='Day 7', size = 8) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + 
  ggtitle('NMDS visualization of Bray-Curtis similarities', subtitle = 'ordination stress = 0.179')

D7

D11 <- ggplot(feces.metanmds, aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') +
  geom_path(data = df_ell11, aes(x=NMDS1, y=NMDS2, color=treatment), size=1.25) + 
  geom_point(data=day11, aes(color = treatment), size = 2) + 
  geom_segment(data = day11, aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) +
  annotate(geom='text', x=-.6, y=.25, label='Day 11', size = 8) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + 
  ggtitle('NMDS visualization of Bray-Curtis similarities', subtitle = 'ordination stress = 0.179')

D11

D14 <- ggplot(feces.metanmds, aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') +
  geom_path(data = df_ell14, aes(x=NMDS1, y=NMDS2, color=treatment), size=1.25) + 
  geom_point(data=day14, aes(color = treatment), size = 2) + 
  geom_segment(data = day14, aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=treatment), alpha=.5) +
  annotate(geom='text', x=-.583, y=.25, label='Day 14', size = 8) + 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        panel.border = element_rect(fill = NA, color = 'grey57'),
        axis.line = element_blank()) + 
  ggtitle('NMDS visualization of Bray-Curtis similarities', subtitle = 'ordination stress = 0.179')

D14

#F8766D
#00BA38
#619CFF

#  + scale_colour_manual(values = c('#F8766D', '#00BA38', '#619CFF'))
# ggplot_build(D14)


##### Tissue Ordination  #####
rownames(share.rare2) == rownames(mastermeta)

mastermeta.tissue <- mastermeta[mastermeta$tissue != 'Feces',]

share.rare.tissue <- share.rare2[rownames(share.rare2) %in% rownames(mastermeta.tissue),]
rownames(share.rare.tissue) == rownames(mastermeta.tissue)


NMDS_all <- NMDS_ellipse(metadata = mastermeta.tissue, OTU_table = share.rare.tissue, grouping_set = 'set')

all.metanmds <- NMDS_all[[1]]
all.df_ell <- NMDS_all[[2]]

all.df_ell$day <- gsub('(D[0-9]+)_([A-Za-z]+)_([A-Za-z]+)', '\\1', all.df_ell$group)
all.df_ell$tissue <- gsub('(D[0-9]+)_([A-Za-z]+)_([A-Za-z]+)', '\\2', all.df_ell$group)
all.df_ell$treatment <- gsub('(D[0-9]+)_([A-Za-z]+)_([A-Za-z]+)', '\\3', all.df_ell$group)

all.df_ell$day <- as.numeric(gsub('D', '', all.df_ell$day))
all.metanmds$day
all.metanmds$set


all.df_ell_4 <- filter(all.df_ell, day == 4)
all.df_ell_7 <- filter(all.df_ell, day == 7)
all.df_ell_14 <- filter(all.df_ell, day == 14)

all.metanmds_4 <- filter(all.metanmds, day == 4)
all.metanmds_7 <- filter(all.metanmds, day == 7)
all.metanmds_14 <- filter(all.metanmds, day == 14)

##########


ggplot(data = filter(all.metanmds, tissue =='Colon'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_4, tissue == 'Colon'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_4, tissue == 'Colon'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_4, tissue == 'Colon'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Colon, Day 4', subtitle = 'stress: 0.169')

ggplot(data = filter(all.metanmds, tissue =='Colon'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_7, tissue == 'Colon'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_7, tissue == 'Colon'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_7, tissue == 'Colon'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Colon, Day 7', subtitle = 'stress: 0.169')

ggplot(data = filter(all.metanmds, tissue =='Colon'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_14, tissue == 'Colon'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_14, tissue == 'Colon'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_14, tissue == 'Colon'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Colon, Day 14', subtitle = 'stress: 0.169')

########

ggplot(data = filter(all.metanmds, tissue =='Ileum'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_4, tissue == 'Ileum'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_4, tissue == 'Ileum'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_4, tissue == 'Ileum'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Ileum, Day 4', subtitle = 'stress: 0.169')

ggplot(data = filter(all.metanmds, tissue =='Ileum'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_7, tissue == 'Ileum'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_7, tissue == 'Ileum'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_7, tissue == 'Ileum'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Ileum, Day 7', subtitle = 'stress: 0.169')

ggplot(data = filter(all.metanmds, tissue =='Ileum'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_14, tissue == 'Ileum'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_14, tissue == 'Ileum'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_14, tissue == 'Ileum'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Ileum, Day 14', subtitle = 'stress: 0.169')



############### END NMDS ###################


##### tissues pairwise adonis #####
rownames(mastermeta) ==  rownames(share.rare2)

rownames(share.rare.tissue) == rownames(mastermeta.tissue)

#share.rare.tissue[1:10,1:10]

PWadon2 <- pairwise.adonis(share.rare.tissue, mastermeta.tissue$set)

colonVScolon <- grep(".*_Colon_.* vs .*_Colon_.*", PWadon2$pairs)
ileumVSileum <- grep(".*_Ileum_.* vs .*_Ileum_.*", PWadon2$pairs)

goodones <- c(colonVScolon, ileumVSileum)
PWadon.tissue <- PWadon2[goodones,]

D4t <- grep("D4_.*_.* vs D4_.*_.*", PWadon.tissue$pairs)
D7t <- grep("D7_.*_.* vs D7_.*_.*", PWadon.tissue$pairs)
D14t <- grep("D14_.*_.* vs D14_.*_.*", PWadon.tissue$pairs)

PWadon.tissue <- PWadon.tissue[c(D4t, D7t, D14t),]


################  Just plotting adonis Fval to control vs day  ##################

PWadon.tocont <- PWadon.tissue[grep('control vs ', PWadon.tissue$pairs),]
PWadon.tocont$day <- gsub('D([0-9]+)_.*', '\\1', PWadon.tocont$pairs)
PWadon.tocont$treatment <- gsub('D[0-9]+_[A-Za-z]+_control vs D[0-9]+_[A-Za-z]+_([A-Za-z]+)', '\\1', PWadon.tocont$pairs)
PWadon.tocont$group <- gsub('D[0-9]+_[A-Za-z]+_control vs (D[0-9]+_[A-Za-z]+_[A-Za-z]+)', '\\1', PWadon.tocont$pairs)
#PWadon.tocont <- PWadon.tocont[order(as.numeric(PWadon.tocont$day)),]
PWadon.tocont$day <- factor(PWadon.tocont$day, levels = c('4', '7', '14'))

PWadon.tocont$group


PWadon.tocont.colon <- PWadon.tocont[grep('Colon', PWadon.tocont$group),]
PWadon.tocont.ileum <- PWadon.tocont[grep('Ileum', PWadon.tocont$group),]

p3 <- ggplot(PWadon.tocont.colon, aes(x=day, y=F.Model, color=treatment)) +
  geom_path(aes(group=treatment), size=1.3) + 
  scale_fill_manual(values = c('#00BA38', '#619CFF')) + 
  geom_label(aes(label=p.value, fill=treatment), color='black', show.legend = FALSE) +
  scale_color_manual(values = c('#00BA38', '#619CFF')) + ylim(1,2.5)
p3 + ggtitle('PERMANOVA comparisions to Control group', subtitle = 'Colon') + 
  ylab('PERMANOVA F vs control\n(difference relative to control)')

p4 <- ggplot(PWadon.tocont.ileum, aes(x=day, y=F.Model, color=treatment)) +
  geom_path(aes(group=treatment), size=1.3) + 
  scale_fill_manual(values = c('#00BA38', '#619CFF')) + 
  geom_label(aes(label=p.value, fill=treatment), color='black', show.legend = FALSE) +
  scale_color_manual(values = c('#00BA38', '#619CFF')) #+ ylim(0,2.5)
p4 + ggtitle('PERMANOVA comparisions to Control group', subtitle = 'Ileum') + 
  ylab('PERMANOVA F vs control\n(difference relative to control)')


################### Deseq2 #####################



library(DESeq2)
library(phyloseq)
otu <- import_mothur(mothur_shared_file = 'FS1.V4.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'FS1.V4.taxonomy')

#meta <- read.table(file = 'V4.metadata.txt', sep = '\t', header = TRUE)

rownames(meta) <- meta$Group
meta$set <- paste(meta$day, meta$tissue, meta$treatment, sep = '_')
phy_meta <- sample_data(meta) 
FS1 <- phyloseq(otu, taxo)
FS1 <- merge_phyloseq(FS1, phy_meta)                       # combines the metadata with this phyloseq object
colnames(tax_table(FS1)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
FS1 <- subset_samples(FS1, experiment == 'X1')


FS1 <- prune_samples(sample_sums(FS1) > 2000, FS1)  # This removes samples that have fewer than 700 sequences associated with them.
FS1 <- prune_taxa(taxa_sums(FS1) > 10, FS1)        # removes OTUs that occur less than 10 times globally

FS1 <- subset_taxa(FS1, !grepl('Mitochondria_ge', Genus))

# If you want to group OTUs uncomment the tax_glom() line and select your desired taxrank
# right now all these analysis are done at the OTU level.

FS1.genus <- FS1
#
FS1.genus <- tax_glom(FS1, taxrank = "Family")
#
##########  IN-FEED VS CONTROL ###########

# FS1.De <- phyloseq_to_deseq2(FS1, ~ set)
# FS1@sam_data
# FS1.De <- DESeq(FS1.De, test = "Wald", fitType = "parametric")

# D0 #

FS1.D0 <- subset_samples(FS1.genus, day == 0)

sample_sums(FS1.D0)
FS1.D0 <- prune_taxa(taxa_sums(FS1.D0) > 1, FS1.D0)
#rowSums(FS1.D0@otu_table)
FS1.D0.De <- phyloseq_to_deseq2(FS1.D0, ~ set)
FS1.D0.De <- DESeq(FS1.D0.De, test = "Wald", fitType = "parametric")


#Deseq.quickplot(DeSeq.object = FS1.D0.De, phyloseq.object = FS1.D0, contrast.vector = c("set","0_Feces_inject","0_Feces_control"))

#
# Day 0 Inject vs Control

FS1.D0.De.results <- results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control"),
                             cooksCutoff = FALSE, pAdjustMethod = 'BH') 
#

sigtab.D0.feces.IvC = FS1.D0.De.results[which(FS1.D0.De.results$padj < .05), ]

length(sigtab.D0.feces.IvC$baseMean)

####### need ggplot code for this #####
#sigtab.D0.feces.IvC = cbind(as(sigtab.D0.feces.IvC, "data.frame"), as(tax_table(FS1.D0)[rownames(sigtab.D0.feces.IvC), ], "matrix"))
#sigtab.D0.feces.IvC$newp <- format(round(sigtab.D0.feces.IvC$padj, digits = 3), scientific = TRUE)
#sigtab.D0.feces.IvC$Treatment <- ifelse(sigtab.D0.feces.IvC$log2FoldChange >=0, "Inject", "Control")

#deseq.D0 <- ggplot(sigtab.D0.feces.IvC, aes(x=reorder(rownames(sigtab.D0.feces.IvC), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
 # geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D0.feces.IvC), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  #theme(axis.text.x=element_text(color = 'black', size = 12),
   #     axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
    #    axis.title.x=element_text(size = 10),
     #   axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 0')+ coord_flip() + scale_fill_manual(values = c('#F8466D', '#619CFF'))
#deseq.D0



# 0 sig diff OTUs at D0 between inj and control
#
# Day 0 Oral vs Control
FS1.D0.De.results <- results(FS1.D0.De, contrast=c("set","0_Feces_oral","0_Feces_control"),
                             cooksCutoff = FALSE, pAdjustMethod = 'BH') 
#FS1.D0.De.results$padj
#
sigtab.D0.feces.OvC = FS1.D0.De.results[which(FS1.D0.De.results$padj < .05), ]
sigtab.D0.feces.OvC

# One OTU sigdiff at D0
# #### All at once? Except D0 and tissues ####
# 
# FS1.fec <- subset_samples(FS1.genus, day != 0 & tissue == 'Feces')
# 
# sample_sums(FS1.fec)
# FS1.fec <- prune_taxa(taxa_sums(FS1.fec) > 1, FS1.fec)
# 
# FS1.fec.De <- phyloseq_to_deseq2(FS1.fec, ~ set)
# 
# FS1.fec.De <- DESeq(FS1.fec.De, test = "Wald", fitType = "parametric")
# 
# Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '4_Feces_oral', '4_Feces_control'), taxlabel = 'Family')
# Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '7_Feces_oral', '7_Feces_control'))
# Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '11_Feces_oral', '11_Feces_control'))
# Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '14_Feces_oral', '14_Feces_control'))
# 
# Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '4_Feces_inject', '4_Feces_control'))
# Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '7_Feces_inject', '7_Feces_control'))
# Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '11_Feces_inject', '11_Feces_control'))
# Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '14_Feces_inject', '14_Feces_control'))
# 
# oral711 <- Deseq.quickplot(FS1.fec.De, FS1.fec, contrast.vector = c('set', '7_Feces_oral', '11_Feces_oral'), taxlabel = 'Class')
# oral711[[2]]
# oral711[[1]]
########## DAY 4 FECES ###########

FS1.D4 <- subset_samples(FS1.genus, day == 4 & tissue == 'Feces')

#sample_sums(FS1.D4)
FS1.D4 <- prune_taxa(taxa_sums(FS1.D4) > 1, FS1.D4)

#rowSums(FS1.D4@otu_table)

FS1.D4.De <- phyloseq_to_deseq2(FS1.D4, ~ set)

FS1.D4.De <- DESeq(FS1.D4.De, test = "Wald", fitType = "parametric")

######### ORAL VS CONTROL #########


res.D4 = results(FS1.D4.De, contrast=c("set","4_Feces_oral","4_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')

sigtab.D4 = res.D4[which(res.D4$padj < .05), ]
sigtab.D4 = cbind(as(sigtab.D4, "data.frame"), as(tax_table(FS1.D4)[rownames(sigtab.D4), ], "matrix"))
#format(sigtab.D4$padj, scientific = TRUE)
sigtab.D4$newp <- format(round(sigtab.D4$padj, digits = 3), scientific = TRUE)
sigtab.D4$Treatment <- ifelse(sigtab.D4$log2FoldChange >=0, "In-feed", "Control")

deseq.D4 <- ggplot(sigtab.D4, aes(x=reorder(rownames(sigtab.D4), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4), y=0, label = paste(Phylum,Class,Family, sep = ' ')), size=3)+ labs(x="Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant Families: feces, Day 4')+ coord_flip() + 
  scale_fill_manual(values = c('#F8766D', '#619CFF'))
deseq.D4



sigtab.D4$OTU <- rownames(sigtab.D4)
sigtab.D4$comp <- 'D4_feces_OvC'

final.sigtab <- sigtab.D4

##### INJ vs Control  #####

####### NEED TO CHANGE COLORS HERE #########

res.D4 = results(FS1.D4.De, contrast=c("set","4_Feces_inject","4_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D4 = res.D4[which(res.D4$padj < .05), ]
sigtab.D4 = cbind(as(sigtab.D4, "data.frame"), as(tax_table(FS1.D4)[rownames(sigtab.D4), ], "matrix"))
format(sigtab.D4$padj, scientific = TRUE)
sigtab.D4$newp <- format(round(sigtab.D4$padj, digits = 3), scientific = TRUE)
sigtab.D4$Treatment <- ifelse(sigtab.D4$log2FoldChange >=0, "inject", "Control")

# deseq.D4 <- ggplot(sigtab.D4, aes(x=reorder(rownames(sigtab.D4), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 4')+ coord_flip() +
#   scale_fill_manual(values = c('#F8766D', '#00BA38'))
# deseq.D4

### Only one OTU dif ab btw the control and injected, OTU 68 Unclassified proteobacteria

sigtab.D4$OTU <- rownames(sigtab.D4)
sigtab.D4$comp <- 'D4_feces_IvC'

final.sigtab <- rbind(final.sigtab, sigtab.D4)


######## Day 7 Feces ###########

FS1.D7 <- subset_samples(FS1.genus, day == 7 & tissue == 'Feces')

sample_sums(FS1.D7)
FS1.D7 <- prune_taxa(taxa_sums(FS1.D7) > 1, FS1.D7)

rowSums(FS1.D7@otu_table)

FS1.D7.De <- phyloseq_to_deseq2(FS1.D7, ~ set)

FS1.D7.De <- DESeq(FS1.D7.De, test = "Wald", fitType = "parametric")

res.D7 = results(FS1.D7.De, contrast=c("set","7_Feces_oral","7_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D7 = res.D7[which(res.D7$padj < .05), ]
sigtab.D7 = cbind(as(sigtab.D7, "data.frame"), as(tax_table(FS1.D7)[rownames(sigtab.D7), ], "matrix"))
format(sigtab.D7$padj, scientific = TRUE)
sigtab.D7$newp <- format(round(sigtab.D7$padj, digits = 3), scientific = TRUE)
sigtab.D7$Treatment <- ifelse(sigtab.D7$log2FoldChange >=0, "In-feed", "Control")

deseq.D7 <- ggplot(sigtab.D7, aes(x=reorder(rownames(sigtab.D7), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7), y=0, label = paste(Class,Family, sep = ' ')), size=3)+ labs(x="Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant Families: feces, Day 7')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
deseq.D7



sigtab.D7$OTU <- rownames(sigtab.D7)
sigtab.D7$comp <- 'D7_feces_OvC'

final.sigtab <- rbind(final.sigtab, sigtab.D7)

######### INSERT Day 7 IvC here  ############

res.D7 = results(FS1.D7.De, contrast=c("set","7_Feces_inject","7_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D7 = res.D7[which(res.D7$padj < .05), ]
sigtab.D7 = cbind(as(sigtab.D7, "data.frame"), as(tax_table(FS1.D7)[rownames(sigtab.D7), ], "matrix"))
#format(sigtab.D7$padj, scientific = TRUE)
sigtab.D7$newp <- format(round(sigtab.D7$padj, digits = 3), scientific = TRUE)
sigtab.D7$Treatment <- ifelse(sigtab.D7$log2FoldChange >=0, "inject", "Control")

deseq.D7 <- ggplot(sigtab.D7, aes(x=reorder(rownames(sigtab.D7), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 7')+ coord_flip() +
  scale_fill_manual(values = c('#F8766D', '#00BA38'))
deseq.D7


sigtab.D7$OTU <- rownames(sigtab.D7)
sigtab.D7$comp <- 'D7_feces_IvC'

final.sigtab <- rbind(final.sigtab, sigtab.D7)

######### D11 ###########
# subset_samples(FS1.genus, day == 11 & tissue == 'Feces') %>% sample_data()

FS1.D11 <- subset_samples(FS1.genus, day == 11 & tissue == 'Feces')

sample_sums(FS1.D11)
FS1.D11 <- prune_taxa(taxa_sums(FS1.D11) > 1, FS1.D11)

rowSums(FS1.D11@otu_table)

FS1.D11.De <- phyloseq_to_deseq2(FS1.D11, ~ set)

FS1.D11.De <- DESeq(FS1.D11.De, test = "Wald", fitType = "parametric")

res.D11 = results(FS1.D11.De, contrast=c("set","11_Feces_oral","11_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D11 = res.D11[which(res.D11$padj < .05), ]
sigtab.D11 = cbind(as(sigtab.D11, "data.frame"), as(tax_table(FS1.D11)[rownames(sigtab.D11), ], "matrix"))
format(sigtab.D11$padj, scientific = TRUE)
sigtab.D11$newp <- format(round(sigtab.D11$padj, digits = 3), scientific = TRUE)
sigtab.D11$Treatment <- ifelse(sigtab.D11$log2FoldChange >=0, "In-feed", "Control")

# deseq.D11 <- ggplot(sigtab.D11, aes(x=reorder(rownames(sigtab.D11), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D11), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 7')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
# deseq.D11
# 


#sigtab.D11$OTU <- rownames(sigtab.D11)
#sigtab.D11$comp <- 'D11_feces_OvC'

# None at D11

#final.sigtab <- rbind(final.sigtab, sigtab.D11)

######### Day 11 IvC here  ############

res.D11 = results(FS1.D11.De, contrast=c("set","11_Feces_inject","11_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D11 = res.D11[which(res.D11$padj < .05), ]
sigtab.D11 = cbind(as(sigtab.D11, "data.frame"), as(tax_table(FS1.D11)[rownames(sigtab.D11), ], "matrix"))
#format(sigtab.D11$padj, scientific = TRUE)
sigtab.D11$newp <- format(round(sigtab.D11$padj, digits = 3), scientific = TRUE)
sigtab.D11$Treatment <- ifelse(sigtab.D11$log2FoldChange >=0, "inject", "Control")

# deseq.D11 <- ggplot(sigtab.D11, aes(x=reorder(rownames(sigtab.D11), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D11), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 7')+ coord_flip() +
#   scale_fill_manual(values = c('#F8766D', '#00BA38'))
# deseq.D11
# 

sigtab.D11$OTU <- rownames(sigtab.D11)
sigtab.D11$comp <- 'D11_feces_IvC'

final.sigtab <- rbind(final.sigtab, sigtab.D11)
### No sig diff at D11...

########### D14 ################

FS1.D14 <- subset_samples(FS1.genus, day == 14 & tissue == 'Feces')

sample_sums(FS1.D14)
FS1.D14 <- prune_taxa(taxa_sums(FS1.D14) > 1, FS1.D14)

rowSums(FS1.D14@otu_table)

FS1.D14.De <- phyloseq_to_deseq2(FS1.D14, ~ set)

FS1.D14.De <- DESeq(FS1.D14.De, test = "Wald", fitType = "parametric")

res.D14 = results(FS1.D14.De, contrast=c("set","14_Feces_oral","14_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')

sigtab.D14 = res.D14[which(res.D14$padj < .05), ]
sigtab.D14 = cbind(as(sigtab.D14, "data.frame"), as(tax_table(FS1.D14)[rownames(sigtab.D14), ], "matrix"))
format(sigtab.D14$padj, scientific = TRUE)
sigtab.D14$newp <- format(round(sigtab.D14$padj, digits = 3), scientific = TRUE)
sigtab.D14$Treatment <- ifelse(sigtab.D14$log2FoldChange >=0, "In-feed", "Control")

deseq.D14 <- ggplot(sigtab.D14, aes(x=reorder(rownames(sigtab.D14), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 14')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
deseq.D14

sigtab.D14$OTU <- rownames(sigtab.D14)
sigtab.D14$comp <- 'D14_feces_OvC'

final.sigtab <- rbind(final.sigtab, sigtab.D14)

##### D14 Inject VS control  #####

res.D14 = results(FS1.D14.De, contrast=c("set","14_Feces_inject","14_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D14 = res.D14[which(res.D14$padj < .05), ]

sigtab.D14
# no sig diff OTUs at D14 between Inject and Control




############ TISSUES  #############

########### D14 ################

# D14 Colon #

FS1.D14.Colon <- subset_samples(FS1.genus, day == 14 & tissue == 'Colon')

sample_sums(FS1.D14.Colon)
FS1.D14.Colon <- prune_taxa(taxa_sums(FS1.D14.Colon) > 1, FS1.D14.Colon)

#rowSums(FS1.D14.Colon@otu_table)

FS1.D14.Colon.De <- phyloseq_to_deseq2(FS1.D14.Colon, ~ set)

FS1.D14.Colon.De <- DESeq(FS1.D14.Colon.De, test = "Wald", fitType = "parametric")

res.D14.Colon = results(FS1.D14.Colon.De, contrast=c("set","14_Colon_oral","14_Colon_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D14.Colon = res.D14.Colon[which(res.D14.Colon$padj < .05), ]
sigtab.D14.Colon = cbind(as(sigtab.D14.Colon, "data.frame"), as(tax_table(FS1.D14.Colon)[rownames(sigtab.D14.Colon), ], "matrix"))
format(sigtab.D14.Colon$padj, scientific = TRUE)
sigtab.D14.Colon$newp <- format(round(sigtab.D14.Colon$padj, digits = 3), scientific = TRUE)
sigtab.D14.Colon$Treatment <- ifelse(sigtab.D14.Colon$log2FoldChange >=0, "In-feed", "Control")

# deseq.D14.Colon <- ggplot(sigtab.D14.Colon, aes(x=reorder(rownames(sigtab.D14.Colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.Colon), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 14')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
# deseq.D14.Colon
# 


#sigtab.D14.Colon$OTU <- rownames(sigtab.D14.Colon)
#sigtab.D14.Colon$comp <- 'D14_colon_OvC'

final.sigtab <- rbind(final.sigtab, sigtab.D14.Colon)

# only 1 OTU diff btw groups #
# none diff btw groups

######## D14 colon INJ vs Control ########

res.D14.Colon = results(FS1.D14.Colon.De, contrast=c("set","14_Colon_inject","14_Colon_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D14.Colon = res.D14.Colon[which(res.D14.Colon$padj < .05), ]
sigtab.D14.Colon = cbind(as(sigtab.D14.Colon, "data.frame"), as(tax_table(FS1.D14.Colon)[rownames(sigtab.D14.Colon), ], "matrix"))
format(sigtab.D14.Colon$padj, scientific = TRUE)
sigtab.D14.Colon$newp <- format(round(sigtab.D14.Colon$padj, digits = 3), scientific = TRUE)
sigtab.D14.Colon$Treatment <- ifelse(sigtab.D14.Colon$log2FoldChange >=0, "In-feed", "Control")

# deseq.D14.Colon <- ggplot(sigtab.D14.Colon, aes(x=reorder(rownames(sigtab.D14.Colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.Colon), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 14')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
# deseq.D14.Colon

# No sig diff OTUs at D14 in the colonic mucosa between inject and control
# none in V4 either


################
# ONLY 9 SAMPLES IN D7 Colon #
# THIS IS GARBAGE #

# D7 Colon #

FS1.D7.Colon <- subset_samples(FS1.genus, day == 7 & tissue == 'Colon')
sample_data(FS1.D7.Colon)

# 3 controls, 3 inject, 3 oral

#sample_sums(FS1.D7.Colon)

FS1.D7.Colon <- prune_taxa(taxa_sums(FS1.D7.Colon) > 1, FS1.D7.Colon)

#rowSums(FS1.D7.Colon@otu_table)

FS1.D7.Colon.De <- phyloseq_to_deseq2(FS1.D7.Colon, ~ set)

FS1.D7.Colon.De <- DESeq(FS1.D7.Colon.De, test = "Wald", fitType = "parametric")

res.D7.Colon = results(FS1.D7.Colon.De, contrast=c("set","7_Colon_oral","7_Colon_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D7.Colon = res.D7.Colon[which(res.D7.Colon$padj < .05), ]
sigtab.D7.Colon = cbind(as(sigtab.D7.Colon, "data.frame"), as(tax_table(FS1.D7.Colon)[rownames(sigtab.D7.Colon), ], "matrix"))
format(sigtab.D7.Colon$padj, scientific = TRUE)
sigtab.D7.Colon$newp <- format(round(sigtab.D7.Colon$padj, digits = 3), scientific = TRUE)
sigtab.D7.Colon$Treatment <- ifelse(sigtab.D7.Colon$log2FoldChange >=0, "In-feed", "Control")

deseq.D7.Colon <- ggplot(sigtab.D7.Colon, aes(x=reorder(rownames(sigtab.D7.Colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.Colon), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: Colonic Mucosa, Day 7')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
deseq.D7.Colon



####### D7 colon Inj vs control   ########
############# HERE #############
############## HERE ############

res.D7.Colon = results(FS1.D7.Colon.De, contrast=c("set","7_Colon_inject","7_Colon_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D7.Colon = res.D7.Colon[which(res.D7.Colon$padj < .05), ]
sigtab.D7.Colon = cbind(as(sigtab.D7.Colon, "data.frame"), as(tax_table(FS1.D7.Colon)[rownames(sigtab.D7.Colon), ], "matrix"))
format(sigtab.D7.Colon$padj, scientific = TRUE)
sigtab.D7.Colon$newp <- format(round(sigtab.D7.Colon$padj, digits = 3), scientific = TRUE)
sigtab.D7.Colon$Treatment <- ifelse(sigtab.D7.Colon$log2FoldChange >=0, "Inject", "Control")
# 
# deseq.D7.Colon <- ggplot(sigtab.D7.Colon, aes(x=reorder(rownames(sigtab.D7.Colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.Colon), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: feces, Day 14')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#00BA38'))
# deseq.D7.Colon

# none
# none in V4 either
# D4 Colon #

FS1.D4.Colon <- subset_samples(FS1.genus, day == 4 & tissue == 'Colon')
sample_data(FS1.D4.Colon)
# 4 inject, 4 oral, 5 control

sample_sums(FS1.D4.Colon)
FS1.D4.Colon <- prune_taxa(taxa_sums(FS1.D4.Colon) > 1, FS1.D4.Colon)

rowSums(FS1.D4.Colon@otu_table)

FS1.D4.Colon.De <- phyloseq_to_deseq2(FS1.D4.Colon, ~ set)

FS1.D4.Colon.De <- DESeq(FS1.D4.Colon.De, test = "Wald", fitType = "parametric")

res.D4.Colon = results(FS1.D4.Colon.De, contrast=c("set","4_Colon_oral","4_Colon_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D4.Colon = res.D4.Colon[which(res.D4.Colon$padj < .05), ]
sigtab.D4.Colon = cbind(as(sigtab.D4.Colon, "data.frame"), as(tax_table(FS1.D4.Colon)[rownames(sigtab.D4.Colon), ], "matrix"))
format(sigtab.D4.Colon$padj, scientific = TRUE)
sigtab.D4.Colon$newp <- format(round(sigtab.D4.Colon$padj, digits = 3), scientific = TRUE)
sigtab.D4.Colon$Treatment <- ifelse(sigtab.D4.Colon$log2FoldChange >=0, "In-feed", "Control")

deseq.D4.Colon <- ggplot(sigtab.D4.Colon, aes(x=reorder(rownames(sigtab.D4.Colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.Colon), y=0, label = paste(Class,Family, sep = ' ')), size=3)+ labs(x="Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant Families: Colon, Day 4')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
deseq.D4.Colon



sigtab.D4.Colon$OTU <- rownames(sigtab.D4.Colon)
sigtab.D4.Colon$comp <- 'D4_colon_OvC'

final.sigtab <- rbind(final.sigtab, sigtab.D4.Colon)

######### D4 colon Inj vs control ##########

res.D4.Colon = results(FS1.D4.Colon.De, contrast=c("set","4_Colon_inject","4_Colon_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D4.Colon = res.D4.Colon[which(res.D4.Colon$padj < .05), ]
sigtab.D4.Colon = cbind(as(sigtab.D4.Colon, "data.frame"), as(tax_table(FS1.D4.Colon)[rownames(sigtab.D4.Colon), ], "matrix"))
format(sigtab.D4.Colon$padj, scientific = TRUE)
sigtab.D4.Colon$newp <- format(round(sigtab.D4.Colon$padj, digits = 3), scientific = TRUE)
sigtab.D4.Colon$Treatment <- ifelse(sigtab.D4.Colon$log2FoldChange >=0, "Inject", "Control")

# deseq.D4.Colon <- ggplot(sigtab.D4.Colon, aes(x=reorder(rownames(sigtab.D4.Colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.Colon), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: Colon, Day 4', subtitle = 'Control VS Inject')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#00BA38')) + 
#   ylim(-8,8)
# deseq.D4.Colon



sigtab.D4.Colon$OTU <- rownames(sigtab.D4.Colon)
sigtab.D4.Colon$comp <- 'D4_colon_IvC'

final.sigtab <- rbind(final.sigtab, sigtab.D4.Colon)


############ D4 colon INJ vs ORAL ##############

# res.D4.Colon = results(FS1.D4.Colon.De, contrast=c("set","4_Colon_inject","4_Colon_oral"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
# #results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
# sigtab.D4.Colon = res.D4.Colon[which(res.D4.Colon$padj < .05), ]
# sigtab.D4.Colon = cbind(as(sigtab.D4.Colon, "data.frame"), as(tax_table(FS1.D4.Colon)[rownames(sigtab.D4.Colon), ], "matrix"))
# format(sigtab.D4.Colon$padj, scientific = TRUE)
# sigtab.D4.Colon$newp <- format(round(sigtab.D4.Colon$padj, digits = 3), scientific = TRUE)
# sigtab.D4.Colon$Treatment <- ifelse(sigtab.D4.Colon$log2FoldChange >=0, "Inject", "In-feed")
# 
# deseq.D4.Colon <- ggplot(sigtab.D4.Colon, aes(x=reorder(rownames(sigtab.D4.Colon), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.Colon), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: Colon, Day 4')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
# deseq.D4.Colon




################# Ileum Oral vs control ################
# D14 ileum #
### This seems very suspicious. There are very few differences between the treatment groups at day 14. Yet, here there are like a hundred or so OTUs that are different....
### I think this should be omitted.  I don't see a lot of value in any of the Ileum data actually.  If part of it is suspect then all of it is.


FS1.D14.ileum <- subset_samples(FS1.genus, day == 14 & tissue == 'Ileum')
sample_data(FS1.D14.ileum)
# 5 control, 8 inject, 7 in-feed

sample_sums(FS1.D14.ileum)
FS1.D14.ileum <- prune_taxa(taxa_sums(FS1.D14.ileum) > 1, FS1.D14.ileum)

rowSums(FS1.D14.ileum@otu_table)

FS1.D14.ileum.De <- phyloseq_to_deseq2(FS1.D14.ileum, ~ set)

FS1.D14.ileum.De <- DESeq(FS1.D14.ileum.De, test = "Wald", fitType = "parametric")

res.D14.ileum = results(FS1.D14.ileum.De, contrast=c("set","14_Ileum_oral","14_Ileum_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D14.ileum = res.D14.ileum[which(res.D14.ileum$padj < .05), ]
sigtab.D14.ileum = cbind(as(sigtab.D14.ileum, "data.frame"), as(tax_table(FS1.D14.ileum)[rownames(sigtab.D14.ileum), ], "matrix"))
format(sigtab.D14.ileum$padj, scientific = TRUE)
sigtab.D14.ileum$newp <- format(round(sigtab.D14.ileum$padj, digits = 3), scientific = TRUE)
sigtab.D14.ileum$Treatment <- ifelse(sigtab.D14.ileum$log2FoldChange >=0, "In-feed", "Control")
# also suspicious that all the pvalues are the same 
# deseq.D14.ileum <- ggplot(sigtab.D14.ileum, aes(x=reorder(rownames(sigtab.D14.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.ileum), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: ileum, Day 14')+ coord_flip() +
#   scale_fill_manual(values = c('#F8766D', '#619CFF')) 
# deseq.D14.ileum


sigtab.D14.ileum$OTU <- rownames(sigtab.D14.ileum)
sigtab.D14.ileum$comp <- 'D14_Ileum_OvC'

final.sigtab <- rbind(final.sigtab, sigtab.D14.ileum)


############## D14 Ileum Inj vs control  ################

res.D14.ileum = results(FS1.D14.ileum.De, contrast=c("set","14_Ileum_inject","14_Ileum_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D14.ileum = res.D14.ileum[which(res.D14.ileum$padj < .05), ]
sigtab.D14.ileum = cbind(as(sigtab.D14.ileum, "data.frame"), as(tax_table(FS1.D14.ileum)[rownames(sigtab.D14.ileum), ], "matrix"))
format(sigtab.D14.ileum$padj, scientific = TRUE)
sigtab.D14.ileum$newp <- format(round(sigtab.D14.ileum$padj, digits = 3), scientific = TRUE)
sigtab.D14.ileum$Treatment <- ifelse(sigtab.D14.ileum$log2FoldChange >=0, "Inject", "Control")

# deseq.D14.ileum <- ggplot(sigtab.D14.ileum, aes(x=reorder(rownames(sigtab.D14.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D14.ileum), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: ileum, Day 14')+ coord_flip() +
#   scale_fill_manual(values = c('#F8766D', '#00BA38')) 
# deseq.D14.ileum

sigtab.D14.ileum$OTU <- rownames(sigtab.D14.ileum)
sigtab.D14.ileum$comp <- 'D14_Ileum_IvC'

final.sigtab <- rbind(final.sigtab, sigtab.D14.ileum)



################
# ONLY 8 SAMPLES IN D7 ileum #
# THIS IS GARBAGE #

# D7 ileum #

FS1.D7.ileum <- subset_samples(FS1.genus, day == 7 & tissue == 'Ileum')
sample_data(FS1.D7.ileum)
# 3 control, 3 inject, 2 in-feed

sample_sums(FS1.D7.ileum)

FS1.D7.ileum <- prune_taxa(taxa_sums(FS1.D7.ileum) > 1, FS1.D7.ileum)

rowSums(FS1.D7.ileum@otu_table)

FS1.D7.ileum.De <- phyloseq_to_deseq2(FS1.D7.ileum, ~ set)

FS1.D7.ileum.De <- DESeq(FS1.D7.ileum.De, test = "Wald", fitType = "parametric")

res.D7.ileum = results(FS1.D7.ileum.De, contrast=c("set","7_Ileum_oral","7_Ileum_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D7.ileum = res.D7.ileum[which(res.D7.ileum$padj < .05), ]
sigtab.D7.ileum = cbind(as(sigtab.D7.ileum, "data.frame"), as(tax_table(FS1.D7.ileum)[rownames(sigtab.D7.ileum), ], "matrix"))
format(sigtab.D7.ileum$padj, scientific = TRUE)
sigtab.D7.ileum$newp <- format(round(sigtab.D7.ileum$padj, digits = 3), scientific = TRUE)
sigtab.D7.ileum$Treatment <- ifelse(sigtab.D7.ileum$log2FoldChange >=0, "In-feed", "Control")

deseq.D7.ileum <- ggplot(sigtab.D7.ileum, aes(x=reorder(rownames(sigtab.D7.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
  geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D7.ileum), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
  theme(axis.text.x=element_text(color = 'black', size = 12),
        axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
        axis.title.x=element_text(size = 10),
        axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: Ileum, Day 14')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
deseq.D7.ileum



# ######## D7 ileum inj vs control  ###########
# res.D7.ileum = results(FS1.D7.ileum.De, contrast=c("set","7_Ileum_inject","7_Ileum_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
# #results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
# sigtab.D7.ileum = res.D7.ileum[which(res.D7.ileum$padj < .05), ]
# sigtab.D7.ileum = cbind(as(sigtab.D7.ileum, "data.frame"), as(tax_table(FS1.D7.ileum)[rownames(sigtab.D7.ileum), ], "matrix"))
# format(sigtab.D7.ileum$padj, scientific = TRUE)
# sigtab.D7.ileum$newp <- format(round(sigtab.D7.ileum$padj, digits = 3), scientific = TRUE)
# sigtab.D7.ileum$Treatment <- ifelse(sigtab.D7.ileum$log2FoldChange >=0, "In-feed", "Control")
# 
# FS1.D4.ileum <- subset_samples(FS1.genus, day == 4 & tissue == 'Ileum')
# sample_data(FS1.D4.ileum)
# 
# 

# sample_sums(FS1.D4.ileum)

FS1.D4.ileum <- prune_taxa(taxa_sums(FS1.D4.ileum) > 1, FS1.D4.ileum)

# rowSums(FS1.D4.ileum@otu_table)

FS1.D4.ileum.De <- phyloseq_to_deseq2(FS1.D4.ileum, ~ set)

FS1.D4.ileum.De <- DESeq(FS1.D4.ileum.De, test = "Wald", fitType = "parametric")

res.D4.ileum = results(FS1.D4.ileum.De, contrast=c("set","4_Ileum_oral","4_Ileum_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
#results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
sigtab.D4.ileum = res.D4.ileum[which(res.D4.ileum$padj < .05), ]
sigtab.D4.ileum = cbind(as(sigtab.D4.ileum, "data.frame"), as(tax_table(FS1.D4.ileum)[rownames(sigtab.D4.ileum), ], "matrix"))
format(sigtab.D4.ileum$padj, scientific = TRUE)
sigtab.D4.ileum$newp <- format(round(sigtab.D4.ileum$padj, digits = 3), scientific = TRUE)
sigtab.D4.ileum$Treatment <- ifelse(sigtab.D4.ileum$log2FoldChange >=0, "In-feed", "Control")

# ############ D4 ileum Inject vs Control #######
# res.D4.ileum = results(FS1.D4.ileum.De, contrast=c("set","4_Ileum_inject","4_Ileum_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
# #results(FS1.D0.De, contrast=c("set","0_Feces_inject","0_Feces_control")) 
# sigtab.D4.ileum = res.D4.ileum[which(res.D4.ileum$padj < .05), ]
# sigtab.D4.ileum = cbind(as(sigtab.D4.ileum, "data.frame"), as(tax_table(FS1.D4.ileum)[rownames(sigtab.D4.ileum), ], "matrix"))
# format(sigtab.D4.ileum$padj, scientific = TRUE)
# sigtab.D4.ileum$newp <- format(round(sigtab.D4.ileum$padj, digits = 3), scientific = TRUE)
# sigtab.D4.ileum$Treatment <- ifelse(sigtab.D4.ileum$log2FoldChange >=0, "In-feed", "Control")
# 
# deseq.D4.ileum <- ggplot(sigtab.D4.ileum, aes(x=reorder(rownames(sigtab.D4.ileum), log2FoldChange), y=log2FoldChange, fill = Treatment)) +
#   geom_bar(stat='identity') + geom_text(aes(x=rownames(sigtab.D4.ileum), y=0, label = paste(Class,Family,Genus, sep = ' ')), size=3)+ labs(x="Genus")+
#   theme(axis.text.x=element_text(color = 'black', size = 12),
#         axis.text.y=element_text(color = 'black', size=12, face = 'italic'), 
#         axis.title.x=element_text(size = 10),
#         axis.title.y=element_text(size = 10))+ ggtitle('Differentially abundant OTUs: ileum, Day 4')+ coord_flip() + scale_fill_manual(values = c('#F8766D', '#619CFF'))
# deseq.D4.ileum
# 
# 
# 

unique(final.sigtab$OTU)

#final.sigtab[final.sigtab$tissue == 'Ileum',]


#final.sigtab[final.sigtab$tissue == 'Ileum',][final.sigtab[final.sigtab$tissue == 'Ileum',]$OTU %in% final.sigtab[final.sigtab$tissue == 'colon',]$OTU,]


### Writing out significantly different OTUs.  Dont really trust the Ileum data, so going to omit that.
### Two OTUs from the Ileum are also seen in the feces and colon data, but very low abundance there.

# going to repeat this with different taxanomic groupings OTU up to phylum
final.sigtab$tissue <- gsub('D[1-9]+_([A-Za-z]+)_[A-Za-z]+','\\1',final.sigtab$comp)

final.sigtab <- final.sigtab[final.sigtab$tissue != 'Ileum',]

#write.csv(final.sigtab, file = 'sig_diff_OTUs.csv')
#write.csv(final.sigtab, file = 'sig_diff_genus.csv')
#write.csv(final.sigtab, file = 'sig_diff_family.csv')
#write.csv(final.sigtab, file = 'sig_diff_order.csv')
#write.csv(final.sigtab, file = 'sig_diff_class.csv')
#write.csv(final.sigtab, file = 'sig_diff_phylum.csv')



read.csv('sig_diff_OTUs.csv')
read.csv('sig_diff_genus.csv')
read.csv('sig_diff_family.csv')
read.csv('sig_diff_order.csv')
read.csv('sig_diff_class.csv')
read.csv('sig_diff_phylum.csv')

#### Add Plots here ###



# #######
# tax <- extract_mothur_tax(filename = 'FS1.V4.taxonomy')
# swap <- otu_tax_labels(tax)
# swap[unique(final.sigtab$OTU)]
# 
# tax


#View(final.sigtab)


############ Stacked BARS  ################

FS1.genus

FS1.phylum <- tax_glom(FS1, 'Phylum')

phyla_tab <- as.data.frame(t(FS1.phylum@otu_table))
colnames(phyla_tab) <- FS1.phylum@tax_table[,2]
phyla_tab <- phyla_tab/rowSums(phyla_tab)



phyla_tab$Group <- rownames(phyla_tab)
fobar <- merge(meta, phyla_tab, by = 'Group')
#write.table(fobar, 'phyla_proportions.txt', quote = FALSE, row.names = FALSE, col.names = TRUE, sep = '\t')

fobar.gather <- fobar %>% gather(phylum, value, -(Group:set))  #converts to long form dataframe, this is handy for using ggplot faceting functions, check out tidyverse tutorials

# ggplot(data=filter(fobar.gather, tissue == 'Feces'), aes(x=day, y=value, fill=phylum)) +
#   geom_bar(stat = 'identity') +
#   facet_grid(treatment~tissue)
# 
# 

phyla_tab <- phyla_tab[,-(length(phyla_tab[1,]))]

phyla_tab <- phyla_tab[,(colSums(phyla_tab)/length(phyla_tab$Bacteroidetes))>0.01] # greater than 1% avg abundance across all samples all tissues

phyla_tab$Group <- rownames(phyla_tab)

fobar2 <- merge(meta, phyla_tab, by = 'Group')

fobar2.gather <- fobar2 %>% gather(phylum, value, -(Group:set))  #converts to long form dataframe, this is handy for using ggplot faceting functions, check out tidyverse tutorials

length(unique(fobar2.gather$phylum))
num_taxa <- length(unique(fobar2.gather$phylum))
fobar2.gather <- fobar2.gather %>% group_by(set) %>% mutate(value2=(value/(length(set)/num_taxa))*100) # This essentially makes it so the bars all add up to 100% even though there are different numbers of animals in each set


#fobar2.gather %>% group_by(set) %>% summarise(mean=mean(value), sd=sd(value), se=sd/n())

#fobar2.gather$value <- (fobar2.gather$value/4)*100

fobar2.gather$day <- as.numeric(fobar2.gather$day)

ggplot(data=filter(fobar2.gather, tissue == 'Feces'), aes(x=day, y=value2, fill=phylum)) +
  geom_bar(stat = 'identity') +
  facet_grid(treatment~tissue)+ scale_fill_brewer(palette = 'Set1') + ylab('Percent of total community') + 
  ggtitle('Phyla over time: Feces')


fo_area <- fobar2.gather %>% mutate(set_Phylum=paste(set, phylum, sep = '_')) %>%
  group_by(set_Phylum) %>% summarise(set_val=mean(value), stddev= sd(value), se=stddev/sqrt(n()))

fo_area$set_Phylum

fo_area$day <- as.numeric(sub('([0-9]+)_([A-Za-z]+)_([A-Za-z]+)_([A-Za-z_]+)', '\\1', fo_area$set_Phylum))
fo_area$tissue <- sub('([0-9]+)_([A-Za-z]+)_([A-Za-z]+)_([A-Za-z_]+)', '\\2', fo_area$set_Phylum)
fo_area$treatment <-  sub('([0-9]+)_([A-Za-z]+)_([A-Za-z]+)_([A-Za-z_]+)', '\\3', fo_area$set_Phylum)
fo_area$phylum <- sub('([0-9]+)_([A-Za-z]+)_([A-Za-z]+)_([A-Za-z_]+)', '\\4', fo_area$set_Phylum)

filter(fo_area, tissue == 'Feces') %>% ggplot(aes(x=day, y=set_sum, color=treatment)) + 
  geom_line(size = 1) + facet_wrap(~phylum, scales = "free") + scale_color_brewer(palette = 'Set1') + 
  ggtitle("Phyla over time: Feces")

filter(fo_area, tissue == 'Feces' & phylum == 'Proteobacteria') %>%
  ggplot(aes(x=day)) + 
  geom_line(aes(y=set_val, color=treatment), size = 1)+ 
  geom_ribbon(aes(ymin = set_val-(se), ymax = set_val+(se), fill = treatment), alpha =.4) + 
  # scale_color_brewer(palette = 'Set1') + 
  # scale_fill_brewer(palette = 'Set1') + 
  ggtitle("Proteobacteria over time: Feces") + ylab('p')


filter(fo_area, tissue == 'Feces') %>%
  ggplot(aes(x=day)) + 
  geom_line(aes(y=set_val, color=treatment), size = 1)+ 
  geom_ribbon(aes(ymin = set_val-(se), ymax = set_val+(se), fill = treatment), alpha =.4) + 
  # scale_color_brewer(palette = 'Set1') + 
  # scale_fill_brewer(palette = 'Set1') + 
  ggtitle("Proteobacteria over time: Feces") + ylab('p') + facet_wrap(~phylum, scales = 'free')


filter(fo_area, tissue == 'Colon') %>% ggplot(aes(x=day, y=set_sum, color=treatment)) + 
  geom_line(size = 2) + facet_grid(phylum~., scales = "free")+ scale_color_brewer(palette = 'Set1') + 
  ggtitle("Phyla over time: Colon")

#filter(fo_area, tissue == 'Ileum') %>% ggplot(aes(x=day, y=set_sum, fill=phylum)) + geom_area(position = 'stack') + facet_grid(treatment~.)

# 
# ggplot(data=filter(fobar2.gather, tissue == 'Feces' & treatment == 'control')) +
#   geom_area(aes(x=day, y=value2), position = 'stack') +
#   scale_fill_brewer(palette = 'Set2') + ylab('Percent of total community')

# fobar2.gather$day
# 
# 
# ggplot(data=filter(fobar2.gather, tissue == 'Feces'), aes(x=day, y=value2, fill=phylum)) +
#   geom_area() +
#   facet_grid(treatment~.)+ scale_fill_brewer(palette = 'Set2') + ylab('Percent of total community')

# fobar2.gather %>% ggplot(aes(x=day, y=value, group=set, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
#   facet_grid(tissue~phylum, scales = 'free')

# Boxplot comparing Phyla between treatments at various tissues

fobar2.gather$percent_tot <- fobar2.gather$value * 100

fobar2.gather %>% filter(tissue == 'Feces') %>% ggplot(aes(x=day, y=percent_tot, group=set, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~phylum, scales = 'free') + ggtitle('Major Phyla over time: Feces')

fobar2.gather %>% filter(tissue == 'Colon') %>% ggplot(aes(x=day, y=percent_tot, group=set, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~phylum, scales = 'free') + ggtitle('Major Phyla over time: Colon')

fobar2.gather %>% filter(tissue == 'Ileum') %>% ggplot(aes(x=day, y=percent_tot, group=set, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~phylum, scales = 'free')+ ggtitle('Major Phyla over time: Ileum')



fobar2.gather %>% filter(tissue == 'Feces' & day == 7) %>% ggplot(aes(x=treatment, y=percent_tot, group=set, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~phylum, scales = 'free') + ggtitle("Day 7: Feces")

fobar2.gather %>% filter(tissue == 'Feces' & day == 4) %>% ggplot(aes(x=treatment, y=percent_tot, group=set, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~phylum, scales = 'free')+ ggtitle("Day 4: Feces")

fobar2.gather %>% filter(tissue == 'Colon' & day == 7) %>% ggplot(aes(x=treatment, y=percent_tot, group=set, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~phylum, scales = 'free')+ ggtitle("Day 7: Colon")

fobar2.gather %>% filter(tissue == 'Colon' & day == 4) %>% ggplot(aes(x=treatment, y=percent_tot, group=set, fill=treatment)) + geom_boxplot(position = position_dodge2(preserve = 'total')) +
  facet_wrap(~phylum, scales = 'free')+ ggtitle("Day 4: Colon")


################### end stacked bars  ################

# 
# 
# share.rare3 <- as.data.frame(t(FS1.genus@otu_table))
# 
# 
# swap[colnames(share.rare3)]
# 
# share.rare3 <- share.rare3/rowSums(share.rare3)*100
# 
# ################################# END DESEQ  #######################
# # 
# 
# share.rare3 <- as.data.frame(t(FS1.genus@otu_table))
# 
# 
# share.rare3 <- share.rare3/rowSums(share.rare3)*100
# feces4.share <- share.rare3
# 
# ###################################################
# 
# abxcor4.meta <- mastermeta[mastermeta$day == 4 & mastermeta$tissue == 'Feces',]
# 
# feces4.share <- feces4.share[rownames(feces4.share) %in% abxcor4.meta$Group,]
# 
# rowSums(feces4.share)
# 
# feces4.share <- feces4.share[,colSums(feces4.share)> 0.1]
# 
# abxcor4.meta$Group == rownames(feces4.share) 
# 
# feces4.share <- feces4.share[match(abxcor4.meta$Group, rownames(feces4.share)),]
# 
# abxcor4.meta$Group == rownames(feces4.share) 
# 
# 
# #
# abxcor4.meta$tissue <- 'feces'
# abxcor4.meta$group2 <- paste(abxcor4.meta$pig_number, abxcor4.meta$day, abxcor4.meta$tissue, sep = '_')
# 
# colnames(feces.abx)[7] <- 'group2'
# 
# 
# feces.abx$group2
# 
# feces4.abx <- merge(x=feces.abx, y=abxcor4.meta, by='group2')
# 
# 
# feces4.share <- feces4.share[match(feces4.abx$group, rownames(feces4.share)),]
# 
# 
# rownames(feces4.share) == feces4.abx$group
# #
# # rownames(feces4.share)[46]
# # feces.abx$group[46]
# 
# 
# 
# 
# feces4.share <- as.matrix(feces4.share)
# colnames(feces4.share) <- swap[colnames(feces4.share)]
# maybe <- rcorr(y=feces4.abx$concentration, x=feces4.share)
# abxcor <- rcorr_to_ggnet(maybe)
# keepers <- c(grep('^y$', abxcor$from), grep('^y$', abxcor$to))
# abxcor <- abxcor[keepers,]
# 
# 
# 
# 
# 
# 
# nodes <- rbind(gather_nodes(feces4.share, '16S'), 
#                gather_nodes(feces4.abx, 'abx'))
# 
# nodes$node <- as.character(nodes$node)
# nodes$node[770] <- 'y'
# 
# 
# all <- fortify(as.edgedf(abxcor), nodes)
# 
# p.butnet1 <- ggplot(all, aes(from_id = from_id, to_id = to_id, label=from, color=type)) + 
#   geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)), 
#            aes(color = type, label = from_id),
#            linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
#            repel = TRUE, fontsize=3, singletons = FALSE,labelcolour="black",
#            labelgeom = 'text') +
#   theme_net()
# 
# p.butnet1
# 
# ########################################### DAY 7 CORRELATIONS #####################
# 
# abxcor4.meta <- mastermeta[mastermeta$day == 7 & mastermeta$tissue == 'Feces',]
# 
# feces7.share <- share.rare3[rownames(share.rare3) %in% abxcor4.meta$group,]
# 
# rowSums(feces7.share)
# 
# feces7.share <- feces7.share[,colSums(feces7.share)> 0.1]
# 
# abxcor4.meta$group == rownames(feces7.share) 
# 
# feces7.share <- feces7.share[match(abxcor4.meta$group, rownames(feces7.share)),]
# 
# abxcor4.meta$group == rownames(feces7.share) 
# 
# 
# #
# abxcor4.meta$tissue <- 'feces'
# abxcor4.meta$group2 <- paste(abxcor4.meta$pig_number, abxcor4.meta$day, abxcor4.meta$tissue, sep = '_')
# 
# colnames(feces.abx)[7] <- 'group2'
# 
# #prob need to change this
# feces.abx <- merge(abxcor4.meta, feces.abx, by='group2')
# 
# feces7.share <- feces7.share[match(feces.abx$group, rownames(feces7.share)),]
# 
# 
# rownames(feces7.share) == feces.abx$group
# #
# # rownames(feces7.share)[46]
# # feces.abx$group[46]
# 
# 
# 
# 
# feces7.share <- as.matrix(feces7.share)
# colnames(feces7.share) <- swap[colnames(feces7.share)]
# maybe <- rcorr(y=feces.abx$concentration, x=feces7.share)
# abxcor <- rcorr_to_ggnet(maybe)
# keepers <- c(grep('^y$', abxcor$from), grep('^y$', abxcor$to))
# abxcor <- abxcor[keepers,]
# 
# 
# 
# 
# 
# 
# nodes <- rbind(gather_nodes(feces7.share, '16S'), 
#                gather_nodes(feces.abx, 'abx'))
# 
# nodes$node <- as.character(nodes$node)
# # nodes$node[1108] <- 'y'
# 
# 
# all <- fortify(as.edgedf(abxcor), nodes)
# 
# p.butnet1 <- ggplot(all, aes(from_id = from_id, to_id = to_id, label=from, color=type)) + 
#   geom_net(layout.alg = 'fruchtermanreingold', #layout.par = list(c('cell.pointpointrad', 30), c('niter', 1000)), 
#            aes(color = type, label = from_id),
#            linewidth = 0.5, size = 5, vjust = 0, alpha = 0.3,
#            repel = TRUE, fontsize=3, singletons = FALSE,labelcolour="black",
#            labelgeom = 'text') +
#   theme_net()
# 
# p.butnet1
# 
# 
# 
# 
