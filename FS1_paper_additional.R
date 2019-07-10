getwd()
setwd('~/FS1/Helping_Nicole/')

library(tidyverse)
library(vegan)
library(dplyr)
library(ggplot2)
library(Hmisc)
library(geomnet)
library(funfuns) 
library(cowplot)

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

tax<-read.table("FS1.V4.taxonomy", header = TRUE, stringsAsFactors = FALSE, sep="\t")  #added by NR because it seems to be missing
colnames(tax)
library(tidyverse)
tax <- tax %>% separate(Taxonomy, into=c("Kingdom", "Phylum", "Class", "Order", "Genus","Species"), sep=";", extra="merge", drop=TRUE)
#still NR - need to get rid of the numbers following each assignment
tax$Genus <- gsub("\\s*\\([^\\)]+\\)", "", tax$Genus)


#NR - can't figure out how to do this for every column at once so just did Genus for now
#JT - this does the gsub command on every column 
apply(X = tax, MARGIN = 2, FUN = gsub, pattern = "\\s*\\([^\\)]+\\)", replace = "", tax$Genus)
#

tax[grep('Mitochondria', tax$Genus),]
shared <- shared[,-grep('Mitochondria', tax$Genus)]
shared <- shared[,-grep('Chloroplast', tax$Genus)]  #NR added this as well

shared <- shared[rowSums(shared) > 2000,] # removes samples with less than 2000 reads
shared <- shared[,colSums(shared) > 10] # removes OTUs occuring less than 10 times globally
mastermeta <- mastermeta[rownames(mastermeta) %in% rownames(shared),]

# 314 samples

all_groups <- mastermeta %>% group_by(set) %>% summarise(no_rows = length(set)) # the number of samples falling into each set
mastermeta <- mastermeta[mastermeta$set != 'D7_Ileum_inject',] # Only 2 samples in the D7_Ileum_inject group, removed both.
# write.csv(all_groups, 'pigs_per_group.csv')
shared <- shared[rownames(shared) %in% rownames(mastermeta),]

min(rowSums(shared))
library(vegan)
share.rare <- rrarefy(shared, min(rowSums(shared)))  # rarrefied for alpha and beta diversity calcs

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

sow_day <- mastermeta.feces %>% group_by(day, Sow, treatment) %>% tally()

sow_day %>% ggplot(aes(x=Sow, y=n, fill=treatment)) + geom_col() + facet_wrap(~day, nrow = 1)

adonis(share.rare ~ day + treatment, data = mastermeta.feces)

D0_met <- mastermeta.feces %>% filter(day == 0)
D0_share <- share.rare[rownames(share.rare) %in% D0_met$Group,]

D4_met <- mastermeta.feces %>% filter(day == 4)
D4_share <- share.rare[rownames(share.rare) %in% D4_met$Group,]

rownames(D4_share) == D4_met$Group


D7_met <- mastermeta.feces %>% filter(day == 7)
D7_share <- share.rare[rownames(share.rare) %in% D7_met$Group,]

rownames(D7_share) == D7_met$Group


D11_met <- mastermeta.feces %>% filter(day == 11)
D11_share <- share.rare[rownames(share.rare) %in% D11_met$Group,]

rownames(D11_share) == D11_met$Group


D14_met <- mastermeta.feces %>% filter(day == 14)
D14_share <- share.rare[rownames(share.rare) %in% D14_met$Group,]

rownames(D14_share) == D14_met$Group

###


D0_sow_adon <- adonis(formula = D0_share ~ treatment + Sow + Sex, data = D0_met, method = 'bray')
D4_sow_adon <- adonis(formula = D4_share ~ treatment + Sow + Sex, data = D4_met, method = 'bray')
D7_sow_adon <- adonis(formula = D7_share ~ treatment + Sow + Sex, data = D7_met, method = 'bray')
D11_sow_adon <- adonis(formula = D11_share ~ treatment + Sow + Sex, data = D11_met, method = 'bray')
D14_sow_adon <- adonis(formula = D14_share ~ treatment + Sow + Sex, data = D14_met, method = 'bray')


D0_sow_adon <- as.data.frame(D0_sow_adon$aov.tab)[1:3,] %>% select(F.Model, `Pr(>F)`) %>% mutate(day=0, var=rownames(.))
D4_sow_adon <- as.data.frame(D4_sow_adon$aov.tab)[1:3,] %>% select(F.Model, `Pr(>F)`) %>% mutate(day=4, var=rownames(.))
D7_sow_adon <- as.data.frame(D7_sow_adon$aov.tab)[1:3,] %>% select(F.Model, `Pr(>F)`) %>% mutate(day=7, var=rownames(.))
D11_sow_adon <- as.data.frame(D11_sow_adon$aov.tab)[1:3,] %>% select(F.Model, `Pr(>F)`) %>% mutate(day=11, var=rownames(.))
D14_sow_adon <- as.data.frame(D14_sow_adon$aov.tab)[1:3,] %>% select(F.Model, `Pr(>F)`) %>% mutate(day=14, var=rownames(.))

sow_adon <- bind_rows(list(D0_sow_adon, 
               D4_sow_adon, 
               D7_sow_adon, 
               D11_sow_adon, 
               D14_sow_adon))


sow_adon <- sow_adon %>% mutate(pval=round(p.adjust(`Pr(>F)`, method = 'BH'),digits = 3))

sow_adon <- sow_adon %>% mutate(plab=ifelse(pval < 0.05, pval, NA))


sow_adon %>% ggplot(aes(x=day, y=F.Model, color=var,fill=var, label=plab)) + 
  geom_point()+ geom_line() + geom_text_repel() + ggtitle('Do covariates effect microbial community structure?', 
                                                          subtitle = 'significant pvalues shown')

# We see that treatment effect again, but now we can see that sow has a significant influence
# on the community structure similarities.  This effect is strongest early on and then wanes.
# We anticipated this and did our best to evenly distribute piglets from each litter into the 
# different treatment groups.

# this makes sense, as the early microbiota is what they got from their mothers.
# then as they adapt to the new diet these initial communities are replaced.



library(ggrepel)

mastermeta.feces$set

mastermeta.feces$set2 <- paste(mastermeta.feces$day,
                               mastermeta.feces$tissue,
                               mastermeta.feces$treatment,
                               mastermeta.feces$pen, sep = '_')

PWadon_pen <- pairwise.adonis(share.rare, mastermeta.feces$set2)

PWadon_pen$pairs

# this grabs all penA vs penB comparisons in each room at all timepoints

test <- grep("(^[0-9]+)_Feces_(.*)_([0-9]+)[AB] vs \\1_Feces_\\2_\\3.*", PWadon_pen$pairs)

pen_comps <- PWadon_pen[test,]
pen_comps$p.adjusted <- p.adjust(pen_comps$p.value, method = 'BH') %>% round(digits = 2)

pen_comps$day <- gsub('([0-9]+)_.*','\\1',pen_comps$pairs)
pen_comps$room <- gsub('.*_.*_.*_([123])[AB] vs .*_.*_.*_(.*)','\\1',pen_comps$pairs)

pen_comps %>% ggplot(aes(x=day, y=F.Model, color=room, group=room)) + geom_line()

pen_res <- pen_comps %>% mutate(comparision=paste(room, 'A',' vs ', room, 'B', sep = ''), 
                     F.Model=round(F.Model, digits = 3), 
                     day= as.numeric(day)) %>% arrange(day) %>% 
  select(day, comparision, F.Model, p.adjusted)

write_csv(pen_res, 'pen_permanova.csv')

# after fdr adjusting pvalues in these comparisons, no sig diffs detected between pens

#### checking for diffs in sow, pen, etc



PWadon <- pairwise.adonis(share.rare, mastermeta.feces$set)



fecesVSfeces <- grep(".*_Feces_.* vs .*_Feces_.*", PWadon$pairs)

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
  scale_color_manual(values = c('#00BA38', '#619CFF'), labels = c('Inject', 'Feed'), name='Treatment')
p2 + 
  ylab('PERMANOVA F vs NM\n(difference relative to NM)')


##### feces ordination ######  

#note that phyloseq now also has a function for plotting ellipses
#note that phyloseq's ellipses are stinky and I havent figured out how to use them well

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


# these are for plotting each day over the total timecourse

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
  scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))+
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
        axis.line = element_blank()) + scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))+
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
        axis.line = element_blank())+scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))+
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
        axis.line = element_blank()) + scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))+
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
        axis.line = element_blank()) + scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))+
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
        axis.line = element_blank()) + scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))+
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
  ggtitle('NMDS: Colon, Day 4', subtitle = 'stress: 0.169')+ scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))

ggplot(data = filter(all.metanmds, tissue =='Colon'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_7, tissue == 'Colon'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_7, tissue == 'Colon'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_7, tissue == 'Colon'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Colon, Day 7', subtitle = 'stress: 0.169')+scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))

ggplot(data = filter(all.metanmds, tissue =='Colon'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_14, tissue == 'Colon'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_14, tissue == 'Colon'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_14, tissue == 'Colon'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Colon, Day 14', subtitle = 'stress: 0.169') + scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))

########

ggplot(data = filter(all.metanmds, tissue =='Ileum'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_4, tissue == 'Ileum'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_4, tissue == 'Ileum'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_4, tissue == 'Ileum'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Ileum, Day 4', subtitle = 'stress: 0.169')+scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))

ggplot(data = filter(all.metanmds, tissue =='Ileum'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_7, tissue == 'Ileum'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_7, tissue == 'Ileum'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_7, tissue == 'Ileum'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Ileum, Day 7', subtitle = 'stress: 0.169') + scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))

ggplot(data = filter(all.metanmds, tissue =='Ileum'), aes(x=MDS1, y=MDS2)) + geom_point(color='grey57') + 
  geom_point(data = filter(all.metanmds_14, tissue == 'Ileum'), aes(color = set)) + 
  geom_path(data=filter(all.df_ell_14, tissue == 'Ileum'), aes(x=NMDS1, y=NMDS2, color = group)) + 
  geom_segment(data=filter(all.metanmds_14, tissue == 'Ileum'), aes(x=MDS1, xend=centroidX, y=MDS2, yend=centroidY, color=set)) + 
  ggtitle('NMDS: Ileum, Day 14', subtitle = 'stress: 0.169') + scale_color_discrete(labels=c('NM', 'Inject', 'Feed'))


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
  scale_color_manual(values = c('#00BA38', '#619CFF'), labels=c('Inject', 'Feed')) + ylim(1,2.5)
p3 + ggtitle('PERMANOVA comparisions to Control group', subtitle = 'Colon') + 
  ylab('PERMANOVA F vs control\n(difference relative to control)')

p4 <- ggplot(PWadon.tocont.ileum, aes(x=day, y=F.Model, color=treatment)) +
  geom_path(aes(group=treatment), size=1.3) + 
  scale_fill_manual(values = c('#00BA38', '#619CFF')) + 
  geom_label(aes(label=p.value, fill=treatment), color='black', show.legend = FALSE) +
  scale_color_manual(values = c('#00BA38', '#619CFF'), labels=c('Inject', 'Feed'))#+ ylim(0,2.5)
p4 + ggtitle('PERMANOVA comparisions to Control group', subtitle = 'Ileum') + 
  ylab('PERMANOVA F vs control\n(difference relative to control)')


################### Deseq2 #####################



library(DESeq2)
library(phyloseq)

getwd()
otu <- import_mothur(mothur_shared_file = 'FS1.V4.shared')
taxo <- import_mothur(mothur_constaxonomy_file = 'FS1.V4.taxonomy')

#meta <- read.table(file = 'V4.metadata.txt', sep = '\t', header = TRUE)

rownames(meta) <- meta$Group
meta$set <- paste(meta$day, meta$tissue, meta$treatment, sep = '_')
meta$day <- factor(meta$day, levels = c(0,4,7,11,14))
meta$treatment2 <- ifelse(meta$treatment == 'oral', 'Feed', 
                          ifelse(meta$treatment == 'inject', 'Inject', 
                                 ifelse(meta$treatment == 'control', 'NM', NA)))

meta$treatment2 <- factor(meta$treatment2, levels = c('NM', 'Inject', 'Feed'))

phy_meta <- sample_data(meta) 
FS1 <- phyloseq(otu, taxo)
FS1 <- merge_phyloseq(FS1, phy_meta)                       # combines the metadata with this phyloseq object
colnames(tax_table(FS1)) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus')
FS1 <- subset_samples(FS1, experiment == 'X1')


FS1 <- prune_samples(sample_sums(FS1) > 2000, FS1)  # This removes samples that have fewer than 700 sequences associated with them.
FS1 <- prune_taxa(taxa_sums(FS1) > 5, FS1)        # removes OTUs that occur less than 10 times globally

###
# Added thing that removes mitochondria and chloroplast reads

prune_vec <- !grepl('Mitochondria',tax_table(FS1)[,6]) & !grepl('Chloroplast',tax_table(FS1)[,6])

FS1 <- prune_taxa(taxa = prune_vec, x = FS1)

##### function attempt #####

# this is a large unweildy function that subsets the data according to DAY and TIS
# runs it through the DESeq2 pipe and generates the types of plots you were looking for earlier

Deseq_plots2 <- function(phyloseq_obj, DAY, TIS, cooks_cut = TRUE,
                         pAdjustMethod = 'BH', shrink_type = 'normal', 
                         tax_lab, alpha=0.05, alpha_filt=0.05){
  tax_lab <- enquo(tax_lab)
  print(DAY)
  print(TIS)
  FS1 <- prune_samples(x = phyloseq_obj, samples = phyloseq_obj@sam_data$day == DAY & phyloseq_obj@sam_data$tissue == TIS)
  
  FS1 <- prune_taxa(taxa_sums(FS1) > 1, FS1)
  print(nsamples(FS1))
  
  FS1@sam_data$treatment2
  
  FS1.De <- phyloseq_to_deseq2(FS1, ~ treatment2)
  
  FS1.De <- DESeq(FS1.De, test = "Wald", fitType = "parametric")
  
  
  print(resultsNames(FS1.De))
  
  ########## pairwise comparisons ##########
  
  res_Feed_v_NM <- results(FS1.De, contrast = c('treatment2', 'Feed', 'NM'), alpha = alpha, cooksCutoff = cooks_cut, pAdjustMethod = pAdjustMethod)
  res_Feed_v_NM <- lfcShrink(coef = 3, res=res_Feed_v_NM, dds = FS1.De, type = shrink_type)
  
  res_Feed_v_NM <- res_Feed_v_NM[which(res_Feed_v_NM$padj < alpha_filt),]
  
  
  
  # res_FvC = results(FS1.De, contrast=c("set","D4_Feces_oral","D4_Feces_control"),cooksCutoff = FALSE, pAdjustMethod = 'BH')
  
  sigtab_Feed = res_Feed_v_NM[which(res_Feed_v_NM$padj < .05), ]
  print(nrow(sigtab_Feed))
  if (nrow(sigtab_Feed) != 0){
    # browser()
    sigtab_Feed = cbind(as(sigtab_Feed, "data.frame"), as(tax_table(FS1)[rownames(sigtab_Feed), ], "matrix"))
    #format(sigtab$padj, scientific = TRUE)
    # sigtab_Feed$newp <- format(round(sigtab_Feed$padj, digits = 5), scientific = TRUE)
    sigtab_Feed$Treatment <- ifelse(sigtab_Feed$log2FoldChange >=0, "Feed", "NM")
    
    sigtab_Feed$Treatment <- ifelse(sigtab_Feed$log2FoldChange >=0, "Enriched in Feed", "Depleted in Feed")
    # sigtab_Feed <- sigtab_Feed %>% select(log2FoldChange, Kingdom, Phylum, Class, Order, Treatment, Genus)
    # # colnames(sigtab_Feed)[1]<-"log2FoldChange"
    sigtab_Feed$OTU<-rownames(sigtab_Feed)
    
    
  }
  
  
  # Inject vs NM
  res_Inject_v_NM = results(FS1.De, contrast=c("treatment2", "Inject", "NM"), alpha=0.05, cooksCutoff = cooks_cut, pAdjustMethod = pAdjustMethod)
  res_Inject_v_NM <- lfcShrink(coef = 2, res=res_Inject_v_NM, dds = FS1.De, type = shrink_type)
  
  
  sigtab_Inject = res_Inject_v_NM[which(res_Inject_v_NM$padj < .05), ]
  print(nrow(sigtab_Inject))
  if (nrow(sigtab_Inject) != 0){
    # browser()
    sigtab_Inject = cbind(as(sigtab_Inject, "data.frame"), as(tax_table(FS1)[rownames(sigtab_Inject), ], "matrix"))
    sigtab_Inject$Treatment <- ifelse(sigtab_Inject$log2FoldChange >=0, "Inject", "NM")
    
    
    #need to pull these into a combined table and label the log2FoldChange by treatment and changing the labels to show which control is which
    
    sigtab_Inject$Treatment <- ifelse(sigtab_Inject$log2FoldChange >=0, "Enriched in Inject", "Depleted in Inject")
    # sigtab_Inject <- sigtab_Inject %>% select(log2FoldChange, Kingdom, Phylum, Class, Order, Treatment, Genus)
    # # colnames(sigtab_Inject)[1]<-"log2FoldChange"
    sigtab_Inject$OTU<-rownames(sigtab_Inject)
    
  }
  
  
  
  if (nrow(sigtab_Inject)!=0 & nrow(sigtab_Feed) != 0){
    # browser()
    D4.sigs <- rbind(sigtab_Feed, sigtab_Inject)
  }
  
  if (nrow(sigtab_Inject) == 0 & nrow(sigtab_Feed) != 0){
    D4.sigs <- sigtab_Feed
  }
  
  if (nrow(sigtab_Inject) !=0 & nrow(sigtab_Feed) == 0){
    D4.sigs <- sigtab_Inject
  }
  
  
  
  if (nrow(sigtab_Inject)==0 & nrow(sigtab_Feed) == 0){
    warning('NO SIG DIF FEATURES IN EITHER GROUP TO CONTROL')
    return(NULL)
  }
  
  D4.sigs$Treatment<-factor(D4.sigs$Treatment, levels=c("Enriched in Feed",
                                                        "Enriched in Inject",
                                                        "Depleted in Feed",
                                                        "Depleted in Inject"), ordered=TRUE)
  
  
  D4.sigs <- D4.sigs[order(D4.sigs$log2FoldChange),]
  #HERE! D4.sigs$OTU doesnt exist
  # D4.sigs$OTU <- rownames(D4.sigs)
  D4.sigs$OTU <- reorder(D4.sigs$OTU, D4.sigs$log2FoldChange)
  D4.sigs$day <- DAY
  # browser()
  D4_labs <- D4.sigs %>% select(OTU, !!tax_lab) %>% unique()
  
  # browser()
  res <- ggplot(D4.sigs, aes(x=OTU,
                             y=log2FoldChange,
                             fill = Treatment)) + 
    geom_bar(stat='identity', position=position_dodge(preserve = 'single')) + coord_flip() + 
    geom_text(data=D4_labs, aes(x=OTU, y=0, label=!!tax_lab),
              inherit.aes = FALSE, size=4, fontface='italic', vjust=-.2) +
    scale_fill_manual(values = c('royalblue', 'skyblue', 'forestgreen', 'lightgreen'), drop=FALSE) +xlab('Order') + 
    theme(axis.text.y=element_blank(), 
          axis.ticks.y = element_blank()) #+ ylim(-5, 5)
  
  return(list(plot=res, result=D4.sigs))
  
}

# feces OTU level, all days

days <- c(0,4,7,11,14)

final_res <- list()
OTU_res <- list()

for (i in 1:length(days)){
  day_lab <- paste('day', days[i], sep = '_')
  print(day_lab)
  OTU_res[[day_lab]] <- Deseq_plots2(FS1, DAY=days[i], TIS='Feces', tax_lab = Family)
  if (!is.null(OTU_res[[day_lab]])){
    final_res[[day_lab]] <- OTU_res[[day_lab]]$result
  }
  
}
# final_res[[1]]$Treatment

feces_diff_OTUs <- bind_rows(final_res)

feces_diff_OTUs <- feces_diff_OTUs %>% mutate(tissue='Feces', 
                           level='OTU')

### colon OTU level all days ###
days <- c(4,7,14)

final_res <- list()
OTU_res <- list()
for (i in 1:length(days)){
  day_lab <- paste('day', days[i], sep = '_')
  print(day_lab)
  OTU_res[[day_lab]] <- Deseq_plots2(FS1, DAY=days[i], TIS='Colon', tax_lab = Family)
  if (!is.null(OTU_res[[day_lab]])){
    final_res[[day_lab]] <- OTU_res[[day_lab]]$result
  }
  
}

colon_diff_OTUs <- bind_rows(final_res)
colon_diff_OTUs <- colon_diff_OTUs %>% mutate(tissue='Colon', 
                           level='OTU')

# colon_diff_OTUs$tissue <- 'Colon'
# colon_diff_OTUs$tissue <- 'Colon'
colon_diff_OTUs %>% group_by(day) %>% tally()

# feces gloms, all days
days <- c(4,7,11,14)
tax_levs <- c('Order', 'Family', 'Genus')

final_res <- list()
OTU_res <- list()
lev_res <- list()

for (x in 1:length(tax_levs)){
  FS1.glom <- tax_glom(physeq = FS1, taxrank = tax_levs[x])
  for (i in 1:length(days)){
    day_lab <- paste('day', days[i], sep = '_')
    print(day_lab)
    OTU_res[[day_lab]] <- Deseq_plots2(FS1.glom, DAY=days[i], TIS='Feces', tax_lab = Family)
    if (!is.null(OTU_res[[day_lab]])){
      lev_res[[day_lab]] <- OTU_res[[day_lab]]$result
    }
    final_res[[tax_levs[[x]]]] <- lev_res
  }
  
}

names(final_res)

feces_order <- bind_rows(final_res[['Order']])
feces_family <- bind_rows(final_res[['Family']])
feces_genus <- bind_rows(final_res[['Genus']])

feces_order <- feces_order %>% mutate(tissue='Feces', 
                                      level='Order')

feces_family <- feces_family %>% mutate(tissue='Feces', 
                                        level='Family')

feces_genus <- feces_genus %>% mutate(tissue='Feces', 
                                      level='Genus')




#### colon gloms, all days ###

days <- c(4,7,14)
tax_levs <- c('Order', 'Family', 'Genus')

final_res <- list()
OTU_res <- list()
lev_res <- list()

for (x in 1:length(tax_levs)){
  FS1.glom <- tax_glom(physeq = FS1, taxrank = tax_levs[x])
  for (i in 1:length(days)){
    day_lab <- paste('day', days[i], sep = '_')
    print(day_lab)
    OTU_res[[day_lab]] <- Deseq_plots2(FS1.glom, DAY=days[i], TIS='Colon', tax_lab = Family)
    if (!is.null(OTU_res[[day_lab]])){
      lev_res[[day_lab]] <- OTU_res[[day_lab]]$result
    }
    final_res[[tax_levs[[x]]]] <- lev_res
  }
  
}

names(final_res)

colon_order <- bind_rows(final_res[['Order']])
colon_family <- bind_rows(final_res[['Family']])
colon_genus <- bind_rows(final_res[['Genus']])

colon_order <- colon_order %>% mutate(tissue='Colon', 
                                      level='Order')

colon_family <- colon_family %>% mutate(tissue='Colon', 
                                      level='Family')

colon_genus <- colon_genus %>% mutate(tissue='Colon', 
                                      level='Genus')



master_difabund <- bind_rows(list(feces_diff_OTUs, feces_genus, feces_family, feces_order, 
               colon_diff_OTUs, colon_genus, colon_family, colon_order))

write_csv(master_difabund, "master_diffabund.csv")
#
#

#Figs you were looking for
FS1.glom <- tax_glom(FS1, taxrank = "Order")

# DESeq_plots2 returns a list, the first item in the list is a ggplot2 object
# the second is the results dataframe

fec_d4 <- Deseq_plots2(FS1.glom, DAY=4, TIS='Feces', tax_lab = Order)
fec_d7 <- Deseq_plots2(FS1.glom, DAY=7, TIS='Feces', tax_lab = Order)


col_d4 <- Deseq_plots2(FS1.glom, DAY=4, TIS='Colon', tax_lab = Order)
# your colon D7 code was failing because there are no differentially abundant orders between
# control and inject at that time
col_d7 <- Deseq_plots2(FS1.glom, DAY=7, TIS='Colon', tax_lab = Order)


fec_d4[[1]] + ylim(-5,5)  # D4 feces order
fec_d7[[1]] + ylim(-6,6)  # D7 feces order

col_d4[[1]] + ylim(-5,5)  # D4 colon order
col_d7[[1]] + ylim(-8,6)  # D7 colon order


# ggplot(D4.sigs, aes(x=OTU,
#                     y=log2FoldChange,
#                     fill = Treatment)) + 
#   geom_bar(stat='identity', position=position_dodge(preserve = 'single')) + coord_flip() + 
#   geom_text(data=D4_labs, aes(x=OTU, y=0, label=tax_lab),
#             inherit.aes = FALSE, size=4, fontface='italic', vjust=-.2) +
#   scale_fill_manual(values = c('royalblue', 'skyblue', 'forestgreen', 'lightgreen'), drop=FALSE) +xlab('Order') + 
#   theme(axis.text.y=element_blank(), 
#         axis.ticks.y = element_blank()) #+ ylim(-5, 5)
