


###### ABX conc. data  #########


getwd()
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

tis.melt$concentration[tis.melt$concentration == 0] <- 1
tis.melt$Treatment <- ifelse(tis.melt$treatment == 'control', 'NM', 
                             ifelse(tis.melt$treatment == 'inject', 'Inject', 
                                    ifelse(tis.melt$treatment == 'oral', 'Feed', NA)))

tis.melt$Tissue <- ifelse(tis.melt$tissue == 'Serum', 'Plasma', 
                             ifelse(tis.melt$tissue == 'feces', 'Feces',
                                    ifelse(tis.melt$tissue == 'Nasal', 'Nasal', 
                                           ifelse(tis.melt$tissue == 'Ileum', 'Ileum', NA))))






tis.melt$Treatment <- factor(tis.melt$Treatment, levels = c('NM', 'Inject', 'Feed'))

tis.melt %>% filter(Treatment %in% c('Inject', 'Feed')) %>%
  filter(Tissue %in% c('Feces', 'Plasma')) %>% 
  filter(day != 0) %>% 
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=Treatment)) +
  geom_boxplot() +
  ylab('concentration ng/mL') + scale_y_log10() +
  facet_wrap(~Tissue, scales = 'free')+
  theme_bw() + scale_fill_manual(values = c('#00BA38', '#619CFF'))

tis.melt %>% filter(Treatment %in% c('Inject', 'Feed')) %>%
  #filter(Tissue %in% c('Feces', 'Plasma')) %>% 
  filter(day != 0) %>% 
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=Treatment)) +
  geom_boxplot() +
  ylab('concentration ng/mL') + scale_y_log10() +
  facet_wrap(~Tissue, scales = 'free')+
  theme_bw() + scale_fill_manual(values = c('#00BA38', '#619CFF'))



#
tis.melt$Day <- as.numeric(as.character(tis.melt$day))

tis.melt %>% filter(Treatment %in% c('Inject', 'Feed')) %>%
  filter(Tissue %in% c('Feces', 'Plasma')) %>% 
  filter(day != 0) %>% 
  ggplot(aes(x=Day, y=concentration, fill=Treatment, color=Treatment)) +
  geom_point(shape=21, alpha=.5) + geom_smooth() +
  ylab('concentration ng/mL') + scale_y_log10() +
  facet_wrap(~Tissue, scales = 'free')+
  theme_bw() + scale_fill_manual(values = c('#00BA38', '#619CFF'))+
  scale_color_manual(values = c('#00BA38', '#619CFF'))









tis.melt %>% filter(tissue %in% c('Serum', 'feces')) %>%
  filter(day != 0) %>% 
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=Treatment)) +
  geom_boxplot() +
  ylab('concentration ng/mL') +# scale_y_log10() +
  facet_wrap(~Tissue, scales = 'free')+
  theme_bw()

tis.melt %>% filter(tissue %in% c('Serum', 'feces')) %>%
  filter(day != 0) %>% 
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=Treatment)) +
  geom_boxplot() +
  ylab('concentration ng/mL') + scale_y_continuous(trans = 'log2') +
  facet_wrap(~Tissue, scales = 'free')+
  theme_bw()



tis.melt %>% group_by(Treatment, day) %>% 
  summarise(mean_conc = mean(concentration)) %>% 
  spread(key=Treatment, value=mean_conc, -day)

scales::trans_new(transform = log2, inverse = exp)

tis.melt %>% filter(tissue == 'Serum') %>%                                 
  ggplot(aes(x=day, y=concentration, group=dayXtreat, fill=treatment)) +
  geom_boxplot() +
  ylab('concentration ng/mL') +
  ggtitle('Oxytet concentrations: Serum') +
  theme_bw()


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
