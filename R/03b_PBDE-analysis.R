#### Load libraries ####

library(MASS)
library(tidyverse)
library(here)
library(lme4)
library(RLRsim)
library(ggdist)
library(tidyquant)
library(dplyr)

# set directories
clean_data_dir <- here("data", "clean")

# get excel file name
totals <- "totals_all.csv"

# load dataframe
mussel_df <- read_csv(here(clean_data_dir, totals), show_col_types = FALSE)

#abs of dry value for now
mussel_df$dry_value <- abs(mussel_df$dry_value)

#years as year from 2010
#mussel_df$year <- mussel_df$year - 2010

#as category
mussel_df$year <- as.factor(mussel_df$year)
mussel_df$wria_nr <- as.factor(mussel_df$wria_nr)
mussel_df$lio_areas <- as.factor(mussel_df$lio_areas)

#remove Penn Cove reference samples
mussel_df <- mussel_df[!is.na(mussel_df$fund_source),]
mussel_df <- subset(mussel_df, site_name!="Penn Cove Reference")

#remove samples from outside Puget sound
mussel_df <- subset(mussel_df, longitude > -123.5)


#create simple regions based on longitude
# mussel_df <- mussel_df %>% mutate(region = case_when(longitude > -122.75 & latitude > 47.6 ~  'Northeast',
#                                                      longitude > -122.75 & latitude < 47.6 ~  'Southeast',
#                                                      longitude < -122.75 ~ 'West'
# ))
# mussel_df$region <- as.factor(mussel_df$region)



#create analyte dataframes
PB_df <- mussel_df %>%
  filter(mussel_df$analyte == "SumPBDEs11")
PB_df <- PB_df[complete.cases(PB_df),]

#### subregion things ####
table(PB_df$year)
table(PB_df$year, PB_df$lio_areas)
table(PB_df$year, PB_df$wria_nr)


#### means by subregions ####
all_mean <- ddply(PB_df, c("year"), summarise,
                  mean = mean(dry_value))
#change is + - -
lio_mean <-ddply(PB_df, c("lio_areas","year"), summarise,
                 mean = mean(dry_value))
wria_mean <- ddply(PB_df, c("wria_nr","year"), summarise,
                   mean = mean(dry_value))

#### raincloud plots ####
##over years
pdf("PBDEs_Over_Years.pdf")
PB_df %>% 
  ggplot(aes(x = year, y = log10(dry_value))) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1)+
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  xlab("Year") +
  ylab("log10(dry value)")+
  ggtitle("PBDEs Across Years")
dev.off()
## over LIOs
pdf("PBDEs_Over_LIOs.pdf")
PB_df %>% 
  ggplot(aes(x = lio_areas, y = log10(dry_value))) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1)+
   theme(axis.text.x=element_text(angle=315,hjust = 0.01)) +
  xlab("LIO Area") +
  ylab("log10(dry value)")+
  ggtitle("PBDEs Across LIOs")
dev.off()
##over wria_nr
pdf("PBDEs_Over_WRIAs.pdf")
PB_df %>% 
  ggplot(aes(x = wria_nr, y = log10(dry_value))) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  ggdist::stat_dots(side = "left", dotsize = .05, justification = 1.1, binwidth = .1)+
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  xlab("WRIA") +
  ylab("log10(dry value)")+
  ggtitle("PBDEs Across WRIAs")
dev.off()
#### WRIA plot loop ####
plot_list = list()
for (i in unique(PB_df$wria_nr)){
  wria_subset <- PB_df[PB_df$wria_nr==i,]
  p = wria_subset %>%
    ggplot(aes(x = year, y = log10(dry_value))) +
    ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .1, outlier.shape = NA) +
    ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
    ggtitle(paste("PBDEs - WRIA #",i)) +
    xlab("Year") +
    ylab("log10(dry value)")
  plot_list[[i]] = p
}
#end plot loop
#start plot to pdf loop

pdf("PBDE_WRIA_plots.pdf")
for (i in unique(PB_df$wria_nr)){
  print(plot_list[[i]])
}
dev.off()

#### LIO plot loop ####
plot_list = list()
for (i in unique(PB_df$lio_areas)){
  wria_subset <- PB_df[PB_df$lio_areas==i,]
  p = wria_subset %>%
    ggplot(aes(x = year, y = log10(dry_value))) +
    ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .1, outlier.shape = NA) +
    ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
    theme(axis.text.x=element_text(angle=270,hjust=1)) +
    xlab("LIO") +
    ylab("log10(dry value)")+
    ggtitle(paste("PBDEs - LIO -",i))
  plot_list[[i]] = p
}
#end plot loop
#start plot to pdf loop

pdf("PBDE_LIO_plots.pdf")
for (i in unique(PB_df$lio_areas)){
  print(plot_list[[i]])
}
dev.off()

#### 2020 only plots ####
PB20_df <- PB_df %>% filter(PB_df$year == 2020)
pdf("2020_PBDEs_Over_WRIAs.pdf")
PB20_df %>% 
  ggplot(aes(x = wria_nr, y = log10(dry_value))) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1)+
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  xlab("WRIA") +
  ylab("log10(dry value)")+
  ggtitle("2020 - PBDEs Across WRIAs")
dev.off()



#### anova ####
PB_aov <- aov(PB_df$dry_value ~ PB_df$wria_nr, data = PB_df)



