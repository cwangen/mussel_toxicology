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
PC_df <- mussel_df %>%
  filter(mussel_df$analyte == "SumPCBs2x17")
PC_df <- PC_df[complete.cases(PC_df),]




#### Fit model with RE for lat and lon and plots ####
PC_LMM <- lmer(log(dry_value) ~
                 time +
                 year +
                 mean_is_au +
                 year:lio_areas +
                 (1|longitude),
               data = PC_df)



#### QQ plots ####
## set plot area
par(mai = c(0.9, 0.9, 0.6, 0.1),
    omi = c(0, 0, 0, 0),
    mfrow = c(1,3), cex.lab = 1.2)

## residuals vs fitted
plot(fitted(PC_LMM), residuals(PC_LMM), las = 1, pch = 1,
     xlab = "Fitted", ylab = "Residuals",
     main = "Residuals vs Fitted")
abline(h=0, lty = "dashed")


## residuals
qqnorm(residuals(PC_LMM), main = "QQ plot (residuals)", las = 1, pch = 1)
qqline(residuals(PC_LMM))

## Random effects
qqnorm(unlist(ranef(PC_LMM)), main = "QQ plot (RE's)", las = 1, pch = 1)
qqline(unlist(ranef(PC_LMM)))


#### subregion things ####
table(PC_df$year)
table(PC_df$year, PC_df$lio_areas)
table(PC_df$year, PC_df$wria_nr)


#### means by subregions ####
all_mean <- ddply(PC_df, c("year"), summarise,
                  mean = mean(dry_value))
#change is + - -
lio_mean <-ddply(PC_df, c("lio_areas","year"), summarise,
                 mean = mean(dry_value))
wria_mean <- ddply(PC_df, c("wria_nr","year"), summarise,
                   mean = mean(dry_value))

#### raincloud plots ####
##over years
pdf("PCBs_Over_Years.pdf")
PC_df %>% 
  ggplot(aes(x = year, y = log10(dry_value))) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
theme(axis.text.x=element_text(angle=270,hjust=1)) +
  xlab("Year") +
  ylab("log10(dry value)")+
ggtitle("PCBs Across Years")
dev.off()
## over LIO
pdf("PCBs_Over_LIOs.pdf")
PC_df %>% 
  ggplot(aes(x = lio_areas, y = log10(dry_value))) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1)+
  theme(axis.text.x=element_text(angle=315,hjust = 0.01)) +
  xlab("LIO Area") +
  ylab("log10(dry value)")+
ggtitle("PCBs Across LIOs")
dev.off()
##over wria_nr
pdf("PCBs_Over_WRIAs.pdf")
PC_df %>% 
  ggplot(aes(x = wria_nr, y = log10(dry_value))) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  ggdist::stat_dots(side = "left", dotsize = 0.5, justification = 1.1, binwidth = .1) +
  theme(axis.text.x=element_text(angle=270,hjust=1)) +
  xlab("WRIA") +
  ylab("log10(dry value)")+
ggtitle("PCBs Across WRIAs")
dev.off()
#### WRIA plot loop ####
plot_list = list()
for (i in unique(PC_df$wria_nr)){
  wria_subset <- PC_df[PC_df$wria_nr==i,]
  main = paste("WRIA #",i)
  p = wria_subset %>%
    ggplot(aes(x = year, y = log10(dry_value))) +
    ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .1, outlier.shape = NA) +
    ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
    ggtitle(paste("PCBs - WRIA #",i)) +
    xlab("Year") +
    ylab("log10(dry value)")
  plot_list[[i]] = p
}
#end plot loop
#start plot to pdf loop

pdf("PCBs_WRIA_plots.pdf")
for (i in unique(PC_df$wria_nr)){
  print(plot_list[[i]])
}
dev.off()

#### LIO plot loop ####
plot_list = list()
for (i in unique(PC_df$lio_areas)){
  wria_subset <- PC_df[PC_df$lio_areas==i,]
  p = wria_subset %>%
    ggplot(aes(x = year, y = log10(dry_value))) +
    ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
    geom_boxplot(width = .1, outlier.shape = NA) +
    ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1) +
    theme(axis.text.x=element_text(angle=270,hjust=1)) +
    xlab("LIO") +
    ylab("log10(dry value)")+
    ggtitle(paste("PCBs - LIO -",i))
  plot_list[[i]] = p
}
#end plot loop
#start plot to pdf loop

pdf("PCBs_LIO_plots.pdf")
for (i in unique(PC_df$lio_areas)){
  print(plot_list[[i]])
}
dev.off()

#### 2020 only plots ####
PC20_df <- PC_df %>% filter(PC_df$year == 2020)
pdf("2020_PCBs_Over_WRIAs.pdf")
PC20_df %>% 
  ggplot(aes(x = wria_nr, y = log10(dry_value))) +
  ggdist::stat_halfeye(adjust = .5, width = .3, .width = 0, justification = -.3, point_colour = NA) + 
  geom_boxplot(width = .1, outlier.shape = NA) +
  ggdist::stat_dots(side = "left", dotsize = .1, justification = 1.1, binwidth = .1)
                    theme(axis.text.x=element_text(angle=270,hjust=1)) +
                      xlab("WRIA") +
                      ylab("log10(dry value)")+
  ggtitle("2020 - PCBs Across WRIAs")
dev.off()


#### anova ####
PC_aov <- aov(PC_df$dry_value ~ PC_df$wria_nr, data = PC_df)



