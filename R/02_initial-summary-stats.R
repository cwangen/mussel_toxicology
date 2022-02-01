##This script performs data summaries on the 2013-2015 mussel data.
##To be added and corrected as needed for additional data.

#### setup ####

#load libraries
library(here)
library(tidyverse)
library(ggplot2)
library(ggmap)
library(sf)
library(mapview)

# set directories 
clean_data_dir <- here("data", "clean")

#get excel file name
totals <- "totals_1316.csv"

#load dataframe
m1316_df <- read_csv(here(clean_data_dir, totals)) #, show_col_types = FALSE) 

#### reorganization as needed ####

#remove test sites 
m1316_exp <- m1316_df[!is.na(m1316_df$latitude),]

#remove rows with missing wet_value (there is only one)
m1316_exp <- m1316_exp[!is.na(m1316_exp$wet_value),]

### tables and numbers ###

#number of sites per year
site_totals <-table(m1316_exp$analyte, m1316_exp$year)
print(site_totals)
site_totals <-data.frame(unclass(site_totals))

#number of overlapping sites
sites2013 <- subset(m1316_exp, m1316_exp$year == "2013")
sites2016 <- subset(m1316_exp, m1316_exp$year == "2016")
overlap <- sites2013$site_name %in% sites2016$site_name
total_overlap <- sum(overlap, na.rm = TRUE) / 3

#days samples spent in water
todiff <- m1316_exp[order(m1316_exp$deployment_date),]
difftime(todiff$retrieval_date, todiff$deployment_date)

#analyte min/max/mean table
by_analyte <- m1316_exp %>% group_by(analyte)
analyte_table <- by_analyte %>% summarise(avg_wv = mean(wet_value),
                         min_wv = min(wet_value),
                         max_wv = max(wet_value))

#### maps ####

#map of sample locations
m1316_l <- m1316_exp %>% filter(m1316_exp$analyte == "lipids")
lon <- m1316_l$longitude
lat <- m1316_l$latitude
df <- as.data.frame(cbind(lon, lat))
qmplot(lon, lat, data = df, maptype = "toner-lite", color = I("blue"), size = I(0.2))

###histograms of wet_values for different analytes ####

#get wet values
m1316_lipids <- subset(m1316_exp, m1316_exp$analyte == "lipids")
m1316_PDBE <- subset(m1316_exp, m1316_exp$analyte == "SumPBDEs11")
m1316_PCB <- subset(m1316_exp, m1316_exp$analyte == "SumPCBs2x17")

#histograms combined years
hist(m1316_lipids$wet_value)
hist(m1316_PDBE$wet_value, breaks = 0:16)
hist(m1316_PCB$wet_value, breaks = 0:80)

#by year both histograms (13 is blue, 16 is pink)
m1316_lipids16 <- subset(m1316_lipids, m1316_lipids$year == "2016")
m1316_lipids13 <- subset(m1316_lipids, m1316_lipids$year == "2013")
bins = seq(0, 1.5, by = 0.1)
hg16 <- hist(m1316_lipids16$wet_value, breaks = bins, plot = FALSE)
hg13 <- hist(m1316_lipids13$wet_value, breaks = bins, plot = FALSE)
c1 <- rgb(173, 216, 230, max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255, 192, 203, max = 255, alpha = 80, names = "lt.pink")
plot(hg13, col = c1)
plot(hg16, col = c2, add = TRUE)

m1316_PDBE16 <- subset(m1316_PDBE, m1316_PDBE$year == "2016")
m1316_PDBE13 <- subset(m1316_PDBE, m1316_PDBE$year == "2013")
bins = seq(0, 16, by = 0.5)
hg16 <- hist(m1316_PDBE16$wet_value, breaks = bins, plot = FALSE)
hg13 <- hist(m1316_PDBE13$wet_value, breaks = bins, plot = FALSE)
c1 <- rgb(173, 216, 230, max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255, 192, 203, max = 255, alpha = 80, names = "lt.pink")
plot(hg13, col = c1)
plot(hg16, col = c2, add = TRUE)

m1316_PCB16 <- subset(m1316_PCB, m1316_PCB$year == "2016")
m1316_PCB13 <- subset(m1316_PCB, m1316_PCB$year == "2013")
bins = seq(0, 80, by = 2)
hg16 <- hist(m1316_PCB16$wet_value, breaks = bins, plot = FALSE)
hg13 <- hist(m1316_PCB13$wet_value, breaks = bins, plot = FALSE)
c1 <- rgb(173, 216, 230, max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255, 192, 203, max = 255, alpha = 80, names = "lt.pink")
plot(hg13, col = c1)
plot(hg16, col = c2, add = TRUE)




