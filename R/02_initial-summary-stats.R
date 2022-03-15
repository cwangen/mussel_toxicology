## This script performs data summaries on the 2013-2020 mussel data.


#### setup ####

# load libraries
library(here)
library(tidyverse)
library(ggplot2)
library(ggmap)
library(sf)
library(mapview)
library(RANN)
library(geosphere)
library(stringi)
library(ggpubr)

# set directories
clean_data_dir <- here("data", "clean")

# get excel file name
totals <- "totals_all.csv"

# load dataframe
mussel_df <- read_csv(here(clean_data_dir, totals), show_col_types = FALSE)

#### reorganization as needed ####

# remove rows with missing wet_value (there is only one)
mussel_exp <- mussel_df[!is.na(mussel_df$wet_value), ]

#### tables and numbers ####

# number of sites per year
site_totals <- table(mussel_exp$analyte, mussel_exp$year)
print(site_totals)

# analyte min/max/mean table
by_analyte <- mussel_exp %>% group_by(analyte)
analyte_table <- by_analyte %>% summarise(
  avg_wv = mean(wet_value),
  min_wv = min(wet_value),
  max_wv = max(wet_value)
)


#### nearest neighbors numbers ####

# remove test sites (known repeats, currently not used) and choose one analyte
basenn_df <- mussel_df
# basenn_df <- basenn_df[!is.na(mussel_df$fund_source),]
# basenn_df <- subset(basenn_df, site_name!="Penn Cove Reference")
basenn_df <- basenn_df %>% filter(basenn_df$analyte == "lipids")


# number of unique site names (226, 111 have duplicates)
uniqnames <- unique(basenn_df$site_name)
uniqnames <- length(uniqnames)

# count number of neighbors in 400m (within a single year)
nn_df <- basenn_df
nn_df <- nn_df[order(nn_df$ID), ]
latlon <- cbind(nn_df$ID, nn_df$latitude, nn_df$longitude, nn_df$year)

latlon13 <- latlon[latlon[, 4] == 2013, ]
latlon16 <- latlon[latlon[, 4] == 2016, ]
latlon18 <- latlon[latlon[, 4] == 2018, ]
latlon20 <- latlon[latlon[, 4] == 2020, ]
nn_13 <- cbind(latlon13, X = rowSums(distm(latlon13[, 3:2],
  fun = distHaversine
) <= 400)) # 400m
nn_16 <- cbind(latlon16, X = rowSums(distm(latlon16[, 3:2],
  fun = distHaversine
) <= 400)) # 400m
nn_18 <- cbind(latlon18, X = rowSums(distm(latlon18[, 3:2],
  fun = distHaversine
) <= 400)) # 400m
nn_20 <- cbind(latlon20, X = rowSums(distm(latlon20[, 3:2],
  fun = distHaversine
) <= 400)) # 400m

nn_all <- rbind(nn_13, nn_16, nn_18, nn_20)
nn_all <- data.frame(nn_all)

names(nn_all) <- c("index", "latitude", "longitude", "year", "n_neighbors")
nn_all <- nn_all[order(nn_all$index), ]
nn_all$n_neighbors <- nn_all$n_neighbors - 1 # minus 1 b/c includes self
# number_n <- sum(nn_df$n_neighbors>0) #sites have neighbor
nnshort_df <- basenn_df[order(basenn_df$ID), ]
nnshort_df$n_neighbors <- round(nn_all$n_neighbors) # df with nn

# count number of neighbors in 5000m (within a single year)
nn_df <- basenn_df
nn_df <- nn_df[order(nn_df$ID), ]
latlon <- cbind(nn_df$ID, nn_df$latitude, nn_df$longitude, nn_df$year)

latlon13 <- latlon[latlon[, 4] == 2013, ]
latlon16 <- latlon[latlon[, 4] == 2016, ]
latlon18 <- latlon[latlon[, 4] == 2018, ]
latlon20 <- latlon[latlon[, 4] == 2020, ]
nn_13 <- cbind(latlon13, X = rowSums(distm(latlon13[, 3:2],
  fun = distHaversine
) <= 5000)) # 5000m
nn_16 <- cbind(latlon16, X = rowSums(distm(latlon16[, 3:2],
  fun = distHaversine
) <= 5000)) # 5000m
nn_18 <- cbind(latlon18, X = rowSums(distm(latlon18[, 3:2],
  fun = distHaversine
) <= 5000)) # 5000m
nn_20 <- cbind(latlon20, X = rowSums(distm(latlon20[, 3:2],
  fun = distHaversine
) <= 5000)) # 5000m

nn_all <- rbind(nn_13, nn_16, nn_18, nn_20)
nn_all <- data.frame(nn_all)

names(nn_all) <- c("index", "latitude", "longitude", "year", "n_neighbors")
nn_all <- nn_all[order(nn_all$index), ]
nn_all$n_neighbors <- nn_all$n_neighbors - 1 # minus 1 b/c includes self
nnlong_df <- basenn_df[order(basenn_df$ID), ]
nnlong_df$n_neighbors <- round(nn_all$n_neighbors) # dataframe with number of neighbors


# 400m bar
shorttable <- table(nnshort_df["n_neighbors"])
barshort <- barplot(shorttable,
  main = "Neighbors in 400m (within year)",
  xlab = "Number of Neighbors",
  ylab = "Count",
  ylim = c(0, 350)
)

# 5000m bar
longtable <- table(nnlong_df["n_neighbors"])
barlong <- barplot(longtable,
  main = "Neighbors in 5000m (within year)",
  xlab = "Number of Neighbors",
  ylab = "Count",
  ylim = c(0, 125)
)



#### maps ####


# PLEASE NOTE: This will not run without a Google API but
# But is the code that makes the Puget Sound map loaded below
# mapPS <- get_map(location= c(-124.793701,46.420957,-121.827393,48.999396),
#                zoom=10,
#                maptype = "roadmap",
#                source='google',
#                color='color',
#                force = TRUE)

# Load Puget Sound map
mapPS <- readRDS("./R/mapPS.rds")


## Sort Data
musselabs_df <- mussel_df
musselabs_df$dry_value <- abs(musselabs_df$dry_value)

BDE13 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPBDEs11", musselabs_df$year == "2013")

BDE16 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPBDEs11", musselabs_df$year == "2016")
BDE18 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPBDEs11", musselabs_df$year == "2018")
BDE20 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPBDEs11", musselabs_df$year == "2020")

CB13 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPCBs2x17", musselabs_df$year == "2013")
CB16 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPCBs2x17", musselabs_df$year == "2016")
CB18 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPCBs2x17", musselabs_df$year == "2018")
CB20 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPCBs2x17", musselabs_df$year == "2020")

## PBDE Maps

BDE13m <- ggmap(mapPS) +
  geom_point(
    data = BDE13,
    size = I(1.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = dry_value
    )
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 50)
  ) +
  ggtitle("2013")

BDE16m <- ggmap(mapPS) +
  geom_point(
    data = BDE16,
    size = I(1.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = dry_value
    )
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 50)
  ) +
  ggtitle("2016")

BDE18m <- ggmap(mapPS) +
  geom_point(
    data = BDE18,
    size = I(1.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = dry_value
    )
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 50)
  ) +
  ggtitle("2018")

BDE20m <- ggmap(mapPS) +
  geom_point(
    data = BDE20,
    size = I(1.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = dry_value
    )
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 50)
  ) +
  ggtitle("2020")

# all plotted

allPBDE <- ggarrange(BDE13m, BDE16m, BDE18m, BDE20m,
  ncol = 2, nrow = 2
)
annotate_figure(allPBDE, top = text_grob("PBDEs over time",
  color = "black",
  face = "bold",
  size = 14 
))

## PCB Maps

CB13m <- ggmap(mapPS) +
  geom_point(
    data = CB13,
    size = I(1.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = dry_value
    )
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 285)
  ) +
  ggtitle("2013")

CB16m <- ggmap(mapPS) +
  geom_point(
    data = CB16,
    size = I(1.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = dry_value
    )
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 285)
  ) +
  ggtitle("2016")

CB18m <- ggmap(mapPS) +
  geom_point(
    data = CB13,
    size = I(1.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = dry_value
    )
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 285)
  ) +
  ggtitle("2018")

CB20m <- ggmap(mapPS) +
  geom_point(
    data = CB13,
    size = I(1.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = dry_value
    )
  ) +
  scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 285)
  ) +
  ggtitle("2020")

# all plotted

allPCB <- ggarrange(CB13m, CB16m, CB18m, CB20m,
  ncol = 2, nrow = 2
)
annotate_figure(allPBDE, top = text_grob("PCBs over time",
  color = "black",
  face = "bold",
  size = 14
))



# map of all sample locations
mussel_l <- mussel_df %>% filter(mussel_df$analyte == "lipids")

lon <- mussel_l$longitude
lat <- mussel_l$latitude
df <- as.data.frame(cbind(lon, lat))
allMap <- qmplot(lon,
  lat,
  data = df,
  maptype = "toner-lite",
  color = I("red"),
  size = I(0.2),
  main = "All Sampled Locations"
)

#### histograms of wet_values for different analytes ####

# get wet values
mussel_PDBE <- subset(
  mussel_exp,
  mussel_exp$analyte == "SumPBDEs11"
)
mussel_PCB <- subset(
  mussel_exp,
  mussel_exp$analyte == "SumPCBs2x17"
)

# histograms combined years
hist(mussel_PDBE$dry_value)
hist(mussel_PCB$dry_value)