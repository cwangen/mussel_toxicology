#### code for maps only - I don't think it runs currently ####

#### Setup ####

library(MASS)
library(tidyverse)
library(here)
library(lme4)
library(RLRsim)
library(ggdist)
library(tidyquant)
library(dplyr)
library(ggpubr)
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



# create analyte dataframes
PA_df <- mussel_df %>%
  filter(mussel_df$analyte == "SumPAHs16")
PA_df <- PA_df[complete.cases(PA_df),]






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
#mapPS <- readRDS("mapPS.rds")
 mapPS <- readRDS("./R/mapPS.rds")

#### base map of Puget Sound ####

## code for Puget Sound map was graciously provided by
## Markus Min (mmin@uw.edu, https://github.com/markusmin)

## load USA shape file
usa_spdf <- readOGR(dsn = here("R","map_files", "USA_adm0.shp"))
## convert to df(ish)
usa_spdf_fort <- tidy(usa_spdf)


## draw Puget Sound
puget_sound <- ggplot(usa_spdf_fort, aes(x = long, y = lat, group = group)) +
  geom_polygon(color = "gray70", fill = rgb(251, 234, 194, max = 255)) +
  ylab("Latitude") +
  xlab("Longitude") +
  coord_fixed(xlim = c(-123.5, -122), ylim = c(47, 49), ratio = 1.3) +
  theme(plot.background = element_rect(fill = "white"),
        panel.background = element_rect(fill="white", color = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) + 
  scale_x_continuous(breaks = c(-123, -122.5),
                     expand = c(0, 0),
                     labels=c(expression(paste(123*degree,"W")),
                              expression(paste(122.5*degree,"W")))) +
  scale_y_continuous(breaks = seq(47.5, 48.5, 0.5),
                     expand = c(0, 0),
                     labels=c(expression(paste(47.5*degree,"N")),
                              expression(paste(48*degree,"N")),
                              expression(paste(48.5*degree,"N"))))# + 
# annotate("segment", x = -122.47, xend = -122.47, y = 47.55,
#         yend = 47.67, lwd = 0.5, color = "black") + 
#annotate("segment", x = -122.3, xend = -122.3, y = 47.55,
#        yend = 47.67, lwd = 0.5, color = "black") + 
#annotate("segment", x = -122.3, xend = -122.47, y = 47.55,
#        yend = 47.55, lwd = 0.5, color = "black") + 
#annotate("segment", x = -122.3, xend = -122.47, y = 47.67,
#        yend = 47.67, lwd = 0.5, color = "black") 



#### Sort Data   ####
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

AH13 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPAHs16", musselabs_df$year == "2013")

AH16 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPAHs16", musselabs_df$year == "2016")
AH18 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPAHs16", musselabs_df$year == "2018")
AH20 <- musselabs_df %>%
  filter(musselabs_df$analyte == "SumPAHs16", musselabs_df$year == "2020")


#### PBDE Maps ####

BDE13m <- ggmap(mapPS) +
  #puget_sound +
  geom_point(
    data = BDE13,
    size = I(2.5),
    alpha = 0.4,
    aes(
      x = longitude,
      y = latitude,
      color = log(dry_value)
    ),
    inherit.aes = FALSE) +
   scale_colour_gradientn(
    colours = c("yellow", "red"),
    limits = c(0, 5),
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
#pdf("PBDEmaps.pdf")
allPBDE <- ggarrange(BDE13m, BDE16m, BDE18m, BDE20m,
                     ncol = 2, nrow = 2
)
annotate_figure(allPBDE, top = text_grob("PBDEs over time",
                                         color = "black",
                                         face = "bold",
                                         size = 14 
))
dev.off

#### PCB Maps ####

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
    data = CB18,
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
    data = CB20,
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
annotate_figure(allPCB, top = text_grob("PCBs over time",
                                        color = "black",
                                        face = "bold",
                                        size = 14
))


#### PAH Maps ####

AH13m <- ggmap(mapPS) +
  geom_point(
    data = AH13,
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
    limits = c(0, 5000)
  ) +
  ggtitle("2013")

AH16m <- ggmap(mapPS) +
  geom_point(
    data = AH16,
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
    limits = c(0, 5000)
  ) +
  ggtitle("2016")

AH18m <- ggmap(mapPS) +
  geom_point(
    data = AH18,
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
    limits = c(0, 5000)
  ) +
  ggtitle("2018")

AH20m <- ggmap(mapPS) +
  geom_point(
    data = AH20,
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
    limits = c(0, 5000)
  ) +
  ggtitle("2020")

# all plotted

allPAH <- ggarrange(AH13m, AH16m, AH18m, AH20m,
                    ncol = 2, nrow = 2
)
annotate_figure(allPAH, top = text_grob("PAHs over time",
                                        color = "black",
                                        face = "bold",
                                        size = 14 
))


#### map of all sample locations ####
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


