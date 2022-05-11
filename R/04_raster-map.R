#### Load libraries  ####

library(here)
library(rgdal)
library(broom)
library(ggplot2)
library(cowplot)
library(patchwork)
library(raster)
library(sp)
library(rgdal)


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
  #  theme(plot.background = element_rect(fill = "white"),
  #       panel.background = element_rect(fill="white", color = "black"),
  #      panel.border = element_rect(colour = "black", fill=NA, size=1),
  #     panel.grid.major = element_blank(),
  #    panel.grid.minor = element_blank()) + 
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

