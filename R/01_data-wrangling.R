##This script performs data cleaning on the 2013-2016 mussel data.
##To be added and corrected as needed for additional data.

#### setup ####

#load libraries
library(here)
library(readxl)
library(writexl)
library(janitor)
library(stringi)
library(tidyverse)

# set directories for raw and soon to be cleaned data
raw_data_dir <- here("data", "raw")
clean_data_dir <- here("data", "clean")

#get excel file name
m1316 <- "2013_2015MusselPOPsPAHs_DRAFT_jw.xlsx"

#load dataframe
m1316_df <- here(raw_data_dir, m1316) %>% 
  read_xlsx(sheet = 2)

#### clean-up ####

#clean column names
m1316_df <- clean_names(m1316_df)

#list for substrate name changes (using named list, for easy changing if needed)
substrates <- list("cobble/mud"	=	"cobble_mud",
                   "cobble/sand"	=	"cobble_sad",
                   "Cobble-gravel mix"	=	"cobble_gravel",
                   "Cobble-gravel mix (with sand)"	=	"cobble_gravel_sand",
                   "Cobble-gravel mix (with underlying sand)"	=	"cobble_gravel_sand",
                   "Cobble-gravel mix; Sand-gravel mix"	=	"cobble_gravel_sand",
                   "Mud, silt"	=	"mud_silt",
                   "Mud, silt (cobble below)"	=	"cobble_mud_silt",
                   "sand/mud"	=	"mud_sand",
                   "Sand-cobble mix"	=	"cobble_sand",
                   "Sand-gravel mix"	=	"gravel_sand",
                   "Sand-mud and clay mix"	=	"clay_mud_sand",
                   "Sand-mud mix"	=	"mud_sand")
                   
#find and replace substrate names
m1316_df$substrate<- stri_replace_all_regex(m1316_df$substrate,
                                            pattern = c(names(substrates)),
                                            replacement = c(unname(substrates)),
                                            vectorize = FALSE)

#### write data ####

m1316_df %>% 
  write_csv(file = here(clean_data_dir, "totals_1316.csv"))
