##This script performs data cleaning on the 2013-2020 mussel data.

#### setup ####

#load libraries
library(here)
library(readxl)
library(writexl)
library(janitor)
library(stringi)
library(tidyverse)
library(dplyr)

#set directories for raw and soon to be cleaned data
raw_data_dir <- here("data", "raw")
clean_data_dir <- here("data", "clean")

#get excel file name
mussel <- "2013-20MusselCagesPOPsPAHs_Cnty_WRIA_LIO_Coverages.xlsx"

##load excel sheet
mussel_df <- here(raw_data_dir, mussel) %>%
  read_xlsx(sheet = 1)

#lipid totals
lipids_df <- here(raw_data_dir, mussel) %>%
  read_xlsx(sheet = 1)
lipids_df <- subset(lipids_df, Analyte == "lipids")
#PBDEs totals
PBDE_df <- here(raw_data_dir, mussel) %>%
  read_xlsx(sheet = 1)
PBDE_df <- subset(PBDE_df, Analyte == "SumPBDEs11")
#PCB totals
PCB_df <- here(raw_data_dir, mussel) %>%
  read_xlsx(sheet = 1)
PCB_df <- subset(PCB_df, Analyte == "SumPCBs2x17")
#PAH totals
PAH_df <- here(raw_data_dir, mussel) %>%
  read_xlsx(sheet = 1)
PAH_df <- subset(PAH_df, Analyte == "SumPAHs16")
#combine to one dataframe
mussel_df <- rbind(lipids_df, PCB_df,PBDE_df,PAH_df)


#### clean-up ####

#absolute value of wet_value for now

#clean column names
mussel_df <- clean_names(mussel_df)

#get parentheses out of substrates to avoid escape
mussel_df$substrate <- gsub("[()]", "", mussel_df$substrate)

#list for substrate name changes (using named list, for easy changing if needed)
substrates <- list("Cobble-gravel mix with sand"	=	"cobble_gravel_sand",
                   "Cobble-gravel mix with underlying sand"	=	"cobble_gravel_sand",
                   "Mud, silt cobble below"	=	"cobble_mud_silt",
                   "Bedrock, hardpan" = "bedrock_hardpan",
                   "Cobble-gravel mix; Sand-gravel mix"	=	"cobble_gravel_sand",
                   "Cobble-gravel mix, Sand-gravel mix"	=	"cobble_gravel_sand",
                   "cobble/mud"	=	"cobble_mud",
                   "cobble/sand"	=	"cobble_sad",
                   "N/A" = "NA",
                   "Riprap rock" = "riprap",
                   "Sand-gravel mix, sand, san-mud mix" = "gravel_mud_sand",
                   "Sand-mud and clay mix"	=	"clay_mud_sand",
                   "Sand-cobble mix"	=	"cobble_sand",
                   "Sand-gravel mix"	=	"gravel_sand",
                   "Sand-mud mix"	=	"mud_sand",
                   "sand/mud"	=	"mud_sand",
                   "Sand" = "sand",
                   "Cobble-gravel mix"	=	"cobble_gravel",
                   "Mud, silt"	=	"mud_silt"
                   )
                   
#find and replace substrate names
mussel_df$substrate <- stri_replace_all_regex(mussel_df$substrate,
                                            pattern = c(names(substrates)),
                                            replacement = c(unname(substrates)),
                                            vectorize = FALSE)

#list for site_name fixes (using named list, for easy changing if needed)
sitenames <- list("Bellingham Bay, Little Squalicum Crk" = "Bellingham Bay, Little Squalicum Creek",
                  "Cavalero Beach Co. Park" = "Calvalero Beach",
                  "Cherry Point North" = "Cherry Point",
                  "Cherry Pt Aq Res, 1 Alcoa-BP" = "Cherry Point",
                  "Gig Harbor - Boat Launch" = "Gig Harbor Boat Launch",
                  "Illahee Crk" = "Illahee Creek",
                  "Manchester, Stormwater Outfall" = "Manchester, SWO",
                  "Shilshole" = "Shilshole Bay"
)

#Fixing site_name
mussel_df$site_name <- stri_replace_all_regex(mussel_df$site_name,
                                              pattern = c(names(sitenames)),
                                              replacement = c(unname(sitenames)),
                                              vectorize = FALSE)

 
#insert Penn Cove lat and lon (later spreadsheets include this already)
#mussel_df$latitude[is.na(mussel_df$latitude)] <-  48.218626
#mussel_df$longitude[is.na(mussel_df$longitude)] <- -122.707972

#### create additional columns ####

#add numeric ID column for each sample
mussel_df <- mussel_df %>%
  group_by(year, sample_id) %>%
  dplyr::mutate(ID = cur_group_id())

#add lipid_weight (wet weight/lipids) column
mussel_df["lipid_weight"] <- NA
mussel_df <- mussel_df %>%
  group_by(ID) %>%
  mutate(lipid_weight = wet_value / wet_value[analyte == "lipids"])

#abs of dry value for now
mussel_df$dry_value <- abs(mussel_df$dry_value)

#add column with days sample spent in field
mussel_df$time <- difftime(mussel_df$retrieval_date, mussel_df$deployment_date)

#label Penn Cove Reference Samples in LIO/WRIA Column

#turn WRIA to character

mussel_df$wria_nr <- as.character(mussel_df$wria_nr)


PennCove <- c("Penn Cove Reference")
a = PennCove[1]

mussel_df[mussel_df$site_name == a, which( colnames(mussel_df)=="lio_areas" )] <- "Penn Cove Reference"


mussel_df[mussel_df$site_name == a, which( colnames(mussel_df)=="wria_nr" )] <- "Penn Cove Reference"


PennCove <- c("Penn Cove, Pre-test Baseline",
              "Penn Cove, Pre-test Baseline 1",
              "Penn Cove, Pre-test Baseline 2",
              "Penn Cove, Pre-test Baseline 3",
              "Penn Cove, Pre-test Baseline 4",
              "Penn Cove, Pre-test Baseline 5",
              "Penn Cove, Pre-test Baseline 6",
              "Penn Cove Baseline #1",
              "Penn Cove Baseline #2",
              "Penn Cove Baseline #3")

for (i in 1:length(PennCove))
{
  a = PennCove[i]
  mussel_df[mussel_df$site_name == a, which( colnames(mussel_df)=="lio_areas" )] <- "Penn Cove Baseline"
}

for (i in 1:length(PennCove))
{
  a = PennCove[i]
  mussel_df[mussel_df$site_name == a, which( colnames(mussel_df)=="wria_nr" )] <- "Penn Cove Baseline"
}


#### write data ####
mussel_df %>%
  write_csv(file = here(clean_data_dir, "totals_all.csv"))