# Read in libraries and source code
setwd("C:/Users/Ben/Desktop/Grad School Research/Cardinals/Plumage/Averaged.Spec.readings")

source("RtcsFunctions.R")	
load("bluetitss.dat")
source("visualization.R")

library(pavo)
library(phytools)
library(dplyr)
library(rgl)

files <- dir(pattern = "*.csv")
reads <- files%>%
  purrr::map(read_csv)


color_model_df <- reads%>%
  future_map(~ gather(., code, values, 2:ncol(.))%>%
               #mutate(values = abs(values))%>%
               separate(code, c("patch", "speices_code"), sep = "_")%>%
               mutate(patch = case_when(patch == "Ventrail.tail" ~ "Ventral.tail",
                                        patch == "Undertail" ~ "Ventral.tail",
                                        patch == "PirFaries" ~ "Pirmaries",
                                        patch == "Wingbars" ~ "Wingbar",
                                        patch == "Beak.2" ~ "Back.2",
                                        patch == "Fantle" ~ "Mantle",
                                        patch == "FediuF.Coverts" ~ "Medium.Coverts",
                                        patch == "Mante.2" ~ "Mantle.2",
                                        patch == "RuFp" ~ "Rump",
                                        TRUE ~ as.character(patch)),
                      speices_code = case_when(speices_code == "GranSalF" ~ "GranSelF",
                                               speices_code == "GranSalM" ~ "GranSelM",
                                               speices_code == "CardsinM" ~ "CardSinM",
                                               TRUE ~ as.character(speices_code)))%>%
              # filter(patch %in% overlap_patches)%>%
               separate(speices_code, c("species", "sex"), sep = -1)%>%
              # filter(!species %in% bad_birds)%>%
               unite(code, c("patch", "sex"), sep = "_")%>%
               spread(code, values)%>%
               group_by(species)%>%
               nest())%>%
  reduce(bind_rows)

  (bind_row)%>%
  mutate(refstims = map(data, ~ stim(reads, ss)%>%
         refcart = map(refstims, ~ cartCoord(.x)),
         refsphere = map(refcart, ~ sphereCoord(.x)),
         rmax = map(refsphere, ~ rMax(.x)),
         acheivedR = map(c(refsphere, rmax), ~ acheivedr(.x))),
          spans = map(c(refchart, ~ colorSpan(.x))), 
reduce(bind_rows)
)



#Read in the example reflectance spectra
refs <- read.csv("Amaur_moe_F_output.csv")
refs2 <- read.csv("Amaur_moe_M_output.csv")

## This file is one that I manually concatenated for this species, so both sexes are in one file. Run this first so you can get an idea of the output. The "spans" is what you should look for
refs <- read.csv("refs.csv")

refstims <- stim(refs,ss)

#Convert the stimulation values to Cartesian Coordinates (for later analyses)
refcart <- cartCoord(refstims)

#Converte the Cartesian Coordinates to Sphereical Coordinates (for later analyses)
refsphere <- sphereCoord(refcart)

#Calculate the maximum possible r-value for each hue (angles from the origin)
rmax <- rMax(refsphere)

#Calculate the acheived chroma for each patch (r/rmax)
acheivedr <- acheivedR(refsphere,rmax)

#Normalized brilliance for each patch (aka. brightness)
normbrill <- normBrill(refs)

#Color volume for all patches submitted (minimum convex polygon of all points in the tetracolorspace)
vol <- colorVolume(refcart)

#Hue disparity matrix (all patches compared to each other)
disp <- hueDisp(refsphere)

#Average, variance and maximum hue disparity
disp.summary <- summary.hueDisp(disp)

#Color span matrix (all patches compared to each other)----> This is the major one that I need! The is calculated the Euclidean distance for each particular patch within tetrahedralcolorspace. You can think of it as the color of each patch represents a 3-Dimensional corrdinate within this multi-demensional space 
spans <- colorSpan(refcart)

#Average, variance and maximum color span
spans.summary <- summary.colorSpan(spans)

#Average chroma
avgChroma(refsphere)

#Average acheived chroma
avgAcheivedChroma(acheivedr)

#Average brilliance
avgBrill(normbrill)

#A number of summary measurements from all patches submitted (Average color span, variance in color span, maximum color span, color volume, average hue disparity, variance in hue disparity, maximum hue disparity, average brilliance, average chroma, and average acheived chroma)
summ <- summaryOfpatches(refs,ss)



###### Attempt two; combing files by like matching names and 

####  Remove brillance values less than 0.05 for all individual plumage comparisons 
#Color span matrix (all patches compared to each other)----> This is the major one that I need! The is calculated the Euclidean distance for each particular patch within tetrahedralcolorspace. You can think of it as the color of each patch represents a 3-Dimensional corrdinate within this multi-demensional space 

refstims <- stim(refs,ss)

#Convert the stimulation values to Cartesian Coordinates (for later analyses)
refcart <- cartCoord(refstims)

#Converte the Cartesian Coordinates to Sphereical Coordinates (for later analyses)
refsphere <- sphereCoord(refcart)

#Calculate the maximum possible r-value for each hue (angles from the origin)
rmax <- rMax(refsphere)

#Calculate the acheived chroma for each patch (r/rmax)
acheivedr <- acheivedR(refsphere,rmax)

#Normalized brilliance for each patch (aka. brightness)
normbrill <- normBrill(refs)

#Color volume for all patches submitted (minimum convex polygon of all points in the tetracolorspace)
vol <- colorVolume(refcart)

#Hue disparity matrix (all patches compared to each other)
disp <- hueDisp(refsphere)

#Average, variance and maximum hue disparity
disp.summary <- summary.hueDisp(disp)

#Color span matrix (all patches compared to each other)----> This is the major one that I need! The is calculated the Euclidean distance for each particular patch within tetrahedralcolorspace. You can think of it as the color of each patch represents a 3-Dimensional corrdinate within this multi-demensional space 
spans <- colorSpan(refcart)

#Average, variance and maximum color span
spans.summary <- summary.colorSpan(spans)

#Average chroma
avgChroma(refsphere)

#Average acheived chroma
avgAcheivedChroma(acheivedr)

#Average brilliance
avgBrill(normbrill)

#A number of summary measurements from all patches submitted (Average color span, variance in color span, maximum color span, color volume, average hue disparity, variance in hue disparity, maximum hue disparity, average brilliance, average chroma, and average acheived chroma)
summ <- summaryOfpatches(refs,ss)


