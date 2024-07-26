################################################################################
## Purpose: Process all the leaf images for leaf area
## Author: Monika Kelley
## Date: 2024/02/09
################################################################################

## Notes -----------------------------------------------------------------------
## helpful videos and documents
# video - https://vimeo.com/471493529
# video - https://www.youtube.com/watch?v=CU7gG3_OF34&t=309s  
# doc - https://imagej.net/ij/docs/pdfs/examples.pdf 
# doc - https://www.petiolepro.com/blog/how-to-measure-leaf-area-with-imagej/


## libraries -------------------------------------------------------------------
library("LeafArea")
library("tidyverse")

################################################################################
## Rapid and bulk leaf area analysis
################################################################################

## Image dimensions ~2000
data_2000s <- run.ij(set.directory = "C:/Users/monik/Pictures/leaf_area/02_scans_cropped_rcode/2000s/",
                     distance.pixel = 1182,
                     known.distance = 10,
                     low.size = 0.1,
                     set.memory = 5) 


## Image dimensions 10200
data_10200 <- run.ij(set.directory = "C:/Users/monik/Pictures/leaf_area/02_scans_cropped_rcode/10200",
                     distance.pixel = 1182,
                     known.distance = 10)


## Image dimensions 3400
data_3400 <- run.ij(set.directory = "C:/Users/monik/Pictures/leaf_area/02_scans_cropped_rcode/3400/",
                    distance.pixel = 1569,
                    known.distance = 10,
                    low.size = 0.1)


### write csv for qc/qa  -------------------------------------------------------
## merging all the data
leaf_area_auto <- rbind(data_10200, data_2000s, data_3400)


## writing the data
write.csv(leaf_area_auto, file = "data/03_rproducts/leaf_area_auto.csv")
