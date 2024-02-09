###############################################################################
## LeafArea
###############################################################################
# video to help - https://vimeo.com/471493529d

## Installing packages
library("LeafArea")
library("tidyverse")


# error "no images in the directory"
# 1 image not procssed
  data_10200 <- run.ij(set.directory = "C:/Users/monik/Pictures/sla_2024/10200",
                     distance.pixel = 1182,
                     known.distance = 10)

# error "no images in the directory"
# 45 images not processed
data_3400 <- run.ij(set.directory = "C:/Users/monik/Pictures/sla_2024/3400",
                    distance.pixel = 1569,
                    known.distance = 10,
                    low.size = 0.1)

# mixed dimensions that were all similar ~2000, pixel at 10cm = 1182
# worked, buy only ran 270 images instead of the full 309.
# 39 images not processed
data_2000s <- run.ij(set.directory = "C:/Users/monik/Pictures/sla_2024/2000s",
                    distance.pixel = 1182,
                    known.distance = 10,
                    low.size = 0.1,
                    set.memory = 5)

## merging all the data
data_merge <- rbind(data_10200, data_2000s, data_3400)


## writing the data
write.csv(data_merge, file = "data_leaf.area.csv")




###############################################################################
## Final notes
###############################################################################
## worked thanks to EVAN AND KELLY
## thank you