getwd()
setwd("C:/Users/monik/Documents/reserach/fieldwork_summer_2023/
      sp_n&cvalues/sp._for_try/")
getwd()

trysp <- read.delim2(file = "TryAccSpecies.txt", header = TRUE)

tryspdf <- as.data.frame(trysp)

tryspdf <- tolower(tryspdf)

names(trysp) <- tolower(names(trysp))


trysp_lower <- tolower (trysp$AccSpeciesName)
head(trysp_lower)


trysp$AccSpeciesName <- tolower(trysp$AccSpeciesName)
head(trysp)

trysp <- tolower(trysp)
head(trysp)

colnames(trysp)

write.csv(trysp, file = "TryAccspecies.csv", header = TRUE)

write.csv(trysp, file = "TryAccSpecies_mkedits.txt")

mksp <- read.csv(file = "species_list_mk.csv", header = TRUE)

head(trysp)

head(mksp)

mkspecies_only <- mksp[ , 3]
head(mkspecies_only)

try.sp.numb <- subset(trysp, AccSpeciesName == mkspecies_only)


length(mkspecies_only)
length(trysp$AccSpeciesName)

## Did this look for all the species in mk's list in the try dataset? 
# I think it worked. 
# https://sparkbyexamples.com/r-programming/r-select-rows-based-on-column-value/ 
test <- subset(trysp, AccSpeciesName %in% mkspecies_only)


