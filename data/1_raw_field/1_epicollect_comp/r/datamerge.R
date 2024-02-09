library("readxl")
setwd("C:/Users/monik/Documents/C_Reserach/fieldwork_summer_2023/data/
      epicollect_comp/20230918_all sites/combined_ds")

data <- read_excel("sp.comp_mksummer_01.xlsx")
setwd("C:/Users/monik/Documents/C_Reserach/fieldwork_summer_2023/data/
      epicollect_comp/20230918_all sites/")
data

install.packages("readr")
library(readr)

metadata <- read_csv("form-1__mk-nutnet2023-5dominantplants.csv")
head(metadata, 5)


arch.us <- read_csv("branch-3__species-composition-arch.csv")
temple.us <- read_csv("branch-1__species-composition-temple.csv")
sevi.us <- read_csv("branch-4__species-composition-sevi.csv")
lubb.us <- read_csv("branch-2__species-composition-lubb.csv")

typeof(arch.us)
typeof(temple.us)

df_metadata <- data.frame(metadata)
df_arch <- data.frame(arch.us)
df_temple <- data.frame(temple.us)
df_sevi <- data.frame(sevi.us)
df_lubb <- data.frame(lubb.us)

mergedfile <- merge(df_metadata, df_arch, df_lubb, df_sevi, df_temple
                    by.x = "epiid", by.y = "id", all.x = TRUE)



##### Attempt indiviudal merges
merged_df <- merge(df_metadata, df_arch, by.x = "ec5_uuid", 
                   by.y = "ec5_branch_owner_uuid", all.x = TRUE)

merged_df <- merge (merged_df, df_lubb, by.x = "ec5_uuid", 
                   by.y = "ec5_branch_owner_uuid", all.x = TRUE)

merged_df <- merge (merged_df, df_sevi, by.x = "ec5_uuid", 
                   by.y = "ec5_branch_owner_uuid", all.x = TRUE)


write.csv(merged_df)

csv_file_path <- "C:/Users/monik/Documents/C_Reserach/
fieldwork_summer_2023/data/epicollect_comp/r"
write.csv(merged_df, file = csv_file_path, row.names = FALSE)


csv_file_path <- "C:/Users/monik/Documents/C_Reserach/
fieldwork_summer_2023/data/epicollect_comp"

write.csv(merged_df, file = csv_file_path, row.names = FALSE)


write.csv(merged_df, file = 'C:/Users/monik/Documents/C_Reserach
          /fieldwork_summer_2023/data/epicollect_comp/', row.names = TRUE)




file_path <- 'C:/Users/monik/Documents/'

write.csv(merged_df, file = file.path(file_path, 'output.csv'), 
          row.names = TRUE)



