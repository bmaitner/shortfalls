#Code taken from https://github.com/pebgroup/Distribution_of_phylogenetic_data/blob/master/03_phylo_know_analysis_2022.R and variable locations changed

library(data.table)
library(rgdal)
library(tidyr)

wcsp.names <- fread("manual_downloads/WCVP/wcvp_names.txt", quote = "") # WCVP names file
wcsp.distribution <- fread("manual_downloads/WCVP/wcvp_distribution.txt", quote = "") # WCVP distribution file

# Exclude fern and moss families
ferns <- fread("manual_downloads/Darwinian_shortfalls/fern_list.txt", quote = "", header = F, col.names = "family") # list of ferns
moss <- fread("manual_downloads/Darwinian_shortfalls/bryophyta_list.txt", quote = "", header = F, col.names = "family") # list of mosses
wcsp.names <- subset(wcsp.names, !family %in% ferns$family)
wcsp.names <- subset(wcsp.names, !family %in% moss$family)

# Include only accepted entries
wcsp.names.species <- wcsp.names[grepl("Accepted", wcsp.names$taxon_status),]

# Synonymize plant name IDs in distribution data to use accepted plant name IDs, and subset to taxa with accepted plant name IDs

# Merge with accepted names:
wcsp.distribution <- merge(wcsp.distribution, wcsp.names.species[,c("accepted_plant_name_id", "plant_name_id", "taxon_rank")],
                           by="plant_name_id", all.x=TRUE)

# Remove all entries with no accepted plant name ID
wcsp.dist.subset <- wcsp.distribution[!is.na(wcsp.distribution$accepted_plant_name_id),]

# Remove subsp., var., f. etc that was introduced to the data after synonymizing from distribution data
wcsp.dist.subset <- wcsp.dist.subset[wcsp.dist.subset$taxon_rank == "Species", ]

# Remove extinct and introduced species
wcsp.dist.subset <- subset(wcsp.dist.subset, wcsp.dist.subset$extinct == 0)
wcsp.dist.subset <- subset(wcsp.dist.subset, wcsp.dist.subset$introduced == 0)

# Remove or correct mistaken Area codes
shape <- readOGR("level3.shp")
shape <- readOGR("manual_downloads/Darwinian_shortfalls/level3.shp")
wcsp.dist.subset$area_code_l3 = toupper(wcsp.dist.subset$area_code_l3)
wcsp.dist.subset <- wcsp.dist.subset[wcsp.dist.subset$area_code_l3 %in% shape$LEVEL_3_CO,]

# produce table of WCSP ID's and corresponding names
name.id <- wcsp.names[grepl("Accepted", wcsp.names$taxon_status),]
name.id <- name.id[!(name.id$species == ""), ]
name.id <- subset(name.id, select = c(plant_name_id, taxon_name))

# Replace WCSP ID's with accepted names using the ID and Names object
wcsp.dist.subset$accepted_plant_name <- name.id$taxon_name[match(wcsp.dist.subset$accepted_plant_name_id, name.id$plant_name_id)]
wcsp.dist.subset <- na.omit(wcsp.dist.subset)

# Replace "Ã-" with "x" to indicate accepted hybrids, as this is the format used in GenBank
wcsp.dist.subset <- as.data.frame(sapply(wcsp.dist.subset, gsub, pattern = "Ã-", replacement = "x"))




# ---- Format the WCSP subset into a dataframe that includes L3 areas as columns ---- #




# Create vector of unique L3 area codes to use in the loop
l3.vector <- c(as.character(unique(wcsp.dist.subset$area_code_l3)))

# Create an empty list to fill
area.list <- list()

# Loop and search for all species associated with an l3 region and arrange them in a list
for (i in l3.vector){
  inventory <- wcsp.dist.subset$accepted_plant_name[wcsp.dist.subset$area_code_l3 == i]
  area.list[[i]] <- inventory
}

# Convert the list into the final dataframe format, which will be used in the "diversity_variables" script, where the response variable and biodiversity-related explanatory variables will be produced.
list.equal.length <- lapply(area.list, 'length<-', max(lengths(area.list)))
wcsp.df <- as.data.frame(list.equal.length) # The table is written in the end of the script




# ---- Create a table of synonyms and their accepted names to be used to translate synonyms from Genbank ---- #




# Retrive records for all not-accepted names and remove entries at genus level
wcsp.synonyms <- wcsp.names[!grepl("Accepted", wcsp.names$taxon_status),]
wcsp.synonyms <- wcsp.synonyms[!(wcsp.synonyms$species == ""), ]

# Combine columns: Genus, Genus hybrid marker, Species hybrid marker, Species, Infraspecific rank and infraspecific epithet to producea new column with the full synonym
wcsp.synonyms <- replace(wcsp.synonyms, wcsp.synonyms == "", NA)

# Remove irrelvant columns
wcsp.synonyms <- subset(wcsp.synonyms, select = c(taxon_name,accepted_plant_name_id))

# Substitute plant ID with plant name from "Names and Id's"

wcsp.synonyms$accepted_plant_name_id <- name.id$taxon_name[match(wcsp.synonyms$accepted_plant_name_id, name.id$plant_name_id)]
wcsp.synonyms <- na.omit(wcsp.synonyms)

# Remove cases where synonym = accepted name to reduce computation time during translation
# Note: These cases occur when the only difference between synonym and accepted name is the "author" columns. This is never included in the name on GenBank, so it is safe to remove these)
wcsp.synonyms <- wcsp.synonyms[wcsp.synonyms$taxon_name != wcsp.synonyms$accepted_plant_name_id,]

# Remove duplicated synonyms
# Note: Again the "author" columns is the culprit. These always refer to the same accepted name, so one is expandable to reduce computation time.
wcsp.synonyms <- wcsp.synonyms[!duplicated(wcsp.synonyms$taxon_name),]

# Remove cases where a synonym is the accepted name of another species.
# Note: This one is complicated. Several species names exists as both synonym and accepted name, but for different species. This is, again, due to differences in the author extension to the synonym, which is not included here, as it does not occur on Genbank. -
# - Including these cases would result in the possible translation of an accepted name into another accepted name, which is undesirable
wcsp.synonyms <- wcsp.synonyms[!(wcsp.synonyms$taxon_name %in% intersect(wcsp.synonyms$accepted_plant_name_id, wcsp.synonyms$taxon_name)),]

# Translate "Ã-" to "x" to adapt to the GenBank download
wcsp.synonyms <- as.data.frame(sapply(wcsp.synonyms, gsub, pattern = "Ã-", replacement = "x"))

# Remove cases where the accepted name is "Unplaced Unplaced"
wcsp.synonyms <- wcsp.synonyms[wcsp.synonyms$accepted_plant_name_id != "Unplaced Unplaced",]


# The final synonym translator file can be written in the end of the script




# ---- Write tables for further use ---- #




# Write a table that includes all accepted species and their ID's as metadata (ID's and Names)
write.table(name.id,"data\\ID_and_Names_2022.csv", sep = ",", row.names = F, quote = F)

# Write a table with the relevant subset of the WCSP distribution data as useful metadata (WCSP long format)
write.table(wcsp.dist.subset,"data\\wcsp_long_2022.csv", sep = ",", row.names = F, col.names = T, quote = F)

# Write a table with the relevant subset of the WCSP distribution data in the format necessary for further analysis (WCSP wide format)
write.table(wcsp.df,"data\\wcsp_wide_2022.csv", sep = ",", row.names = F, col.names = T, quote = F)

# Write a table with synonyms and their accepted names to be used to translate the GenBank download into only accepted names (Synonym translator)
write.table(wcsp.synonyms,"data\\synonym_translator_2022.csv", sep = ",", row.names = F, col.names = T, quote = F)




# ---- Additional information of interest included in the paper ---- #



# The total number of species with distribution data according to WCVP
length(unique(wcsp.dist.subset$accepted_plant_name))

# The total number of species in total according to WCVP
wcsp.total <- wcsp.names.species[wcsp.names.species$taxon_rank == "Species", ]
length(unique(wcsp.total$accepted_plant_name_id))

# % of total number of species with distribution data available
length(unique(wcsp.dist.subset$accepted_plant_name))/length(unique(wcsp.total$accepted_plant_name_id))

# Biodiversity variable creation script
#
# This script uses three datasets produced in the "wcvp_subset_2021" script, two subsets of the GenBank download produced in "genbank_download.R, and BIEN data:
# WCVP:
# 1) ID's and Names
# 2) WCVP distribution data in wide format with countries as columns
# 3) WCVP distribution data in a long format with two columns: species and botanical country
# GenBank:
# 4) a list of species that have at least ONE relevant sequence available
# 5) a list of species that preserves duplicate species if they appear with different molecular markers
# BIEN:
# 6) a list of species with BIEN occurrences for each botanical country, translated into WCSP ID's (BIEN_in_WCSP_regions_Sept21_2020download.RData)
# Botanical countries:
# 7) a table of Level 3 botanical counrty codes and their area size
#
# This script produces multiple explanatory variables related to biodiversity and four response variables at the scale of botanical contries, which is all gathered in one table. These variables include:
# 1) Total number of species with relevent molecular data on GenBank (response)
# 2) Inventory completeness of flora with at least one entry with relevant sequence data on GenBank (response)
# 3) Inventory completeness of flora including all relevant molecular data on GenBank (response)
# 4) Inventory completeness of flora with a record occurrence in BIEN (response)
# 5) Species richness
# 6) Mean number of countries a species of a country is present in
# 7) Median number of countries a species of a country is present in
# 8) Mean range size of species
# 9) Total number of endemic species
# 10) Endemics relative to species richness
# 11) Proporption of endemics sequenced

# Total runtime: ~5 minutes


# ---- Read data from Genbank, WCVP and BIEN ---- #


library(data.table)
gen.data <- fread("manual_downloads/Darwinian_shortfalls/genbank_entries_2022.csv", header = T, quote = "", sep = NULL) # As produced in the genbank_download.R script
load("manual_downloads/Darwinian_shortfalls/BIEN_in_WCSP_regions_sep2021.RData") # R workspace with objects "spec.list" and "res". Spec.list is a list of all L3 regions with the species recorded there in BIEN, translated to WCSP ID's. Res is a df with the lengths of the elements of spec.list.
L3.area <- fread("manual_downloads/Darwinian_shortfalls/L3_and_area.csv", header = T) # Botanical countries and their area size
gen.data.duplicate <- fread("manual_downloads/Darwinian_shortfalls/genbank_entries_w_duplicates_2022.csv", quote = "", header = T, sep = NULL) # As produced in the genbank_download.R script


wcsp.data <- fread("data/wcsp_wide_2022.csv", header = T) # As produced in the wcvp_subset_2021.R script
name.id <- fread("data\\ID_and_Names_2022.csv", header = T) # As produced in the wcvp_subset_2021.R script
wcsp.long <- fread("data\\wcsp_long_2022.csv", quote = "") # As produced in the wcvp_subset_2021.R script




# ---- Prepare vectors and dataframes for the variable creation loops ---- #




# Create vector of area codes and taxa with genetic data
area.vector <- c(colnames(wcsp.data))
genbank.vector <- c(as.character(gen.data$species))
genbank.vector.duplicate <- c(as.character(gen.data.duplicate$species))

# Create dataframe with data for distribution range
dist.range <- table(unlist(wcsp.data))
dist.range.df <- data.frame(dist.range)
dist.range.df.endemics <- dist.range.df[dist.range.df$Freq == 1,]

# Create dataframe of species and area, aggregate and sum it, and create a wide-style dataframe with sqkm replacing species.
wcsp.long$area_sq <- L3.area$AREA_SQKM[match(wcsp.long$area_code_l3, L3.area$LEVEL_3_CO)]
L3.area.only <- subset(wcsp.long, select = c("accepted_plant_name" , "area_sq"))
L3.area.sum <- aggregate(. ~ accepted_plant_name, L3.area.only, sum)
L3.area.sum <- as.data.frame(L3.area.sum)
wide.sq <- wcsp.data
wide.sq[] <- lapply(wcsp.data, function(x) L3.area.sum$area_sq[match(x, L3.area.sum$accepted_plant_name)])

# Convert spec.list and wcsp.data into df, as the loop can't handle it as a list
spec.list.equal.length <- lapply(spec.list, 'length<-', max(lengths(spec.list)))
spec.df <- as.data.frame(spec.list.equal.length)
wcsp.data <- as.data.frame(wcsp.data)
wide.sq <- as.data.frame(wide.sq)


# Create vectors to store results
result.list.total <- c()
result.list.relative <- c()
result.list.diversity <- c()
result.list.distrange.mean <- c()
result.list.distrange.median <- c()
result.list.range.area <- c()
result.list.endemism.total <- c()
result.list.endemism.relative <- c()
result.list.bien <- c()
result.list.endemics.list <- c()
result.list.endemics.seq <- c()
result.list.gen.dup <- c()




# ---- Fill vectors with data for desired variables ---- #




# Number of species in each region with genetic data available (Phylogenetic effort)
for (i in area.vector){
  vector01 <- c(as.character(wcsp.data[,i]))
  result.list.total[[i]] <- length(intersect(genbank.vector, vector01))
}

# Number of species with genetic data relative to total number of species in each region (Phylogenetic knowledge, main response variable)
for (i2 in area.vector){
  vector02 <- c(as.character(wcsp.data[,i2]))
  result.list.relative[[i2]] <- length(intersect(genbank.vector, vector02))/length(vector02[!is.na(vector02)])
}

# Total species richness
for (i3 in area.vector){
  vector03 <- c(as.character(wcsp.data[,i3]))
  result.list.diversity[[i3]] <- length(vector03[!is.na(vector03)])
}

# Mean number of regions in which a species occur for each region (mean range)
for (i4 in area.vector){
  vector04 <- c(na.omit(as.character(wcsp.data[,i4])))
  result.list.distrange.mean[[i4]] <- mean(dist.range.df$Freq[match(vector04, dist.range.df$Var1)])
}

# Median number of regions in which a species occur for each region (median range)
for (i5 in area.vector){
  vector05 <- c(na.omit(as.character(wcsp.data[,i5])))
  result.list.distrange.median[[i5]] <- median(dist.range.df$Freq[match(vector05, dist.range.df$Var1)])
}

# Total area of the regions in which a species occur (range area)
for (i8 in area.vector){
  vector08 <- c(na.omit(as.numeric(wide.sq[,i8])))
  result.list.range.area[[i8]] <- mean(vector08)
}

# Number of species that only occur in the specific region (total endemics)
for (i6 in area.vector){
  vector06 <- c(na.omit(as.character(wcsp.data[,i6])))
  result.list.endemism.total[[i6]] <- sum((dist.range.df$Freq[match(vector06, dist.range.df$Var1)]) == 1)
}

# Number of endemics relative to total number of species in the region (endemics relative to SR)
for (i7 in area.vector){
  vector07 <- c(na.omit(as.character(wcsp.data[,i7])))
  result.list.endemism.relative[[i7]] <- (sum((dist.range.df$Freq[match(vector07, dist.range.df$Var1)]) == 1))/length(dist.range.df$Freq[match(vector07, dist.range.df$Var1)])
}


# Number of species with available BIEN data in a region relative to SR (Distribution Knowledge)
for (i in area.vector){
  wcsp.acc.name <- name.id$taxon_name[match(spec.df[,i], name.id$plant_name_id)]
  wcsp.total <- wcsp.data[,i]
  #result.list.bien[[i]] <- length(intersect(wcsp.acc.name, wcsp.total))/length(wcsp.total[!is.na(wcsp.total)])#version from Rudbeck et al.
  result.list.bien[[i]] <- length(na.omit(intersect(wcsp.acc.name, wcsp.total)))/length(wcsp.total[!is.na(wcsp.total)]) #corrected version

}


# List of all species that are endemic to a single botanical country
for (i9 in area.vector){
  vector09 <- c(na.omit(as.character(wcsp.data[,i9])))
  result.list.endemics.list[[i9]] <- dist.range.df.endemics$Var1[match(vector09, dist.range.df.endemics$Var1)]
}

# Format the list created above to become fit for looping
endemism.equal.length <- lapply(result.list.endemics.list, 'length<-', max(lengths(result.list.endemics.list)))
endemism.as.df <- as.data.frame(endemism.equal.length)

# % of endemic species sequenced
for (i10 in area.vector){
  vector10 <- c(as.character(endemism.as.df[,i10]))
  result.list.endemics.seq[[i10]] <- length(intersect(genbank.vector, vector10))/length(vector10[!is.na(vector10)])
}


# Inventory completeness of phylogenetic knowledge, when completeness is achieved by having available sequences for all 128 relevant markers for each species
for (i11 in area.vector){
  vector11 <- c(as.character(wcsp.data[,i11]))
  result.list.gen.dup[[i11]] <- length(genbank.vector.duplicate[genbank.vector.duplicate %in% intersect(genbank.vector.duplicate, vector11)])/((length(vector11[!is.na(vector11)]))*128)
}


# Aggregate results into table

# convert vectors to df's
total.df <- data.frame(TOTAL=unlist(result.list.total),
                       RELATIVE=unlist(result.list.relative),
                       RICHNESS=unlist(result.list.diversity),
                       RANGE_MEAN=unlist(result.list.distrange.mean),
                       RANGE_MEDIAN=unlist(result.list.distrange.median),
                       RANGE_AREA=unlist(result.list.range.area),
                       ENDEMISM_T=unlist(result.list.endemism.total),
                       ENDEMISM_R=unlist(result.list.endemism.relative),
                       BIEN_OCCUR=unlist(result.list.bien),
                       ENDEMIC_SEQ=unlist(result.list.endemics.seq),
                       GEN_DUP=unlist(result.list.gen.dup)
)
total.df <- data.frame(LEVEL_3_CO = row.names(total.df), total.df)



# ---- Create table to be used for logistic regression of the sequencing status of species (1 = sequenced, 0 = not sequenced) on their range size, measured as the sum of the areas (in km^2) of all the botanical countries a species occurs in (logit.Rmd)



L3.area.sum$sequenced <- ifelse(L3.area.sum$accepted_plant_name %in% gen.data$species, 1, 0)



# ---- Write results to be mapped and used in the statistical analysis ---- #




write.table(total.df,"data\\variables_2022.csv",sep = ",", col.names = T, row.names = F, quote = F)

write.table(endemics.seq.df,"data\\endemic_completeness_2022.csv",sep = ",", col.names = F, row.names = T, quote = F)

write.table(L3.area.sum,"data\\seq_area_2022.csv", sep = ",", col.names = T, row.names = F, quote = F)



# ---- Additional information of interest included in the paper ---- #


# % of species with phylogentically relevant data

length(intersect(unique(wcsp.long$accepted_plant_name) , unique(gen.data$species))) # total number

length(intersect(unique(wcsp.long$accepted_plant_name) , unique(gen.data$species)))/length(unique(wcsp.long$accepted_plant_name)) # %

