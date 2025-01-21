################################################################################
####                        Read in counts data                             ####
################################################################################

counts <- read.csv("Pan_china-AllCounts.csv")
View(counts)

################################################################################
####                           Exploring Data                               ####
################################################################################

summary(counts) # BJ subjects have a much lower maximum count number in general. This could be because of the AS (distribution of isoforms).
class(counts) # 
str(counts)

################################################################################
####                          Creating Metadata                             ####
################################################################################

metadata <- read.csv("Metadata_Pan_tidied.csv", header = T)
#View(metadata)

################################################################################
####                        Creating BJ Metadata                            ####
################################################################################

BJ_metadata <- metadata[1:10,]
#View(BJ_metadata)
write.csv(BJ_metadata, "BJ_metadata.csv", row.names = F)

################################################################################
####                        Creating CD Metadata                            ####
################################################################################

CD_metadata <- metadata[11:20,]
#View(CD_metadata)
write.csv(CD_metadata, "CD_metadata.csv", row.names = F)
