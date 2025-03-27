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

metadata <- read.csv("Metadata_Pan_tidied.csv", stringsAsFactors = T)
write.table(metadata, "DESeq_all_data/Metadata_Pan_tidied.txt", row.names = F, quote = F, sep = "\t")
#View(metadata)

################################################################################
####                        Creating BJ Metadata                            ####
################################################################################

BJ_metadata <- as.data.frame(metadata[1:10,])
#View(BJ_metadata)
write.table(BJ_metadata, "DESeq_BJ/BJ_metadata.txt", row.names = F, quote = F, sep = "\t")

################################################################################
####                        Creating CD Metadata                            ####
################################################################################

CD_metadata <- as.data.frame(metadata[11:20,])
#View(CD_metadata)
write.table(CD_metadata, "DESeq_CD/CD_metadata.txt", row.names = F, quote = F, sep = "\t")





