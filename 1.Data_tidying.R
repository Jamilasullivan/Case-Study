## Read in counts data #########################################################

counts <- read.csv("Pan_china-AllCounts.csv")
View(counts)

## Explore data ################################################################

summary(counts) # BJ subjects have a much lower maximum count number in general. This could be because of the AS (distribution of isoforms).
class(counts) # 
str(counts)
