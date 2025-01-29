#relatedness calculation - CL - 20250129 
# install "related" with all the steps from their github page - important!
library(related)
library(readxl)
library(tidyverse)
genotypedata <- read_excel("~/Downloads/genotypedata.xlsx")
genotypedata %>% distinct(Locus)
gtd <- genotypedata %>% dplyr::select(BirdID,Locus,AlleleA,AlleleB) %>% pivot_wider(names_from = "Locus",values_from = c("AlleleA","AlleleB")) 

gtd2 <- gtd[,c(1,2,43,3,44,4,45,5,46,6,47,7,48,8,49,9,50,10,51,11,52,12,53,13,54,14,55,15,56,16,57,17,58,18,59,19,60,20,61,21,62,22,63,23,64,24,65,25,66,26,67,27,68,28,69,29,70,30,71,31,72,32,73,33,74,34,75,35,76,36,77,37,78,38,79,39,80,40,81,41,82,42,83)]

#write_delim(gtd2,"gtd.txt",delim=" ", col_names = F) # remove comment to write

swgtd <- readgenotypedata("gtd.txt")

relatednessperbird <- coancestry(swgtd$gdata, quellergt=2)

relatednessbi <- relatednessperbird$relatedness %>% dplyr::select(ind1.id,ind2.id,quellergt)
#write_csv(relatednessbi,"relatedness.csv") # remove comment to write

ggplot(relatednessbi, aes(quellergt)) + geom_histogram()
