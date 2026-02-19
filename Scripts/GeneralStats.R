physeqAdistance4groupalpha

# number of fpid

distinct(physeqAdistance4groupalpha, SampleYear.x,season.x) %>% arrange(SampleYear.x)

physeqAdistance4groupalpha %>% group_by(SampleYear.x,season.x) %>% summarise(n=n()) %>% ungroup %>%
  summarise(mean=mean(n), se=sd(n)/sqrt(n()))

disperbirdseason <- physeqAdistance4groupalpha  %>% distinct(season.x,SampleYear.x,ItemID1_,ItemID2_,.keep_all = T) %>% filter(Pairs %in% "sameGroup" & Typeofpair2 %in% "Dominantpair") %>% group_by(ItemID1_,ItemID2_)  


disperbirdseason %>% summarise(n=n()) %>% arrange(desc(n))


physeqAdistance4group68 <- physeqAdistance4group %>% distinct(season.x,SampleYear.x,ItemID1_,ItemID2_,.keep_all = T) %>% filter(Typeofpair2 %in% "Dominantpair") %>% group_by(ItemID1_,ItemID2_)   %>% summarise(n=n()) %>% arrange(desc(n))

phychangestatus <- physeqAdistance4groupalpha %>% group_by(BirdID.x) %>%
  mutate(n=n_distinct(Status.x)) %>%
  distinct(BirdID.x,Status.x,.keep_all = T)


phystatus1 <- physeqAdistance4groupalpha %>%
  dplyr::select(BirdID.x,Status.x,season.x,SampleYear.x,TerritoryID.x) 

phystatus2 <- physeqAdistance4groupalpha %>%
  dplyr::select(BirdID.y,Status.y,season.y,SampleYear.y,TerritoryID.y) %>%
  dplyr::rename(BirdID.x=BirdID.y,Status.x=Status.y,season.x=season.y,SampleYear.x=SampleYear.y,TerritoryID.x=TerritoryID.y)

phystatus3 <- rbind(phystatus1,phystatus2) %>%
  mutate(Status2 = case_when(Status.x %in% c("BrM","BrF") ~ "Dom", Status.x %in% c("AB","ABX","H") ~ "Sub", TRUE ~ "Other"))

phystatus3 %>% filter(Status2 %in% "Other") %>% distinct(Status.x)

phystatus4 <- phystatus3 %>% group_by(BirdID.x, TerritoryID.x) %>%
  mutate(n=n_distinct(Status2)) %>%
  distinct(BirdID.x,Status2,.keep_all = T) %>% 
  filter(n>1) 

n_distinct(phystatus4$BirdID.x)

physeqAdistance4groupalpha %>% group_by(SampleYear.x) %>% summarise(n=n())
physeqAdistance4groupalpha %>% group_by(season.x) %>% summarise(n=n())
physeqAdistance4groupalpha %>% group_by(sexdiff) %>% summarise(n=n())
physeqAdistance4groupalpha %>% group_by(samenatal) %>% summarise(n=n())
physeqAdistance4groupalpha %>% summarise(min(agediff), max(agediff), mean(agediff),sd(agediff))
physeqAdistance4groupalpha %>% summarise(min(Timeinfridgediff), max(Timeinfridgediff), mean(Timeinfridgediff),sd(Timeinfridgediff))
physeqAdistance4groupalpha %>% summarise(min(sampledatediff), max(sampledatediff), mean(sampledatediff),sd(sampledatediff))
physeqAdistance4groupalpha %>% summarise(min(R.ped), max(R.ped), mean(R.ped),sd(R.ped))

physeq4asnewsdotu <- otu_table(physeq4asnewsd) %>% as.data.frame()
physeq4asnewsdtax <- tax_table(physeq4asnewsd) %>% as.data.frame()
physeq4asnewsdsd <- sample_data(physeq4asnewsd) %>% as.data.frame()

write_tsv(physeq4asnewsdotu, "InputTables/physeq4asnewsdotu.tsv")
write_tsv(physeq4asnewsdtax, "InputTables/physeq4asnewsdtax.tsv")
write_tsv(physeq4asnewsdsd, "InputTables/physeq4asnewsdsd.tsv")
