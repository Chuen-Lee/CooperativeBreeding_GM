library("readxl")
library("data.table")
library("ggdist")
library("emmeans")
library("tidyverse")
library("ggsignif")
library("tidyverse")
library("microbiome")
library("speedyseq")
library("vegan")
library("lmerTest")
library("DHARMa")
library("car")
library("ggeffects")
library("ggthemes")
library("ggpubr")

# prepare data for models


BreedStatus20240827_2 <- read_excel("InputTables/BreedStatus20240827.xlsx")
BreedGroupLocation20240827 <- read_excel("InputTables/BreedGroupLocation20240827.xlsx")
BreedStatusWide <- BreedStatus20240827_2 %>% dplyr::select(BreedGroupID,BirdID,Status) %>% merge(.,BreedGroupLocation20240827, by="BreedGroupID") %>% group_by(BirdID,FieldPeriodID) %>% dplyr::arrange(factor(Status, levels = c("BrM","BrF","H","ABX","AB","OFL","FL","CH","BrU","B","TBRF","TBRM","NSA","SEEN2","SEEN1")), .by_group = T) %>% distinct(FieldPeriodID,.keep_all = T) %>% filter(Status %in% c("BrM","BrF","H","CH","ABX","AB","OFL","FL")) %>% ungroup %>% group_by(BreedGroupID) %>% pivot_wider(values_from = "BirdID",names_from = "Status",names_sep = "_",names_vary = "fastest") %>% unnest_wider(BrM,names_sep = "_") %>% unnest_wider(BrF,names_sep = "_") %>% unnest_wider(H,names_sep = "_") %>% unnest_wider(ABX,names_sep = "_") %>% unnest_wider(AB,names_sep = "_") %>% unnest_wider(CH,names_sep = "_") %>% unnest_wider(OFL,names_sep = "_")  %>% unnest_wider(FL,names_sep = "_") %>% ungroup  %>% mutate(n_members = rowSums( !is.na( . [,2:28]))) %>% dplyr::arrange(BreedGroupID) %>% group_by(BrM_1,BrF_1) %>% dplyr::mutate(seasonstgt = case_when(is.na(BrM_1) | is.na(BrF_1) | BrM_1 %in% "0" | BrF_1 %in% "0" ~ 0, TRUE ~ 1:n())) %>% mutate(HelperPresent = case_when(!is.na(H_1) ~ "yes", TRUE ~ "no")) %>% mutate(FledglingPresent = case_when(!is.na(FL_1) | !is.na(OFL_1) ~ "yes", TRUE ~ "no"))

ggplot(BreedStatusWide,aes(n_members)) + geom_histogram()
ggplot(BreedStatusWide,aes(seasonstgt)) + geom_histogram()

BreedStatusWide2 <- BreedStatusWide %>% ungroup() %>% dplyr::select(BreedGroupID,n_members,seasonstgt,HelperPresent,FledglingPresent,TerritoryID)

firststatus <- BreedStatus20240827_2 %>% arrange(StatusID) %>% distinct(BirdID,.keep_all = T) %>% merge(.,BreedGroupLocation20240827,by="BreedGroupID") #%>% dplyr::select(BirdID,TerritoryID) %>% rename(firstTerritoryID = TerritoryID)

firststatusbreedgroups <- firststatus %>% dplyr::select(BreedGroupID) %>% merge(.,BreedStatus20240827_2,by="BreedGroupID",all.x=T) %>% filter(Status %in% c("BrM","BrF","H")) %>% dplyr::select(BreedGroupID,Status,BirdID) %>% group_by(BreedGroupID,Status) %>% distinct(BirdID,.keep_all = T) %>% pivot_wider(names_from = "Status",values_from = "BirdID",names_sep = "_",names_vary = "fastest") %>% unnest_wider(BrM,names_sep = "_") %>% unnest_wider(BrF,names_sep = "_") %>% unnest_wider(H,names_sep = "_")

firststatusindividuals <- merge(firststatus,firststatusbreedgroups,by="BreedGroupID",all.x=T) %>% dplyr::select(BirdID,TerritoryID,contains("BrF"),contains("BrM"),contains("H")) %>% rename(firstTerritoryID = TerritoryID)

Offspring20240530 <- read_csv("InputTables/Offspring20240530.csv",  col_types = cols(BirthDate = col_date(format = "%d/%m/%Y")))
Offspring20240530ind <- Offspring20240530 %>% filter(Confidence>80) %>% dplyr::select(Parent,OffID,SpouseID) %>% group_by(OffID) %>% distinct(OffID,.keep_all = T) %>% filter(!is.na(Parent))

relatednessbi <- read_tsv("InputTables/Rped.old")


## functions ----
vegan_otu <- function(physeq) {
  OTU <- otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}

alpha_estimate <- function(physeq) {
  ttt <- transform_sample_counts(physeq, function(x) round(x, 0))
  richnessEstRare<-estimate_richness(ttt, split=TRUE, measures= c("Chao1", "Shannon", "Observed"))
  RareMeta <- as.data.frame(sample_data(ttt))
  RareMeta$Chao1 <- richnessEstRare$Chao1
  RareMeta$Shannon <- richnessEstRare$Shannon
  RareMeta$Observed <- richnessEstRare$Observed
  AlphaMeta <- as_tibble(RareMeta) %>% mutate(TerminalYear = as.factor(TerminalYear)) %>% mutate(Timeinfridge = as.numeric(Timeinfridge)) %>% 
    mutate(CatchTime = tidyr::replace_na(CatchTime, 43292)) %>% mutate(CatchTime = as.numeric(CatchTime))
  return(AlphaMeta)
}

# read physeq object from Sarah's senescence paper ----
physeq3 <- readRDS("InputTables/physeq3.rds") # generated from "BetaAnalysis.R" in Sarah Worsley's paper
physeq4as <- physeq3 %>% filter_sample_data(Ageclass %in% c("A","SA","OFL"))


physeq4assd <- data.frame(sample_data(physeq4as)) %>% rownames_to_column("Rownames") %>% dplyr::rename(season=Season,Timeinfridge=TimeAt4Degrees,CatchTime=MinutesSinceSunrise, SamplingAge = AgeYears) %>% mutate(Timeinfridge_scaled = arm::rescale(Timeinfridge),CatchTime_scaled = arm::rescale(CatchTime), TQcorrected_scaled = arm::rescale(TQcorrected) ) %>% merge(., BreedStatusWide2, by=c("BreedGroupID"),all.x=T) %>% column_to_rownames("Rownames")
physeq4asnewsd <- phyloseq(tax_table(physeq4as),otu_table(physeq4as),sample_data(physeq4assd), phy_tree(physeq4as)) 
physeqrarefied <- rarefy_even_depth(physeq4asnewsd, rngseed=88,sample.size = 8000, replace=F )

physeqrarefiedalpha <- alpha_estimate(physeqrarefied)
physeqrarefiedalpha2 <- physeqrarefiedalpha %>% dplyr::select(sample.id,Observed,Shannon) %>% dplyr::rename(SampleID=sample.id) 

# beta diversity ----
physeqrare <- transform(physeq4as, "compositional") %>% core_members(.,detection = 0.0001, prevalence= 0.05)
physeqclr <- prune_taxa(physeqrare,physeq4as) %>% transform(.,"clr")
physeq_clr_Asddd <- data.frame(sample_data(physeqclr)) %>% rownames_to_column("Rownames") %>% dplyr::rename(season=Season,Timeinfridge=TimeAt4Degrees,CatchTime=MinutesSinceSunrise, SamplingAge = AgeYears) %>% mutate(Timeinfridge_scaled = arm::rescale(Timeinfridge),CatchTime_scaled = arm::rescale(CatchTime), TQcorrected_scaled = arm::rescale(TQcorrected) ) %>% merge(., BreedStatusWide2, by=c("BreedGroupID"),all.x=T) %>% column_to_rownames("Rownames")

PhyseqAbreeddd <- phyloseq(tax_table(physeqclr),otu_table(physeqclr),sample_data(physeq_clr_Asddd), phy_tree(physeqclr)) 

# paired distances ----
physeqAdistance <- phyloseq::distance(PhyseqAbreeddd, method = "euclidean", type="samples") %>% as.matrix() 
physeqAdistance[lower.tri(physeqAdistance)] <- 0

physeqAdistance <- physeqAdistance %>% data.frame() %>% rownames_to_column("SampleID") %>% pivot_longer(.,cols=-c(SampleID)) %>% filter(value > 0) %>% mutate(value=-value)
physeq_Abreedclr_st2 <- data.frame(sample_data(PhyseqAbreeddd)) %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID,BreedGroupID,BirdID,season,SampleYear,Status,n_members,seasonstgt,HelperPresent,FledglingPresent,CatchTime,Timeinfridge,SampleDate,SamplingAge,SexEstimate,TerritoryID)

physeqAdistance2 <- merge(physeqAdistance, physeq_Abreedclr_st2, by.x=("SampleID"), by.y="SampleID") %>% merge(.,physeq_Abreedclr_st2, by.x="name",by.y="SampleID") %>% merge(.,physeqrarefiedalpha2, by="SampleID",all.x=T) %>% merge(.,physeqrarefiedalpha2, by.x="name",by.y="SampleID",all.x=T) 

physeqAdistance3 <- physeqAdistance2 %>% mutate(Pairs = as.factor(case_when(name == SampleID ~ "samesample" ,as.numeric(BirdID.x) == as.numeric(BirdID.y) ~ "sameBird", as.numeric(BreedGroupID.x) == as.numeric(BreedGroupID.y) ~ "sameGroup", as.numeric(BreedGroupID.x) != as.numeric(BreedGroupID.y) ~ "diffGroup"))) %>% filter(!(Pairs %in% "samesample")) %>% mutate(seasonpairs = as.factor(case_when(season.x == season.y & SampleYear.x == SampleYear.y ~ "sameyearsameseason",season.x == season.y ~ "sameseasondiffyear",SampleYear.x == SampleYear.y ~ "sameyeardiffseason", TRUE ~ "diffyeardiffseason")), noseasonstgt = case_when(Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF") & seasonstgt.x == seasonstgt.y & Pairs %in% "sameGroup" ~ seasonstgt.x)) %>% mutate(SampleYear.x = as.factor(SampleYear.x)) %>% filter(Status.x %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H") & Status.y %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H")) %>% mutate(CatchTimediff = abs(CatchTime.x - CatchTime.y), Timeinfridgediff = abs(Timeinfridge.x - Timeinfridge.y), sampledatediff = abs(as.numeric(as.Date(SampleDate.x,"%d/%m/%Y") - as.Date(SampleDate.y,"%d/%m/%Y"))), agediff = abs(SamplingAge.x-SamplingAge.y), sexdiff = case_when(SexEstimate.x == SexEstimate.y ~ "0", TRUE ~ "1"))

physeqAdistance4 <- physeqAdistance3 %>% filter(seasonpairs %in% "sameyearsameseason")  %>% mutate(Typeofpair2 = as.factor(case_when( (Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF")) & Pairs %in% "sameGroup" ~ "Dominantpair", (Status.x %in% c("BrM","BrF") & Status.y %in% c("H") | Status.y %in% c("BrM","BrF") & Status.x %in% c("H")) & Pairs %in% "sameGroup" ~ "DomHelp" ,  (Status.x %in% c("BrM","BrF") & Status.y %in% c("AB","ABX") | Status.y %in% c("BrM","BrF") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "DomSub", (Status.x %in% c("H") & Status.y %in% c("AB","ABX") | Status.y %in% c("H") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "HelpSub", (Status.x %in% c("H","AB","ABX") & Status.y %in% c("H","AB","ABX")) & Pairs %in% "sameGroup"  ~ "SubSub",Pairs %in% "sameGroup" ~ "sameTerrothers", Pairs %in% "sameBird" ~ "sameBird",TRUE ~ Pairs))) %>% mutate(Typeofpair2 = factor(Typeofpair2,level=c("diffGroup","Dominantpair","DomHelp","DomSub","HelpSub","SubSub","sameTerrothers","sameBird")), Pairs=factor(Pairs, level=c("diffGroup","sameGroup","sameBird"))) %>% filter(!(Typeofpair2 %in% "sameTerrothers")) %>% filter(!(Pairs %in% "sameBird"))

physeqAdistance4group <- physeqAdistance4 %>%
  mutate(ItemID1_ = pmin(BirdID.x  ,BirdID.y),
         ItemID2_ = pmax(BirdID.x  ,BirdID.y)) %>%
  group_by(ItemID1_,ItemID2_) %>% mutate(BirdGroupID = cur_group_id()) %>% 
  merge(.,relatednessbi,by.x=c("ItemID1_","ItemID2_"),by.y = c("IID1","IID2"),all.x=T) %>% 
  merge(.,firststatusindividuals,by.x = "BirdID.x",by.y="BirdID",all.x=T) %>% 
  merge(.,firststatusindividuals,by.x = "BirdID.y",by.y="BirdID",all.x=T) %>% 
  mutate(samenatal = case_when( (BirdID.y == BrF_1.x | BirdID.y== BrM_1.x | BirdID.y == H_1.x | BirdID.y == BrF_2.x | BirdID.y== BrM_2.x | BirdID.y == H_2.x |BirdID.y == H_3.x | BirdID.x == BrF_1.y | BirdID.x== BrM_1.y | BirdID.y == H_1.y | BirdID.y == BrF_2.y | BirdID.y== BrM_2.y | BirdID.y == H_2.y |BirdID.y == H_3.y ) ~ "same", TRUE ~ "different" )) %>% merge(., Offspring20240530ind,by.x="BirdID.x", by.y = "OffID",all.x=T) 

#write_csv(physeqAdistance4group,"Output/physeqAdistance4group.csv")
#physeqAdistance4group <- read_csv("Output/physeqAdistance4group.csv") %>% mutate(Typeofpair2 = factor(Typeofpair2,level=c("diffGroup","Dominantpair","DomHelp","DomSub","HelpSub","SubSub","sameTerrothers","sameBird")), Pairs=factor(Pairs, level=c("diffGroup","sameGroup","sameBird")))


physeqAdistance4groupsamepair <- physeqAdistance4group %>% filter(Pairs%in%"sameGroup") 

physeqAdistance4groupsamepair %>% summarise(mean(R.ped), min(R.ped),max(R.ped), sd=sd(R.ped))

physeqAdistance4groupsamepair %>% group_by(samenatal,Typeofpair2) %>% summarise(n=n())

physeqAdistance4groupsamepairhelp <- physeqAdistance4groupsamepair %>% filter(Typeofpair2%in%"DomHelp" & samenatal%in%"different")

physeqAdistance4sss <- physeqAdistance4group %>% group_by(Pairs) %>% dplyr::summarise(count=n()) %>% mutate(value = 0)
physeqAdistance4ss <- physeqAdistance4group %>% filter(Pairs %in% "sameGroup") %>% group_by(Typeofpair2) %>% dplyr::summarise(count=n()) %>% mutate(value = 0) 

## descriptive stats ----
# number of samples / number of bird id
physeqAdistance4group %>% ungroup() %>% reframe(n=n())
physeqAdistance4group %>% ungroup() %>% dplyr::select(SampleID,name) %>% pivot_longer(cols = everything()) %>% distinct(value) %>% summarise(n()) # 648 samples
physeqAdistance4group %>% ungroup() %>% dplyr::select(BirdID.x,BirdID.y) %>% pivot_longer(cols = everything()) %>% distinct(value) %>% summarise(n()) # 345 birds
physeqAdistance4group %>% ungroup() %>% dplyr::select(SampleID,name,Pairs) %>% pivot_longer(cols = c("SampleID","name")) %>% group_by(Pairs) %>% distinct(value,.keep_all = T) %>% summarise(n()) # 648 samples diff Group, 322 samples same group, 188 samples samebird
physeqAdistance4group %>% ungroup() %>% dplyr::select(BirdID.x,BirdID.y,Pairs) %>% pivot_longer(cols = c("BirdID.x","BirdID.y")) %>% group_by(Pairs) %>% distinct(value,.keep_all = T) %>% summarise(n()) # 345 birds diff group, 204 birds same group, 86 birds same bird

physeqAdistance4groupsamepair %>% ungroup() %>% reframe(n=n()) # 279 pairwise comparisons
physeqAdistance4groupsamepair %>% ungroup() %>% dplyr::select(SampleID,name) %>% pivot_longer(cols = c("SampleID","name"))%>% distinct(value,.keep_all = T) %>% summarise(n()) # 322 samples
physeqAdistance4groupsamepair %>% ungroup() %>% dplyr::select(BirdID.x,BirdID.y) %>% pivot_longer(cols = c("BirdID.x","BirdID.y"))%>% distinct(value,.keep_all = T) %>% summarise(n()) # 204 birds
physeqAdistance4groupsamepair %>% ungroup() %>% dplyr::select(SampleID,name,Typeofpair2) %>% pivot_longer(cols = c("SampleID","name")) %>% group_by(Typeofpair2) %>% distinct(value,.keep_all = T) %>% summarise(n()) # 121 samples dompair, 69 samples domhelp, 156 samples domsub, 21 samples helpsub, 46 samples helpsub
physeqAdistance4groupsamepair %>% ungroup() %>% dplyr::select(BirdID.x,BirdID.y,Typeofpair2) %>% pivot_longer(cols = c("BirdID.x","BirdID.y")) %>% group_by(Typeofpair2) %>% distinct(value,.keep_all = T) %>% summarise(n()) # 93 birds dompair, 54 birds domhelp, 107 birds domsub, 15 birds helpsub, 35 birds subsub
#physeqAdistance4group %>% ungroup() %>% dplyr::select(SamplingAge.x,SamplingAge.y) %>% pivot_longer(cols = everything()) %>% summarise(minage = min(value),maxage = max(value), meanage = mean(value), sdage = sd(value), se = sd(value)/sqrt(n()))

physeqAdistance4group %>% ungroup() %>% distinct(BirdGroupID,.keep_all = T) %>% summarise(minagediff = min(agediff), maxagediff = max(agediff), meanagediff = mean(agediff), sdagediff = sd(agediff), sediff = sd(agediff)/sqrt(n())) # mean age diff = 3.07, sd = 3.01, se = 0.0225
physeqAdistance4group %>% ungroup() %>% distinct(BirdGroupID,.keep_all = T) %>% summarise(mindays = min(Timeinfridgediff), maxagediff = max(Timeinfridgediff), meanagediff = mean(Timeinfridgediff), sdagediff = sd(Timeinfridgediff), sediff = sd(Timeinfridgediff)/sqrt(n())) # timein fridge diff, mean = 23.7, sd = 19.3, se = 0.144
physeqAdistance4group %>% ungroup() %>% distinct(BirdGroupID,.keep_all = T) %>% summarise(mindays = min(CatchTimediff), maxagediff = max(CatchTimediff), meanagediff = mean(CatchTimediff), sdagediff = sd(CatchTimediff), sediff = sd(CatchTimediff)/sqrt(n()), n = n()) # catch time diff, mean = 193, sd = 147, se = 1.09

# paired alpha diversity ----
physeqAdistance4groupalpha <- physeqAdistance4group %>% mutate(Observeddiff = -abs(Observed.x-Observed.y), Shannondiff = -abs(Shannon.x-Shannon.y))
physeqAdistance4groupalpha2 <- physeqAdistance4groupalpha %>% filter(Pairs %in% "sameGroup")



# paired beta diversity ----
# within group vs between group ----



# type of pair ----

physeqAdistance4groupsamepair %>% group_by(Typeofpair2) %>% summarise(n())



# genus abundance and aerobic/anaerobic ----
physeq4asnewsdgen <- tax_glom(physeq4asnewsd, taxrank="Genus", NArm=T) 
taxphyseqclr <- tax_table(physeq4asnewsdgen) %>% data.frame() 

taxphyseqclrgenus <- taxphyseqclr %>% distinct(Genus)
taxphyseqclrgenus %>% summarise(n=n())
#write_csv(taxphyseqclrgenus, "taxphyseqclrgenus.csv")
taxaassignchuen <- read_csv("InputTables/taxphyseqclrgenus.csv")
chatgptassign <- read_csv("InputTables/genus_chatgpt.csv")
knowleslabannote <- read_csv("InputTables/41559_2024_2381_MOESM3_ESM.csv")

checkchuenchat <- merge(chatgptassign,taxaassignchuen,by="Genus",all=T)%>% filter(!is.na(Aerotolerance) & !is.na(Oxygen)) %>% mutate(conflict = case_when((Oxygen %in% c("Aerobic","Facultative Anaerobic") & Aerotolerance %in% c( "Aerobic","facultatively anaerobic","mixed")) | (Oxygen %in% "Anaerobic" & Aerotolerance %in% "Anaerobe") | Aerotolerance%in% "unknown" ~ "0", TRUE ~ "1"))
checkchuenchat %>% group_by(conflict) %>% summarise(nconflict = n())
47/(17+47) # 73.4%
checkknowleschat <- merge(chatgptassign,knowleslabannote,by="Genus") %>% mutate(conflict = case_when((Oxygen %in% c("Aerobic","Facultatively anaerobic") & Aerotolerance_binary %in% "aerotolerant") | (Oxygen %in% "Anaerobic" & Aerotolerance_binary %in% "anaerobic") | (Oxygen %in% "Unknown" & Aerotolerance_binary %in% "unknown") | Aerotolerance_binary %in% "unknown" ~ "0", TRUE ~ "1"))
checkknowleschat %>% group_by(conflict) %>% summarise(nconflict = n())
98/(98+34) #74.2%

geminiassign <- read_csv("InputTables/genus_gemini.csv")
geminiassign %>% group_by(Aero) %>% summarise(n=n()) # Aerobic 651, Anaerobic 224, Facultative Anaerobic 65, Not applicable 2, Unknown 133, Variable 36

checkchuengemini <- merge(geminiassign,taxaassignchuen,by="Genus",all=T) %>% filter(!is.na(Aerotolerance)) %>% mutate(conflict = case_when((Aero %in% c("Aerobic","Facultative Anaerobic") & Aerotolerance %in% c( "Aerobic","facultatively anaerobic","mixed")) | (Aero %in% "Anaerobic" & Aerotolerance %in% "Anaerobe") ~ "0", TRUE ~ "1"))
checkknowlesgemini <- merge(geminiassign,knowleslabannote,by="Genus") %>% mutate(conflict = case_when((Aero %in% c("Aerobic","Facultative Anaerobic") & Aerotolerance_binary %in% c("aerotolerant","mixed")) | (Aero %in% "Anaerobic" & Aerotolerance_binary %in% "anaerobic")  ~ "0", TRUE ~ "1")) %>% filter(!(Aerotolerance_binary %in% "unknown") & !(Aero %in% "Unknown"))

checkchuengemini %>% group_by(conflict) %>% summarise(nconflict = n())
77/(77+3) # 96.3% accuracy in 80 genera
checkknowlesgemini %>% group_by(conflict) %>% summarise(nconflict = n())
160/(160+13) # 92.5% accuracy in 173 genera
161/(161+10)

taxphyseqclrgenusannote <- merge(taxaassignchuen,knowleslabannote, by="Genus",all.x=T) %>% filter(!is.na(Aerotolerance) & !is.na(Aerotolerance_binary)) %>% mutate(conflict = case_when((Aerotolerance %in% c("Aerobic","facultatively anaerobic") & Aerotolerance_binary %in% "aerotolerant") | (Aerotolerance %in% "Anaerobe" & Aerotolerance_binary %in% "anaerobic") | (Aerotolerance %in% "Unknown" & Aerotolerance_binary %in% "unknown") | Aerotolerance_binary %in% "unknown" ~ "0", TRUE ~ "1"))
# i got two wrong and aura got one wrong, out of 29

taxphyseqclrgenusannote %>% filter(is.na(Aerotolerance_binary)) %>% nrow()

taxaphyseqclr <- transform(physeq4asnewsdgen,"clr") %>% psmelt() %>% merge(.,geminiassign, by="Genus", all.x=T ) 

taxaphyseqclr %>% distinct(Genus) %>% nrow()
taxaphyseqclr %>% distinct(Aero) 
taxaphyseqclr %>% distinct(Genus,.keep_all = T) %>% group_by(Aero) %>% summarise(n())

taxaphyseqclrsum <- taxaphyseqclr%>% group_by(Sample,Aero) %>% summarise(sumgenusoxy = sum(Abundance))
ggplot(taxaphyseqclrsum, aes(x=Aero,y=sumgenusoxy)) + geom_boxplot()

taxaiassign <- merge(taxphyseqclrgenus, chatgptassign, by="Genus",all.x=T) %>%
  dplyr::rename(ChatGPT = Oxygen) %>%
  dplyr::select(-SporeForming) %>%
  merge(., geminiassign, by="Genus",all.x=T) %>%
  dplyr::rename(Gemini = Aero) %>%
  dplyr::select(-SporeForming) %>% 
  merge(., taxaassignchuen, by="Genus",all.x=T) %>%
  dplyr::rename(Chuen = Aerotolerance) %>%
  dplyr::select(-Spore_formation,-Reference,-Indirect,-Comments) %>% 
  merge(., knowleslabannote, by="Genus",all.x=T) %>%
  dplyr::rename(Aura = Aerotolerance_binary) %>%
  dplyr::select(-Spore_formation,-Reference,-indirect_inference_aero,-Comments, -indirect_inference_spor, -Detected_in_mouse_gut, -`Detected_in_local_soil;`,-...1) %>%
  mutate(ChatGPT = case_when(ChatGPT %in% "Anaerobic" ~ "Anaerobic", ChatGPT %in% c("Facultatively anaerobic", "Aerobic") ~ "Aerotolerant", !is.na(ChatGPT) ~ "Unknown", TRUE ~ "Unassigned")) %>%
  mutate(Gemini = case_when(Gemini %in% "Anaerobic" ~ "Anaerobic", Gemini %in% c("Facultative Anaerobic", "Aerobic") ~ "Aerotolerant", !is.na(Gemini) ~ "Unknown", TRUE ~ "Unassigned")) %>%
  mutate(Chuen = case_when(Chuen %in% "Anaerobe" ~ "Anaerobic", Chuen %in% c("facultatively anaerobic", "Aerobic","mixed") ~ "Aerotolerant", !is.na(Chuen) ~ "Unknown", TRUE ~ "Unassigned")) %>%
  mutate(Aura = case_when(Aura %in% "anaerobic" ~ "Anaerobic", Aura %in% c("aerotolerant", "Aerobic","mixed") ~ "Aerotolerant", !is.na(Aura) ~ "Unknown", TRUE ~ "Unassigned")) %>%
  mutate(GeminiCL = case_when(Gemini %in% "Anaerobic" & Chuen %in% "Anaerobic" ~ "Correct",Gemini %in% "Aerotolerant" & Chuen %in% "Aerotolerant" ~ "Correct", Gemini %in% "Anaerobic" & Chuen %in% "Aerotolerant" ~ "Wrong", Gemini %in% "Aerotolerant" & Chuen %in% "Anaerobic" ~ "Wrong", TRUE ~ "Not tested" )) %>%
  mutate(GeminiAura = case_when(Gemini %in% "Anaerobic" & Aura %in% "Anaerobic" ~ "Correct",Gemini %in% "Aerotolerant" & Aura %in% "Aerotolerant" ~ "Correct", Gemini %in% "Anaerobic" & Aura %in% "Aerotolerant" ~ "Wrong", Gemini %in% "Aerotolerant" & Aura %in% "Anaerobic" ~ "Wrong", TRUE ~ "Not tested" ))

taxaiassign %>% group_by(GeminiCL) %>% summarise(n=n())
taxaiassign %>% group_by(GeminiAura) %>% summarise(n=n())

write_csv(taxaiassign, "taxaiassign2.csv")

# tax table
physeq4asnewsdtaxtable <- data.frame(tax_table(physeq4asnewsdgen)) %>% rownames_to_column("OTU") %>% merge(.,geminiassign, by="Genus", all.x=T ) 

physeq4asnewsdtaxtable %>% group_by(Aero) %>% summarise(n=n())

physeqrareaero <- transform(physeq4asnewsdgen, "compositional") %>% core_members(.,detection = 0.0001, prevalence= 0.05)
physeq4asnewsdgenclr <-  transform(physeq4asnewsd, "clr")
# aerobic only  ----
aerobicasv <- physeq4asnewsdtaxtable %>% filter(Aero %in% c("Aerobic","Facultative Anaerobic")) %>% distinct(OTU)
taxaphyseqclraerobic <- prune_taxa(aerobicasv$OTU, physeq4asnewsd) 
taxaphyseqclraerobicrare <- transform(taxaphyseqclraerobic, "compositional") %>% core_members(.,detection = 0.00001, prevalence= 0.01)
taxaphyseqclraerobeta <- prune_taxa(taxaphyseqclraerobicrare, taxaphyseqclraerobic)%>% transform(.,"clr")


# paired distances -
physeqAdistanceaero <- phyloseq::distance(taxaphyseqclraerobeta, method = "euclidean", type="samples") %>% as.matrix() 
physeqAdistanceaero[lower.tri(physeqAdistanceaero)] <- 0

physeqAdistanceaero <- physeqAdistanceaero %>% data.frame() %>% rownames_to_column("SampleID") %>% pivot_longer(.,cols=-c(SampleID)) %>% filter(value > 0) %>% mutate(value=-value)
physeq_Abreedclr_st2aero <- data.frame(sample_data(PhyseqAbreeddd)) %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID,BreedGroupID,BirdID,season,SampleYear,Status,n_members,seasonstgt,HelperPresent,FledglingPresent,CatchTime,Timeinfridge,SampleDate,SamplingAge,SexEstimate)

physeqAdistanceaero2 <- merge(physeqAdistanceaero, physeq_Abreedclr_st2aero, by.x=("SampleID"), by.y="SampleID") %>% merge(.,physeq_Abreedclr_st2, by.x="name",by.y="SampleID") %>% merge(.,physeqrarefiedalpha2, by="SampleID",all.x=T) %>% merge(.,physeqrarefiedalpha2, by.x="name",by.y="SampleID",all.x=T) 

physeqAdistanceaero3 <- physeqAdistanceaero2 %>% mutate(Pairs = as.factor(case_when(name == SampleID ~ "samesample" ,as.numeric(BirdID.x) == as.numeric(BirdID.y) ~ "sameBird", as.numeric(BreedGroupID.x) == as.numeric(BreedGroupID.y) ~ "sameGroup", as.numeric(BreedGroupID.x) != as.numeric(BreedGroupID.y) ~ "diffGroup"))) %>% filter(!(Pairs %in% "samesample")) %>% mutate(seasonpairs = as.factor(case_when(season.x == season.y & SampleYear.x == SampleYear.y ~ "sameyearsameseason",season.x == season.y ~ "sameseasondiffyear",SampleYear.x == SampleYear.y ~ "sameyeardiffseason", TRUE ~ "diffyeardiffseason")), noseasonstgt = case_when(Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF") & seasonstgt.x == seasonstgt.y & Pairs %in% "sameGroup" ~ seasonstgt.x)) %>% mutate(SampleYear.x = as.factor(SampleYear.x)) %>% filter(Status.x %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H") & Status.y %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H")) %>% mutate(CatchTimediff = abs(CatchTime.x - CatchTime.y), Timeinfridgediff = abs(Timeinfridge.x - Timeinfridge.y), sampledatediff = abs(as.numeric(as.Date(SampleDate.x,"%d/%m/%Y") - as.Date(SampleDate.y,"%d/%m/%Y"))), agediff = SamplingAge.x-SamplingAge.y, sexdiff = case_when(SexEstimate.x == SexEstimate.y ~ "0", TRUE ~ "1"))

physeqAdistanceaero4 <- physeqAdistanceaero3 %>% filter(seasonpairs %in% "sameyearsameseason")  %>% mutate(Typeofpair2 = as.factor(case_when( (Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF")) & Pairs %in% "sameGroup" ~ "Dominantpair", (Status.x %in% c("BrM","BrF") & Status.y %in% c("H") | Status.y %in% c("BrM","BrF") & Status.x %in% c("H")) & Pairs %in% "sameGroup" ~ "DomHelp" ,  (Status.x %in% c("BrM","BrF") & Status.y %in% c("AB","ABX") | Status.y %in% c("BrM","BrF") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "DomSub", (Status.x %in% c("H") & Status.y %in% c("AB","ABX") | Status.y %in% c("H") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "HelpSub", (Status.x %in% c("H","AB","ABX") & Status.y %in% c("H","AB","ABX")) & Pairs %in% "sameGroup"  ~ "SubSub",Pairs %in% "sameGroup" ~ "sameTerrothers", Pairs %in% "sameBird" ~ "sameBird",TRUE ~ Pairs))) %>% mutate(Typeofpair2 = factor(Typeofpair2,level=c("diffGroup","Dominantpair","DomHelp","DomSub","HelpSub","SubSub","sameTerrothers","sameBird")), Pairs=factor(Pairs, level=c("diffGroup","sameGroup","sameBird"))) %>% filter(!(Typeofpair2 %in% "sameTerrothers"))%>% filter(!(Pairs %in% "sameBird"))

physeqAdistance4groupaero <- physeqAdistanceaero4 %>%
  mutate(ItemID1_ = pmin(BirdID.x  ,BirdID.y),
         ItemID2_ = pmax(BirdID.x  ,BirdID.y)) %>%
  group_by(ItemID1_,ItemID2_) %>% mutate(BirdGroupID = cur_group_id())%>% 
  merge(.,relatednessbi,by.x=c("ItemID1_","ItemID2_"),by.y = c("IID1","IID2"),all.x=T) %>% 
  merge(.,firststatusindividuals,by.x = "BirdID.x",by.y="BirdID",all.x=T) %>% 
  merge(.,firststatusindividuals,by.x = "BirdID.y",by.y="BirdID",all.x=T) %>% 
  mutate(samenatal = case_when( (BirdID.y == BrF_1.x | BirdID.y== BrM_1.x | BirdID.y == H_1.x | BirdID.y == BrF_2.x | BirdID.y== BrM_2.x | BirdID.y == H_2.x |BirdID.y == H_3.x | BirdID.x == BrF_1.y | BirdID.x== BrM_1.y | BirdID.y == H_1.y | BirdID.y == BrF_2.y | BirdID.y== BrM_2.y | BirdID.y == H_2.y |BirdID.y == H_3.y ) ~ "same", TRUE ~ "different" )) %>% merge(., Offspring20240530ind,by.x="BirdID.x", by.y = "OffID",all.x=T)



physeqAdistance4groupaero2 <- physeqAdistance4groupaero %>% filter(Pairs %in% "sameGroup")



# anaerobes only ----
anaasv <- physeq4asnewsdtaxtable %>% filter(Aero %in% c("Anaerobic")) %>% distinct(OTU)
taxaphyseqclranae <- prune_taxa(anaasv$OTU, physeq4asnewsd) #%>% transform(.,"clr")
taxaphyseqclranaerare <- transform(taxaphyseqclranae, "compositional") %>% core_members(.,detection = 0.00001, prevalence= 0.01)
taxaphyseqclranaebeta <- prune_taxa(taxaphyseqclranaerare, taxaphyseqclranae)%>% transform(.,"clr")

# alpha 
taxaphyseqclranaerared <- prune_taxa(anaasv$OTU,physeqrarefied)
taxaphyseqclranaeraredalpha <- alpha_estimate(taxaphyseqclranaerared)
taxaphyseqclranaeraredalpha2 <- taxaphyseqclranaeraredalpha %>% dplyr::select(sample.id,Observed,Shannon) %>% dplyr::rename(SampleID=sample.id)
# paired distances -
physeqAdistanceanae <- phyloseq::distance(taxaphyseqclranaebeta, method = "euclidean", type="samples") %>% as.matrix() 
physeqAdistanceanae[lower.tri(physeqAdistanceanae)] <- 0

physeqAdistanceanae <- physeqAdistanceanae %>% data.frame() %>% rownames_to_column("SampleID") %>% pivot_longer(.,cols=-c(SampleID))%>% filter(value > 0) %>% mutate(value=-value)
physeq_Abreedclr_st2anae <- data.frame(sample_data(PhyseqAbreeddd)) %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID,BreedGroupID,BirdID,season,SampleYear,Status,n_members,seasonstgt,HelperPresent,FledglingPresent,CatchTime,Timeinfridge,SampleDate,SamplingAge,SexEstimate)

physeqAdistanceanae2 <- merge(physeqAdistanceanae, physeq_Abreedclr_st2anae, by.x=("SampleID"), by.y="SampleID") %>% merge(.,physeq_Abreedclr_st2, by.x="name",by.y="SampleID") %>% merge(.,taxaphyseqclranaeraredalpha2, by="SampleID",all.x=T) %>% merge(.,taxaphyseqclranaeraredalpha2, by.x="name",by.y="SampleID",all.x=T) 

physeqAdistanceanae3 <- physeqAdistanceanae2 %>% mutate(Pairs = as.factor(case_when(name == SampleID ~ "samesample" ,as.numeric(BirdID.x) == as.numeric(BirdID.y) ~ "sameBird", as.numeric(BreedGroupID.x) == as.numeric(BreedGroupID.y) ~ "sameGroup", as.numeric(BreedGroupID.x) != as.numeric(BreedGroupID.y) ~ "diffGroup"))) %>% filter(!(Pairs %in% "samesample")) %>% mutate(seasonpairs = as.factor(case_when(season.x == season.y & SampleYear.x == SampleYear.y ~ "sameyearsameseason",season.x == season.y ~ "sameseasondiffyear",SampleYear.x == SampleYear.y ~ "sameyeardiffseason", TRUE ~ "diffyeardiffseason")), noseasonstgt = case_when(Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF") & seasonstgt.x == seasonstgt.y & Pairs %in% "sameGroup" ~ seasonstgt.x)) %>% mutate(SampleYear.x = as.factor(SampleYear.x)) %>% filter(Status.x %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H") & Status.y %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H")) %>% mutate(CatchTimediff = abs(CatchTime.x - CatchTime.y), Timeinfridgediff = abs(Timeinfridge.x - Timeinfridge.y), sampledatediff = abs(as.numeric(as.Date(SampleDate.x,"%d/%m/%Y") - as.Date(SampleDate.y,"%d/%m/%Y"))), agediff = SamplingAge.x-SamplingAge.y, sexdiff = case_when(SexEstimate.x == SexEstimate.y ~ "0", TRUE ~ "1"))

physeqAdistanceanae4 <- physeqAdistanceanae3 %>% filter(seasonpairs %in% "sameyearsameseason")  %>% mutate(Typeofpair2 = as.factor(case_when( (Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF")) & Pairs %in% "sameGroup" ~ "Dominantpair", (Status.x %in% c("BrM","BrF") & Status.y %in% c("H") | Status.y %in% c("BrM","BrF") & Status.x %in% c("H")) & Pairs %in% "sameGroup" ~ "DomHelp" ,  (Status.x %in% c("BrM","BrF") & Status.y %in% c("AB","ABX") | Status.y %in% c("BrM","BrF") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "DomSub", (Status.x %in% c("H") & Status.y %in% c("AB","ABX") | Status.y %in% c("H") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "HelpSub", (Status.x %in% c("H","AB","ABX") & Status.y %in% c("H","AB","ABX")) & Pairs %in% "sameGroup"  ~ "SubSub",Pairs %in% "sameGroup" ~ "sameTerrothers", Pairs %in% "sameBird" ~ "sameBird",TRUE ~ Pairs))) %>% mutate(Typeofpair2 = factor(Typeofpair2,level=c("diffGroup","Dominantpair","DomHelp","DomSub","HelpSub","SubSub","sameTerrothers","sameBird")), Pairs=factor(Pairs, level=c("diffGroup","sameGroup","sameBird"))) %>% filter(!(Typeofpair2 %in% "sameTerrothers"))%>% filter(!(Pairs %in% "sameBird"))

physeqAdistance4groupanae <- physeqAdistanceanae4 %>%
  mutate(ItemID1_ = pmin(BirdID.x  ,BirdID.y),
         ItemID2_ = pmax(BirdID.x  ,BirdID.y)) %>%
  group_by(ItemID1_,ItemID2_) %>% mutate(BirdGroupID = cur_group_id())%>% 
  merge(.,relatednessbi,by.x=c("ItemID1_","ItemID2_"),by.y = c("IID1","IID2"),all.x=T) %>% 
  merge(.,firststatusindividuals,by.x = "BirdID.x",by.y="BirdID",all.x=T) %>% 
  merge(.,firststatusindividuals,by.x = "BirdID.y",by.y="BirdID",all.x=T) %>% 
  mutate(samenatal = case_when( (BirdID.y == BrF_1.x | BirdID.y== BrM_1.x | BirdID.y == H_1.x | BirdID.y == BrF_2.x | BirdID.y== BrM_2.x | BirdID.y == H_2.x |BirdID.y == H_3.x | BirdID.x == BrF_1.y | BirdID.x== BrM_1.y | BirdID.y == H_1.y | BirdID.y == BrF_2.y | BirdID.y== BrM_2.y | BirdID.y == H_2.y |BirdID.y == H_3.y ) ~ "same", TRUE ~ "different" )) %>% merge(., Offspring20240530ind,by.x="BirdID.x", by.y = "OffID",all.x=T)


#alpha anaerobes
physeqAdistance4groupanaealpha <- physeqAdistance4groupanae %>% mutate(Observeddiff = abs(Observed.x-Observed.y), Shannondiff = abs(Shannon.x-Shannon.y))
physeqAdistance4groupanaealpha2 <- physeqAdistance4groupanaealpha %>% filter(Pairs %in% "sameGroup")

# beta anaerobes
physeqAdistance4groupanae2 <- physeqAdistance4groupanae %>% filter(Pairs %in% "sameGroup")

