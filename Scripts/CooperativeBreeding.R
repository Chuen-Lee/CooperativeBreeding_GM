library("readxl")
library("data.table")
library("ggdist")
library("lmerMultiMember")
library("emmeans")
library("tidyverse")
library("algatr")
library("ggsignif")


BreedStatus20240827_2 <- read_excel("~/OneDrive - University of East Anglia/Random Tables from SWDB/BreedStatus20240827.xlsx")
BreedGroupLocation20240827 <- read_excel("~/OneDrive - University of East Anglia/Random Tables from SWDB/BreedGroupLocation20240827.xlsx")
BreedStatusWide <- BreedStatus20240827_2 %>% dplyr::select(BreedGroupID,BirdID,Status) %>% merge(.,BreedGroupLocation20240827, by="BreedGroupID") %>% group_by(BirdID,FieldPeriodID) %>% dplyr::arrange(factor(Status, levels = c("BrM","BrF","H","ABX","AB","OFL","FL","CH","BrU","B","TBRF","TBRM","NSA","SEEN2","SEEN1")), .by_group = T) %>% distinct(FieldPeriodID,.keep_all = T) %>% filter(Status %in% c("BrM","BrF","H","CH","ABX","AB","OFL","FL")) %>% ungroup %>% group_by(BreedGroupID) %>% pivot_wider(values_from = "BirdID",names_from = "Status",names_sep = "_",names_vary = "fastest") %>% unnest_wider(BrM,names_sep = "_") %>% unnest_wider(BrF,names_sep = "_") %>% unnest_wider(H,names_sep = "_") %>% unnest_wider(ABX,names_sep = "_") %>% unnest_wider(AB,names_sep = "_") %>% unnest_wider(CH,names_sep = "_") %>% unnest_wider(OFL,names_sep = "_")  %>% unnest_wider(FL,names_sep = "_") %>% ungroup  %>% mutate(n_members = rowSums( !is.na( . [,2:28]))) %>% dplyr::arrange(BreedGroupID) %>% group_by(BrM_1,BrF_1) %>% dplyr::mutate(seasonstgt = case_when(is.na(BrM_1) | is.na(BrF_1) | BrM_1 %in% "0" | BrF_1 %in% "0" ~ 0, TRUE ~ 1:n())) %>% mutate(HelperPresent = case_when(!is.na(H_1) ~ "yes", TRUE ~ "no")) %>% mutate(FledglingPresent = case_when(!is.na(FL_1) | !is.na(OFL_1) ~ "yes", TRUE ~ "no"))

ggplot(BreedStatusWide,aes(n_members)) + geom_histogram()
ggplot(BreedStatusWide,aes(seasonstgt)) + geom_histogram()

BreedStatusWide2 <- BreedStatusWide %>% ungroup() %>% dplyr::select(BreedGroupID,n_members,seasonstgt,HelperPresent,FledglingPresent)

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
physeq3 <- readRDS("~/Documents/PhD/R_analysis/HostGenetics/InputTables/physeq3.rds")
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

physeqAdistance <- physeqAdistance %>% data.frame() %>% rownames_to_column("SampleID") %>% pivot_longer(.,cols=-c(SampleID)) %>% filter(value > 0)
physeq_Abreedclr_st2 <- data.frame(sample_data(PhyseqAbreeddd)) %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID,BreedGroupID,BirdID,season,SampleYear,Status,n_members,seasonstgt,HelperPresent,FledglingPresent,CatchTime,Timeinfridge,SampleDate,SamplingAge,SexEstimate)

physeqAdistance2 <- merge(physeqAdistance, physeq_Abreedclr_st2, by.x=("SampleID"), by.y="SampleID") %>% merge(.,physeq_Abreedclr_st2, by.x="name",by.y="SampleID") %>% merge(.,physeqrarefiedalpha2, by="SampleID",all.x=T) %>% merge(.,physeqrarefiedalpha2, by.x="name",by.y="SampleID",all.x=T) 

physeqAdistance3 <- physeqAdistance2 %>% mutate(Pairs = as.factor(case_when(name == SampleID ~ "samesample" ,as.numeric(BirdID.x) == as.numeric(BirdID.y) ~ "sameBird", as.numeric(BreedGroupID.x) == as.numeric(BreedGroupID.y) ~ "sameGroup", as.numeric(BreedGroupID.x) != as.numeric(BreedGroupID.y) ~ "diffGroup"))) %>% filter(!(Pairs %in% "samesample")) %>% mutate(seasonpairs = as.factor(case_when(season.x == season.y & SampleYear.x == SampleYear.y ~ "sameyearsameseason",season.x == season.y ~ "sameseasondiffyear",SampleYear.x == SampleYear.y ~ "sameyeardiffseason", TRUE ~ "diffyeardiffseason")), noseasonstgt = case_when(Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF") & seasonstgt.x == seasonstgt.y & Pairs %in% "sameGroup" ~ seasonstgt.x)) %>% mutate(SampleYear.x = as.factor(SampleYear.x)) %>% filter(Status.x %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H") & Status.y %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H")) %>% mutate(CatchTimediff = abs(CatchTime.x - CatchTime.y), Timeinfridgediff = abs(Timeinfridge.x - Timeinfridge.y), sampledatediff = abs(as.numeric(as.Date(SampleDate.x,"%d/%m/%Y") - as.Date(SampleDate.y,"%d/%m/%Y"))), agediff = abs(SamplingAge.x-SamplingAge.y), sexdiff = case_when(SexEstimate.x == SexEstimate.y ~ "0", TRUE ~ "1"))

physeqAdistance4 <- physeqAdistance3 %>% filter(seasonpairs %in% "sameyearsameseason")  %>% mutate(Typeofpair2 = as.factor(case_when( (Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF")) & Pairs %in% "sameGroup" ~ "Dominantpair", (Status.x %in% c("BrM","BrF") & Status.y %in% c("H") | Status.y %in% c("BrM","BrF") & Status.x %in% c("H")) & Pairs %in% "sameGroup" ~ "DomHelp" ,  (Status.x %in% c("BrM","BrF") & Status.y %in% c("AB","ABX") | Status.y %in% c("BrM","BrF") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "DomSub", (Status.x %in% c("H") & Status.y %in% c("AB","ABX") | Status.y %in% c("H") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "HelpSub", (Status.x %in% c("H","AB","ABX") & Status.y %in% c("H","AB","ABX")) & Pairs %in% "sameGroup"  ~ "SubSub",Pairs %in% "sameGroup" ~ "sameTerrothers", Pairs %in% "sameBird" ~ "sameBird",TRUE ~ Pairs))) %>% mutate(Typeofpair2 = factor(Typeofpair2,level=c("diffGroup","Dominantpair","DomHelp","DomSub","HelpSub","SubSub","sameTerrothers","sameBird")), Pairs=factor(Pairs, level=c("diffGroup","sameGroup","sameBird"))) %>% filter(!(Typeofpair2 %in% "sameTerrothers"))

physeqAdistance4group <- physeqAdistance4 %>%
  mutate(ItemID1_ = pmin(BirdID.x  ,BirdID.y),
         ItemID2_ = pmax(BirdID.x  ,BirdID.y)) %>%
  group_by(ItemID1_,ItemID2_) %>% mutate(BirdGroupID = cur_group_id())

physeqAdistance4groupsamepair <- physeqAdistance4group %>% filter(Pairs%in%"sameGroup")

physeqAdistance4sss <- physeqAdistance4 %>% group_by(Pairs) %>% dplyr::summarise(count=n()) %>% mutate(value = 125)
physeqAdistance4ss <- physeqAdistance4group %>% filter(Pairs %in% "sameGroup") %>% group_by(Typeofpair2) %>% dplyr::summarise(count=n()) %>% mutate(value = 125) 

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
physeqAdistance4groupalpha <- physeqAdistance4group %>% mutate(Observeddiff = abs(Observed.x-Observed.y), Shannondiff = abs(Shannon.x-Shannon.y))
physeqAdistance4groupalpha2 <- physeqAdistance4groupalpha %>% filter(Pairs %in% "sameGroup")
ggplot(physeqAdistance4groupalpha,aes(Observeddiff)) + geom_histogram()
ggplot(physeqAdistance4groupalpha,aes(Shannondiff)) + geom_histogram()

# observed 
physeqAdistance4groupalphaobslm <- lmer(sqrt(Observeddiff) ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4groupalpha)
simulateResiduals(physeqAdistance4groupalphaobslm,plot=T)
summary(physeqAdistance4groupalphaobslm)
car::Anova(physeqAdistance4groupalphaobslm,type="III")
emmeans(physeqAdistance4groupalphaobslm, pairwise ~ Pairs)
#physeqAdistance4groupalphaobslmdata <- ggpredict(physeqAdistance4groupalphaobslm, terms="Pairs")
#plot(physeqAdistance4groupalphaobslmdata)

physeqAdistance4groupalphaobslm2 <- lmer(sqrt(Observeddiff) ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4groupalpha2)
summary(physeqAdistance4groupalphaobslm2)
car::Anova(physeqAdistance4groupalphaobslm2,type="III")
emmeans(physeqAdistance4groupalphaobslm2, pairwise ~ Typeofpair2)
#physeqAdistance4groupalphaobslm2data <- ggpredict(physeqAdistance4groupalphaobslm2, terms="Typeofpair2")
#plot(physeqAdistance4groupalphaobslm2data)

# shannon
physeqAdistance4groupalphashanlm <- lmer(sqrt(Shannondiff) ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4groupalpha)
summary(physeqAdistance4groupalphashanlm)
simulateResiduals(physeqAdistance4groupalphashanlm,plot=T)
physeqAdistance4groupalphashanlmdata <- ggpredict(physeqAdistance4groupalphashanlm, terms="Pairs")
plot(physeqAdistance4groupalphashanlmdata)

physeqAdistance4groupalphashanlmdata2 <- physeqAdistance4groupalphashanlmdata %>% data.frame() %>% dplyr::rename(value = predicted, Pairs = x) 

ggplot(physeqAdistance4groupalphashanlmdata2, aes(x=Pairs, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4groupalpha, aes(x=Pairs,y=Shannondiff,fill=Pairs),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4sss, aes(x=Pairs, y=6,label=count), inherit.aes = F)+ 
  xlab("Groups") +
  ylab("Shannon diversity difference") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")

physeqAdistance4groupalphashanlm2 <- lmer(sqrt(Shannondiff) ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4groupalpha2)
summary(physeqAdistance4groupalphashanlm2)
car::Anova(physeqAdistance4groupalphashanlm2, type="III")
emmeans(physeqAdistance4groupalphashanlm2, pairwise ~ Typeofpair2)
physeqAdistance4groupalphashanlm2data <- ggpredict(physeqAdistance4groupalphashanlm2, terms="Typeofpair2")
plot(physeqAdistance4groupalphashanlm2data)

physeqAdistance4groupalphashanlm2data2 <- data.frame(physeqAdistance4groupalphashanlm2data) %>% dplyr::rename(value = predicted, Typeofpair2 = x) 

ggplot(physeqAdistance4groupalphashanlm2data2, aes(x=Typeofpair2, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4groupalpha2, aes(x=Typeofpair2,y=Shannondiff,fill=Typeofpair2),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4ss, aes(x=Typeofpair2, y=6,label=count), inherit.aes = F)+ 
  xlab("Groups") +
  ylab("Shannon diversity difference") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")

# paired beta diversity ----
# within group vs between group ----
ggplot(physeqAdistance4group, aes(value))+geom_histogram()

physeqAdistance4grouplm <- lmer(value ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4group)
summary(physeqAdistance4grouplm)
physeqAdistance4grouplmdata <- ggpredict(physeqAdistance4grouplm, terms="Pairs")
plot(physeqAdistance4grouplmdata)

physeqAdistance4grouplmdatadf <- data.frame(physeqAdistance4grouplmdata) %>% dplyr::rename(value = predicted, Pairs = x) 

pairplot1 <- ggplot(physeqAdistance4grouplmdatadf, aes(x=Pairs, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4, aes(x=Pairs,y=value,fill=Pairs),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4sss, aes(x=Pairs, y=value,label=count), inherit.aes = F)+ 
  xlab("Groups") +
  ylab("Microbiome dissimilarity") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")

# type of pair ----
physeqAdistance4group2 <- physeqAdistance4group %>% filter(Pairs %in% "sameGroup")
physeqAdistance4group2 %>% group_by(Typeofpair2) %>% summarise(n())
physeqAdistance4grouplm2 <- lmer(value ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4group2)
summary(physeqAdistance4grouplm2)
car::Anova(physeqAdistance4grouplm2,type="III")
physeqAdistance4grouplm2data <- ggpredict(physeqAdistance4grouplm2, terms="Typeofpair2")
plot(physeqAdistance4grouplm2data)

emmeans(physeqAdistance4grouplm2, pairwise ~ Typeofpair2)

physeqAdistance4grouplm2datadf <- data.frame(physeqAdistance4grouplm2data) %>% dplyr::rename(value = predicted, Typeofpair2 = x) 

typeplot1 <- ggplot(physeqAdistance4grouplm2datadf, aes(x=Typeofpair2, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4group2, aes(x=Typeofpair2,y=value,fill=Typeofpair2),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4ss, aes(x=Typeofpair2, y=value,label=count), inherit.aes = F)+ 
  xlab("Groups") +
  ylab("Microbiome dissimilarity") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")

ggarrange(pairplot1, typeplot1, widths = c(1.2,2))

# genus abundance and aerobic/anaerobic ----
physeq4asnewsdgen <- tax_glom(physeq4asnewsd, taxrank="Genus", NArm=T) 

PhyseqAbreedddcoremembers <- core_members(physeq4asnewsdgen, detection = 0.01, prevalence = 0.5)
PhyseqAbreedddcore <- prune_taxa(PhyseqAbreedddcoremembers, physeq4asnewsdgen)
taxphyseqclr <- tax_table(physeq4asnewsdgen) %>% data.frame() 

taxphyseqclrgenus <- taxphyseqclr %>% distinct(Genus)
taxphyseqclrgenus %>% summarise(n=n())
#write_csv(taxphyseqclrgenus, "taxphyseqclrgenus.csv")
taxaassignchuen <- read_csv("~/Documents/PhD/R_analysis/CooperativeBreeding/taxphyseqclrgenus.csv")
chatgptassign <- read_csv("~/Documents/PhD/R_analysis/CooperativeBreeding/genus_chatgpt.csv")
knowleslabannote <- read_csv("~/Downloads/41559_2024_2381_MOESM3_ESM.csv")

checkchuenchat <- merge(chatgptassign,taxaassignchuen,by="Genus",all=T)%>% filter(!is.na(Aerotolerance) & !is.na(Oxygen)) %>% mutate(conflict = case_when((Oxygen %in% c("Aerobic","Facultative Anaerobic") & Aerotolerance %in% c( "Aerobic","facultatively anaerobic","mixed")) | (Oxygen %in% "Anaerobic" & Aerotolerance %in% "Anaerobe") | Aerotolerance%in% "unknown" ~ "0", TRUE ~ "1"))
checkchuenchat %>% group_by(conflict) %>% summarise(nconflict = n())
47/(17+47) # 73.4%
checkknowleschat <- merge(chatgptassign,knowleslabannote,by="Genus") %>% mutate(conflict = case_when((Oxygen %in% c("Aerobic","Facultatively anaerobic") & Aerotolerance_binary %in% "aerotolerant") | (Oxygen %in% "Anaerobic" & Aerotolerance_binary %in% "anaerobic") | (Oxygen %in% "Unknown" & Aerotolerance_binary %in% "unknown") | Aerotolerance_binary %in% "unknown" ~ "0", TRUE ~ "1"))
checkknowleschat %>% group_by(conflict) %>% summarise(nconflict = n())
98/(98+34) #74.2%

geminiassign <- read_csv("~/Documents/PhD/R_analysis/CooperativeBreeding/genus_gemini.csv")
geminiassign %>% group_by(Aero) %>% summarise(n=n()) # Aerobic 651, Anaerobic 224, Facultative Anaerobic 65, Not applicable 2, Unknown 133, Variable 36

checkchuengemini <- merge(geminiassign,taxaassignchuen,by="Genus",all=T) %>% filter(!is.na(Aerotolerance)) %>% mutate(conflict = case_when((Aero %in% c("Aerobic","Facultative Anaerobic") & Aerotolerance %in% c( "Aerobic","facultatively anaerobic","mixed")) | (Aero %in% "Anaerobic" & Aerotolerance %in% "Anaerobe") | Aerotolerance%in% "unknown" ~ "0", TRUE ~ "1"))
checkknowlesgemini <- merge(geminiassign,knowleslabannote,by="Genus") %>% mutate(conflict = case_when((Aero %in% c("Aerobic","Facultative Anaerobic") & Aerotolerance_binary %in% "aerotolerant") | (Aero %in% "Anaerobic" & Aerotolerance_binary %in% "anaerobic") | (Aero %in% "Unknown" & Aerotolerance_binary %in% "unknown") | Aerotolerance_binary %in% "unknown" ~ "0", TRUE ~ "1"))

taxgenusassign <- merge(taxphyseqclrgenus,taxaassignchuen, by="Genus",all.x=T )
taxgenusassign %>% filter(is.na(Aerotolerance)) %>% nrow()

checkchuengemini %>% group_by(conflict) %>% summarise(nconflict = n())
77/(77+3) # 96.3% accuracy in 80 genera
checkknowlesgemini %>% group_by(conflict) %>% summarise(nconflict = n())
160/(160+13) # 92.5% accuracy in 173 genera
161/(161+10)

taxphyseqclrgenusannote <- merge(taxphyseqclrgenus,knowleslabannote, by="Genus",all.x=T)

taxphyseqclrgenusannote %>% filter(is.na(Aerotolerance_binary)) %>% nrow()

taxgenusassign2 <- merge(taxgenusassign,taxphyseqclrgenusannote,by="Genus",all=T) %>% mutate(oxy = case_when( (Aerotolerance %in% c("Aerobic","facultatively anaerobic") & Aerotolerance_binary %in% "aerotolerant") | (Aerotolerance %in% c("Aerobic","facultatively anaerobic") & is.na(Aerotolerance_binary)) | (is.na(Aerotolerance) & Aerotolerance_binary %in% "aerotolerant") ~ "aerotolerant", (Aerotolerance %in% c("Anaerobe") & Aerotolerance_binary %in% "anaerobic") | (Aerotolerance %in% c("Anaerobe") & is.na(Aerotolerance_binary)) | (is.na(Aerotolerance) & Aerotolerance_binary %in% "anaerobic") ~ "anaerobic"  ))

#write.csv(taxgenusassign2, "~/Downloads/taxgenusassign2.csv")

taxgenusassign3 <- read_csv("~/Downloads/taxgenusassign2.csv") %>% mutate(spore= case_when((Spore_formation.x %in% "Y" & Spore_formation.y %in% "SF") | (is.na(Spore_formation.x) & Spore_formation.y %in% "SF") | (Spore_formation.x %in% "Y" & is.na(Spore_formation.y))  ~ "Y", (Spore_formation.x %in% "N" & Spore_formation.y %in% "NSF") | (is.na(Spore_formation.x) & Spore_formation.y %in% "NSF") | (Spore_formation.x %in% "N" & is.na(Spore_formation.y))  ~ "N")) %>% dplyr::select(Genus, spore, oxy)

taxaphyseqclr <- transform(physeq4asnewsdgen,"clr") %>% psmelt() %>% merge(.,geminiassign, by="Genus", all.x=T ) 

taxaphyseqclr %>% distinct(Genus) %>% nrow()
taxaphyseqclr %>% distinct(Aero) 
taxaphyseqclr %>% distinct(Genus,.keep_all = T) %>% group_by(Aero) %>% summarise(n())

taxaphyseqclrsum <- taxaphyseqclr%>% group_by(Sample,Aero) %>% summarise(sumgenusoxy = sum(Abundance))
ggplot(taxaphyseqclrsum, aes(x=Aero,y=sumgenusoxy)) + geom_boxplot()

physeqrareaero <- transform(physeq4asnewsdgen, "compositional") %>% core_members(.,detection = 0.0001, prevalence= 0.05)
physeq4asnewsdgenclr <-  transform(physeq4asnewsdgen, "clr")
# aerobic only  ----
aerobicasv <- taxaphyseqclr %>% filter(Aero %in% c("Aerobic","Facultative Anaerobic")) %>% distinct(OTU)
taxaphyseqclraerobic <- prune_taxa(aerobicasv$OTU, physeq4asnewsdgenclr)
# paired distances -
physeqAdistanceaero <- phyloseq::distance(taxaphyseqclraerobic, method = "euclidean", type="samples") %>% as.matrix() 
physeqAdistanceaero[lower.tri(physeqAdistanceaero)] <- 0

physeqAdistanceaero <- physeqAdistanceaero %>% data.frame() %>% rownames_to_column("SampleID") %>% pivot_longer(.,cols=-c(SampleID)) %>% filter(value > 0)
physeq_Abreedclr_st2aero <- data.frame(sample_data(PhyseqAbreeddd)) %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID,BreedGroupID,BirdID,season,SampleYear,Status,n_members,seasonstgt,HelperPresent,FledglingPresent,CatchTime,Timeinfridge,SampleDate,SamplingAge,SexEstimate)

physeqAdistanceaero2 <- merge(physeqAdistanceaero, physeq_Abreedclr_st2aero, by.x=("SampleID"), by.y="SampleID") %>% merge(.,physeq_Abreedclr_st2, by.x="name",by.y="SampleID") %>% merge(.,physeqrarefiedalpha2, by="SampleID",all.x=T) %>% merge(.,physeqrarefiedalpha2, by.x="name",by.y="SampleID",all.x=T) 

physeqAdistanceaero3 <- physeqAdistanceaero2 %>% mutate(Pairs = as.factor(case_when(name == SampleID ~ "samesample" ,as.numeric(BirdID.x) == as.numeric(BirdID.y) ~ "sameBird", as.numeric(BreedGroupID.x) == as.numeric(BreedGroupID.y) ~ "sameGroup", as.numeric(BreedGroupID.x) != as.numeric(BreedGroupID.y) ~ "diffGroup"))) %>% filter(!(Pairs %in% "samesample")) %>% mutate(seasonpairs = as.factor(case_when(season.x == season.y & SampleYear.x == SampleYear.y ~ "sameyearsameseason",season.x == season.y ~ "sameseasondiffyear",SampleYear.x == SampleYear.y ~ "sameyeardiffseason", TRUE ~ "diffyeardiffseason")), noseasonstgt = case_when(Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF") & seasonstgt.x == seasonstgt.y & Pairs %in% "sameGroup" ~ seasonstgt.x)) %>% mutate(SampleYear.x = as.factor(SampleYear.x)) %>% filter(Status.x %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H") & Status.y %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H")) %>% mutate(CatchTimediff = abs(CatchTime.x - CatchTime.y), Timeinfridgediff = abs(Timeinfridge.x - Timeinfridge.y), sampledatediff = abs(as.numeric(as.Date(SampleDate.x,"%d/%m/%Y") - as.Date(SampleDate.y,"%d/%m/%Y"))), agediff = SamplingAge.x-SamplingAge.y, sexdiff = case_when(SexEstimate.x == SexEstimate.y ~ "0", TRUE ~ "1"))

physeqAdistanceaero4 <- physeqAdistanceaero3 %>% filter(seasonpairs %in% "sameyearsameseason")  %>% mutate(Typeofpair2 = as.factor(case_when( (Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF")) & Pairs %in% "sameGroup" ~ "Dominantpair", (Status.x %in% c("BrM","BrF") & Status.y %in% c("H") | Status.y %in% c("BrM","BrF") & Status.x %in% c("H")) & Pairs %in% "sameGroup" ~ "DomHelp" ,  (Status.x %in% c("BrM","BrF") & Status.y %in% c("AB","ABX") | Status.y %in% c("BrM","BrF") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "DomSub", (Status.x %in% c("H") & Status.y %in% c("AB","ABX") | Status.y %in% c("H") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "HelpSub", (Status.x %in% c("H","AB","ABX") & Status.y %in% c("H","AB","ABX")) & Pairs %in% "sameGroup"  ~ "SubSub",Pairs %in% "sameGroup" ~ "sameTerrothers", Pairs %in% "sameBird" ~ "sameBird",TRUE ~ Pairs))) %>% mutate(Typeofpair2 = factor(Typeofpair2,level=c("diffGroup","Dominantpair","DomHelp","DomSub","HelpSub","SubSub","sameTerrothers","sameBird")), Pairs=factor(Pairs, level=c("diffGroup","sameGroup","sameBird"))) %>% filter(!(Typeofpair2 %in% "sameTerrothers"))

physeqAdistance4groupaero <- physeqAdistanceaero4 %>%
  mutate(ItemID1_ = pmin(BirdID.x  ,BirdID.y),
         ItemID2_ = pmax(BirdID.x  ,BirdID.y)) %>%
  group_by(ItemID1_,ItemID2_) %>% mutate(BirdGroupID = cur_group_id())


physeqAdistance4groupaerolm <- lmer(value ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4groupaero)
summary(physeqAdistance4groupaerolm)
physeqAdistance4groupaerolmdata <- ggpredict(physeqAdistance4groupaerolm,terms = "Pairs")
plot(physeqAdistance4groupaerolmdata)

physeqAdistance4groupaero2 <- physeqAdistance4groupaero %>% filter(Pairs %in% "sameGroup")
physeqAdistance4groupaerolm2 <- lmer(value ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4groupaero2)
summary(physeqAdistance4groupaerolm2)
car::Anova(physeqAdistance4groupaerolm2, type="III")
physeqAdistance4groupaerolm2data <- ggpredict(physeqAdistance4groupaerolm2,terms = "Typeofpair2")
plot(physeqAdistance4groupaerolm2data)

physeqAdistance4groupaerolm2data2 <- data.frame(physeqAdistance4groupaerolm2data) %>% dplyr::rename(value = predicted, Typeofpair2 = x) 

aerotypeplot1 <- ggplot(physeqAdistance4groupaerolm2data2, aes(x=Typeofpair2, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4groupaero2, aes(x=Typeofpair2,y=value,fill=Typeofpair2),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4ss, aes(x=Typeofpair2, y=80,label=count), inherit.aes = F)+ 
  xlab("Groups") +
  ylab("Aerobic microbiome dissimilarity") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")

# anaerobes only ----
anaasv <- taxaphyseqclr %>% filter(Aero %in% c("Anaerobic")) %>% distinct(OTU)
taxaphyseqclranae <- prune_taxa(anaasv$OTU, physeq4asnewsdgenclr)

# paired distances -
physeqAdistanceanae <- phyloseq::distance(taxaphyseqclranae, method = "euclidean", type="samples") %>% as.matrix() 
physeqAdistanceanae[lower.tri(physeqAdistanceanae)] <- 0

physeqAdistanceanae <- physeqAdistanceanae %>% data.frame() %>% rownames_to_column("SampleID") %>% pivot_longer(.,cols=-c(SampleID)) %>% filter(value > 0)
physeq_Abreedclr_st2anae <- data.frame(sample_data(PhyseqAbreeddd)) %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID,BreedGroupID,BirdID,season,SampleYear,Status,n_members,seasonstgt,HelperPresent,FledglingPresent,CatchTime,Timeinfridge,SampleDate,SamplingAge,SexEstimate)

physeqAdistanceanae2 <- merge(physeqAdistanceanae, physeq_Abreedclr_st2anae, by.x=("SampleID"), by.y="SampleID") %>% merge(.,physeq_Abreedclr_st2, by.x="name",by.y="SampleID") %>% merge(.,physeqrarefiedalpha2, by="SampleID",all.x=T) %>% merge(.,physeqrarefiedalpha2, by.x="name",by.y="SampleID",all.x=T) 

physeqAdistanceanae3 <- physeqAdistanceanae2 %>% mutate(Pairs = as.factor(case_when(name == SampleID ~ "samesample" ,as.numeric(BirdID.x) == as.numeric(BirdID.y) ~ "sameBird", as.numeric(BreedGroupID.x) == as.numeric(BreedGroupID.y) ~ "sameGroup", as.numeric(BreedGroupID.x) != as.numeric(BreedGroupID.y) ~ "diffGroup"))) %>% filter(!(Pairs %in% "samesample")) %>% mutate(seasonpairs = as.factor(case_when(season.x == season.y & SampleYear.x == SampleYear.y ~ "sameyearsameseason",season.x == season.y ~ "sameseasondiffyear",SampleYear.x == SampleYear.y ~ "sameyeardiffseason", TRUE ~ "diffyeardiffseason")), noseasonstgt = case_when(Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF") & seasonstgt.x == seasonstgt.y & Pairs %in% "sameGroup" ~ seasonstgt.x)) %>% mutate(SampleYear.x = as.factor(SampleYear.x)) %>% filter(Status.x %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H") & Status.y %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H")) %>% mutate(CatchTimediff = abs(CatchTime.x - CatchTime.y), Timeinfridgediff = abs(Timeinfridge.x - Timeinfridge.y), sampledatediff = abs(as.numeric(as.Date(SampleDate.x,"%d/%m/%Y") - as.Date(SampleDate.y,"%d/%m/%Y"))), agediff = SamplingAge.x-SamplingAge.y, sexdiff = case_when(SexEstimate.x == SexEstimate.y ~ "0", TRUE ~ "1"))

physeqAdistanceanae4 <- physeqAdistanceanae3 %>% filter(seasonpairs %in% "sameyearsameseason")  %>% mutate(Typeofpair2 = as.factor(case_when( (Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF")) & Pairs %in% "sameGroup" ~ "Dominantpair", (Status.x %in% c("BrM","BrF") & Status.y %in% c("H") | Status.y %in% c("BrM","BrF") & Status.x %in% c("H")) & Pairs %in% "sameGroup" ~ "DomHelp" ,  (Status.x %in% c("BrM","BrF") & Status.y %in% c("AB","ABX") | Status.y %in% c("BrM","BrF") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "DomSub", (Status.x %in% c("H") & Status.y %in% c("AB","ABX") | Status.y %in% c("H") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "HelpSub", (Status.x %in% c("H","AB","ABX") & Status.y %in% c("H","AB","ABX")) & Pairs %in% "sameGroup"  ~ "SubSub",Pairs %in% "sameGroup" ~ "sameTerrothers", Pairs %in% "sameBird" ~ "sameBird",TRUE ~ Pairs))) %>% mutate(Typeofpair2 = factor(Typeofpair2,level=c("diffGroup","Dominantpair","DomHelp","DomSub","HelpSub","SubSub","sameTerrothers","sameBird")), Pairs=factor(Pairs, level=c("diffGroup","sameGroup","sameBird"))) %>% filter(!(Typeofpair2 %in% "sameTerrothers"))

physeqAdistance4groupanae <- physeqAdistanceanae4 %>%
  mutate(ItemID1_ = pmin(BirdID.x  ,BirdID.y),
         ItemID2_ = pmax(BirdID.x  ,BirdID.y)) %>%
  group_by(ItemID1_,ItemID2_) %>% mutate(BirdGroupID = cur_group_id())


physeqAdistance4groupanaelm <- lmer(value ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4groupanae)
summary(physeqAdistance4groupanaelm)
physeqAdistance4groupanaelmdata <- ggpredict(physeqAdistance4groupanaelm,terms = "Pairs")
plot(physeqAdistance4groupanaelmdata)

physeqAdistance4groupanae2 <- physeqAdistance4groupanae %>% filter(Pairs %in% "sameGroup")
physeqAdistance4groupanaelm2 <- lmer(value ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + (1 | BirdGroupID) + (1|SampleYear.x), data=physeqAdistance4groupanae2)
summary(physeqAdistance4groupanaelm2)
car::Anova(physeqAdistance4groupanaelm2, type="III")
emmeans(physeqAdistance4groupanaelm2, pairwise~Typeofpair2)
physeqAdistance4groupanaelm2data <- ggpredict(physeqAdistance4groupanaelm2,terms = "Typeofpair2")
plot(physeqAdistance4groupanaelm2data)
plot(physeqAdistance4groupanaelm2data,show_data=T)

physeqAdistance4groupanaelm2data2 <- data.frame(physeqAdistance4groupanaelm2data) %>% dplyr::rename(value = predicted, Typeofpair2 = x) 

anaerotypeplot1 <-ggplot(physeqAdistance4groupanaelm2data2, aes(x=Typeofpair2, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4groupanae2, aes(x=Typeofpair2,y=value,fill=Typeofpair2),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4ss, aes(x=Typeofpair2, y=50,label=count), inherit.aes = F)+ geom_signif(comparisons = list(c("Dominantpair","DomSub")),annotations="*") +
  xlab("Groups") +
  ylab("Anaerobic microbiome dissimilarity") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")

ggarrange(aerotypeplot1,anaerotypeplot1)

# facultative anaerobes only
fanaasv <- taxaphyseqclr %>% filter(Aero %in% "Facultative Anaerobic") %>% distinct(OTU)
taxaphyseqclrfanae <- prune_taxa(fanaasv$OTU, PhyseqAbreeddd)

# not enough metagnenomic samples ----
#load sample table (needs editing)
st <- read.csv("~/Documents/PhD/R_analysis/SampleMetadata/InputTables/st.csv")
# sample data 
st$TubeNumber <- sub("^","SW",st$TubeNumber)
# add TubeNumber as rownames 
stA <- st %>% filter(!duplicated(TubeNumber)) %>% mutate(deathdate = LatestEligibleDate) %>% mutate(SampleYear = format(as.Date(SampleDate, format="%Y-%m-%d"), "%Y")) 
st3<-as.data.frame(stA)
st4 <- st3[,-1]
rownames(st4) <- st3[,1]
stdata <- sample_data(st4)
st4 #from taxo

# mpa ----
mpamat <- read.csv("~/Documents/PhD/R_analysis/Phyloseq/InputTables/MFF_Feb2024/merged_reads_abundance.txt", sep = "\t") %>% pivot_wider(names_from = TubeNo, values_from = Reads) %>% dplyr::mutate(OTU = paste0("OTU", 1:nrow(.))) %>% dplyr::rename(SW370=SW370_, SW29=SW29_,SW131=SW131_,SW55=SW55_)
# taxa table
mpataxa <- mpamat %>% dplyr::select(OTU,OUT) %>% separate_wider_delim(OUT, delim = "|", names = c("Kingdom","Phylum","Order","Class","Family","Genus","Species","SGB"), too_few = "align_start") %>% column_to_rownames(., var="OTU")
mpataxo7 <- mpataxa %>% filter(!is.na(SGB))
# otu table
mpaL7 <- mpamat  %>% dplyr::select(-OUT) %>% column_to_rownames(., var="OTU") %>% replace(is.na(.), 0)
mpaL7 <- mpaL7[rownames(mpaL7) %in% row.names(mpataxo7), ]

#matrix and phyloseq
mpaOTUL7 <- otu_table(as.matrix(mpaL7), taxa_are_rows = TRUE)

mpaTAXL7 = tax_table(as.matrix(mpataxo7))

mpapL7 <- phyloseq(mpaOTUL7, mpaTAXL7, stdata) 
#mpapL7 %>% filter(Type %in% "C") %>% distinct(Sample)
mpapL7sd <- data.frame(sample_data(mpapL7)) %>% rownames_to_column("TubeNumber") %>% mutate(negpos = case_when(TubeNumber %in% c("SW2021","SW2022","SW2017","SW2019","SW2020","SW1421") ~ TRUE, TRUE ~ FALSE)) %>% column_to_rownames("TubeNumber")

# paired distances ----
# mpaAdistance <- phyloseq::distance(mpaAbreeddd, method = "euclidean", type="samples") %>% as.matrix() 
# mpaAdistance[lower.tri(mpaAdistance)] <- 0
# 
# mpaAdistance <- mpaAdistance %>% data.frame() %>% rownames_to_column("SampleID") %>% pivot_longer(.,cols=-c(SampleID)) %>% filter(value > 0)
# mpa_Abreedclr_st2 <- data.frame(sample_data(mpaAbreeddd)) %>% rownames_to_column("SampleID") %>% dplyr::select(SampleID,BreedGroupID,BirdID,season,SampleYear,Status,n_members,seasonstgt,HelperPresent,FledglingPresent,CatchTime,Timeinfridge,SampleDate,SamplingAge,SexEstimate)
# 
# mpaAdistance2 <- merge(mpaAdistance, mpa_Abreedclr_st2, by.x=("SampleID"), by.y="SampleID") %>% merge(.,mpa_Abreedclr_st2, by.x="name",by.y="SampleID") %>% merge(.,mpararefiedalpha2, by="SampleID",all.x=T) %>% merge(.,mpararefiedalpha2, by.x="name",by.y="SampleID",all.x=T) 
# 
# mpaAdistance3 <- mpaAdistance2 %>% mutate(Pairs = as.factor(case_when(name == SampleID ~ "samesample" ,as.numeric(BirdID.x) == as.numeric(BirdID.y) ~ "sameBird", as.numeric(BreedGroupID.x) == as.numeric(BreedGroupID.y) ~ "sameGroup", as.numeric(BreedGroupID.x) != as.numeric(BreedGroupID.y) ~ "diffGroup"))) %>% filter(!(Pairs %in% "samesample")) %>% mutate(seasonpairs = as.factor(case_when(season.x == season.y & SampleYear.x == SampleYear.y ~ "sameyearsameseason",season.x == season.y ~ "sameseasondiffyear",SampleYear.x == SampleYear.y ~ "sameyeardiffseason", TRUE ~ "diffyeardiffseason")), noseasonstgt = case_when(Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF") & seasonstgt.x == seasonstgt.y & Pairs %in% "sameGroup" ~ seasonstgt.x)) %>% mutate(SampleYear.x = as.factor(SampleYear.x)) %>% filter(Status.x %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H") & Status.y %in% c("BrM","BrF","H","AB","ABX","OFL","FL","H")) %>% mutate(CatchTimediff = abs(CatchTime.x - CatchTime.y), Timeinfridgediff = abs(Timeinfridge.x - Timeinfridge.y), sampledatediff = abs(as.numeric(as.Date(SampleDate.x,"%d/%m/%Y") - as.Date(SampleDate.y,"%d/%m/%Y"))), agediff = SamplingAge.x-SamplingAge.y, sexdiff = case_when(SexEstimate.x == SexEstimate.y ~ "0", TRUE ~ "1"))
# 
# mpaAdistance4 <- mpaAdistance3 %>% filter(seasonpairs %in% "sameyearsameseason")  %>% mutate(Typeofpair2 = as.factor(case_when( (Status.x %in% c("BrM","BrF") & Status.y %in% c("BrM","BrF")) & Pairs %in% "sameGroup" ~ "Dominantpair", (Status.x %in% c("BrM","BrF") & Status.y %in% c("H") | Status.y %in% c("BrM","BrF") & Status.x %in% c("H")) & Pairs %in% "sameGroup" ~ "DomHelp" ,  (Status.x %in% c("BrM","BrF") & Status.y %in% c("AB","ABX") | Status.y %in% c("BrM","BrF") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "DomSub", (Status.x %in% c("H") & Status.y %in% c("AB","ABX") | Status.y %in% c("H") & Status.x %in% c("AB","ABX")) & Pairs %in% "sameGroup" ~ "HelpSub", (Status.x %in% c("H","AB","ABX") & Status.y %in% c("H","AB","ABX")) & Pairs %in% "sameGroup"  ~ "SubSub",Pairs %in% "sameGroup" ~ "sameTerrothers", Pairs %in% "sameBird" ~ "sameBird",TRUE ~ Pairs))) %>% mutate(Typeofpair2 = factor(Typeofpair2,level=c("diffGroup","Dominantpair","DomHelp","DomSub","HelpSub","SubSub","sameTerrothers","sameBird")), Pairs=factor(Pairs, level=c("diffGroup","sameGroup","sameBird"))) %>% filter(!(Typeofpair2 %in% "sameTerrothers"))
# 
# mpaAdistance4group <- mpaAdistance4 %>%
#   mutate(ItemID1_ = pmin(BirdID.x  ,BirdID.y),
#          ItemID2_ = pmax(BirdID.x  ,BirdID.y)) %>%
#   group_by(ItemID1_,ItemID2_) %>% mutate(BirdGroupID = cur_group_id())