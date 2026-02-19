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
library("gtsummary")


# multimembership models ----
library(lmerMultiMember)

physeqAdistance4groupbirds <- physeqAdistance4group %>% dplyr::select(BirdID.x,BirdID.y)
physeqAdistance4groupmembers <- weights_from_columns(physeqAdistance4groupbirds)

physeqAdistance4groupsamepairbirds <- physeqAdistance4groupsamepair %>% dplyr::select(BirdID.x,BirdID.y)
physeqAdistance4groupsamepairmembers <- weights_from_columns(physeqAdistance4groupsamepairbirds)

classmodel <- lmer(value ~  + R.ped + samenatal+ (1 | BirdGroupID) + (1|SampleYear.x) , data=physeqAdistance4groupalpha )

# alpha models ----
## obs pairs
physeqAdistance4groupalphaobslmmulti <- lmerMultiMember::lmer(Observeddiff ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + R.ped + samenatal + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupalpha, memberships = list(Birdmultimembers = physeqAdistance4groupmembers))
class(physeqAdistance4groupalphaobslmmulti) <- class(classmodel)
summary(physeqAdistance4groupalphaobslmmulti)
tbl_regression(physeqAdistance4groupalphaobslmmulti, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff","Pairs", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 
vif(physeqAdistance4groupalphaobslmmulti)
# shannon pairs
ggplot(physeqAdistance4groupalpha, aes(Shannondiff)) + geom_histogram()
physeqAdistance4groupalphashanlmmulti <- lmerMultiMember::lmer(Shannondiff ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + R.ped + samenatal + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupalpha, memberships = list(Birdmultimembers = physeqAdistance4groupmembers))
class(physeqAdistance4groupalphashanlmmulti) <- class(classmodel)
summary(physeqAdistance4groupalphashanlmmulti)
tbl_regression(physeqAdistance4groupalphashanlmmulti, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff","Pairs", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 

#obs type
physeqAdistance4groupalphaobslm2multi <- lmerMultiMember::lmer(Observeddiff ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff + R.ped + samenatal + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupalpha2, memberships = list(Birdmultimembers = physeqAdistance4groupsamepairmembers))
class(physeqAdistance4groupalphaobslm2multi) <- class(classmodel)
summary(physeqAdistance4groupalphaobslm2multi)
tbl_regression(physeqAdistance4groupalphaobslm2multi, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 

# shannon type
physeqAdistance4groupalphashanlm2multi <- lmerMultiMember::lmer(Shannondiff ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff + R.ped + samenatal + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupalpha2, memberships = list(Birdmultimembers = physeqAdistance4groupsamepairmembers))
class(physeqAdistance4groupalphashanlm2multi) <- class(classmodel)
summary(physeqAdistance4groupalphashanlm2multi)
tbl_regression(physeqAdistance4groupalphashanlm2multi, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 

emmeans(physeqAdistance4groupalphashanlm2multi, pairwise ~ Typeofpair2)

# beta models ----
## beta pairs ----
physeqAdistance4grouplmmulti <- lmerMultiMember::lmer(value ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff + R.ped + samenatal  + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4group, memberships = list(Birdmultimembers = physeqAdistance4groupmembers)) 
class(physeqAdistance4grouplmmulti) <- class(classmodel)
summary(physeqAdistance4grouplmmulti)
vif(physeqAdistance4grouplmmulti)
emmeans(physeqAdistance4grouplmmulti, pairwise ~ Pairs, adjust = "none")
tbl_regression(physeqAdistance4grouplmmulti, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff","Pairs", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 

## beta type ----
physeqAdistance4grouplm2multi <- lmerMultiMember::lmer(value ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff  + R.ped + samenatal + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupsamepair, memberships = list(Birdmultimembers = physeqAdistance4groupsamepairmembers))
class(physeqAdistance4grouplm2multi) <- class(classmodel)
summary(physeqAdistance4grouplm2multi)
emmeans(physeqAdistance4grouplm2multi, pairwise ~ Typeofpair2, adjust = "none")
tbl_regression(physeqAdistance4grouplm2multi, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 


# aerobic model ----
# pairs
physeqAdistance4grouppairsbirdsaero <- physeqAdistance4groupaero %>% dplyr::select(BirdID.x,BirdID.y)
physeqAdistance4grouppairsbirdsmembersaero <- weights_from_columns(physeqAdistance4grouppairsbirdsaero)
physeqAdistance4groupaerolmmulti <- lmerMultiMember::lmer(value ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff + R.ped + samenatal  + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupaero, memberships = list(Birdmultimembers = physeqAdistance4grouppairsbirdsmembersaero))
class(physeqAdistance4groupaerolmmulti) <- class(classmodel)
summary(physeqAdistance4groupaerolmmulti)
tbl_regression(physeqAdistance4groupaerolmmulti, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff","Pairs", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 

# type
physeqAdistance4groupsamepairbirdsanae <- physeqAdistance4groupanae2 %>% dplyr::select(BirdID.x,BirdID.y)
physeqAdistance4groupsamepairmembersanae <- weights_from_columns(physeqAdistance4groupsamepairbirdsanae)

physeqAdistance4groupsamepairbirdsaero <- physeqAdistance4groupaero2 %>% dplyr::select(BirdID.x,BirdID.y)
physeqAdistance4groupsamepairmembersaero <- weights_from_columns(physeqAdistance4groupsamepairbirdsaero)

physeqAdistance4groupaerolm2multi <-  lmerMultiMember::lmer(value ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff + R.ped + samenatal + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupaero2, memberships = list(Birdmultimembers = physeqAdistance4groupsamepairmembersaero))
class(physeqAdistance4groupaerolm2multi) <- class(classmodel)
summary(physeqAdistance4groupaerolm2multi)
emmeans(physeqAdistance4groupaerolm2multi, pairwise ~ Typeofpair2)
tbl_regression(physeqAdistance4groupaerolm2multi, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 
car::Anova(physeqAdistance4groupaerolm2multi, type="3")


# anaerobic model ----
#pairs
physeqAdistance4grouppairsbirdsanae <- physeqAdistance4groupanae %>% dplyr::select(BirdID.x,BirdID.y)
physeqAdistance4grouppairsbirdsmembersanae <- weights_from_columns(physeqAdistance4grouppairsbirdsanae)

physeqAdistance4groupanaelmmulti <- lmerMultiMember::lmer(value ~ Pairs + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff + R.ped + samenatal  + (1 |Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupanae, memberships = list(Birdmultimembers = physeqAdistance4grouppairsbirdsmembersanae))
class(physeqAdistance4groupanaelmmulti) <- class(classmodel)
summary(physeqAdistance4groupanaelmmulti)
tbl_regression(physeqAdistance4groupanaelmmulti, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff","Pairs", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 

physeqAdistance4groupsamepairbirdsanae <- physeqAdistance4groupanae2 %>% dplyr::select(BirdID.x,BirdID.y)
physeqAdistance4groupsamepairmembersanae <- weights_from_columns(physeqAdistance4groupsamepairbirdsanae)

#type
physeqAdistance4groupanaelm2multi <- lmerMultiMember::lmer(value ~ Typeofpair2 + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff + R.ped + samenatal  + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupanae2, memberships = list(Birdmultimembers = physeqAdistance4groupsamepairmembersanae))
class(physeqAdistance4groupanaelm2multi) <- class(classmodel)
summary(physeqAdistance4groupanaelm2multi)
emmeans(physeqAdistance4groupanaelm2multi, pairwise ~ Typeofpair2)
# Dominantpair - DomSub   -2.8284 1.34 214  -2.107  0.0362
tbl_regression(physeqAdistance4groupanaelm2multi, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff", "samenatal")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 
car::Anova(physeqAdistance4groupanaelm2multi,type="3")

# interacting pairs vs not interacting pairs 
physeqAdistance4groupanae3 <- physeqAdistance4groupanae2 %>% mutate(grouppair = case_when(Typeofpair2 %in% c("Dominantpair","DomHelp") ~ "interacting", TRUE ~ "not"))
physeqAdistance4groupanae3pairbirdsanae <- physeqAdistance4groupanae3 %>% dplyr::select(BirdID.x,BirdID.y)
physeqAdistance4groupsame3pairmembersanae <- weights_from_columns(physeqAdistance4groupanae3pairbirdsanae)

physeqAdistance4groupanaelm2multi3 <- lmerMultiMember::lmer(value ~ grouppair + agediff + sexdiff + season.x + CatchTimediff + Timeinfridgediff + R.ped + samenatal  + (1 | Birdmultimembers) + (1|SampleYear.x), data=physeqAdistance4groupanae3, memberships = list(Birdmultimembers = physeqAdistance4groupsame3pairmembersanae))
class(physeqAdistance4groupanaelm2multi3) <- class(classmodel)
summary(physeqAdistance4groupanaelm2multi3)
tbl_regression(physeqAdistance4groupanaelm2multi3, intercept = T, estimate_fun = purrr::partial(style_ratio, digits = 3),pvalue_fun = purrr::partial(style_sigfig, digits = 3),tidy_fun = broom.mixed::tidy, show_single_row=c("season.x","sexdiff", "samenatal", "grouppair")) %>%
  modify_column_hide(column = conf.low) %>%
  modify_column_unhide(column = c(std.error,df,statistic)) %>%
  add_glance_table(include = nobs) 

n_distinct(rbind(data.frame(idx=physeqAdistance4groupanae3$BirdID.x),data.frame(idx=physeqAdistance4groupanae3$BirdID.y)))
n_distinct(rbind(data.frame(idx=physeqAdistance4groupanae3$name),data.frame(idx=physeqAdistance4groupanae3$SampleID)))

# Plots ----
# alpha diversity
physeqAdistance4groupalphashanlm2multidata <- ggpredict(physeqAdistance4groupalphashanlm2multi, terms="Typeofpair2")

physeqAdistance4groupalphashanlm2multidata2 <- data.frame(physeqAdistance4groupalphashanlm2multidata) %>% dplyr::rename(value = predicted, Typeofpair2 = x) 

shannontype2 <- ggplot(physeqAdistance4groupalphashanlm2multidata2, aes(x=Typeofpair2, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4groupalpha2, aes(x=Typeofpair2,y=Shannondiff,fill=Typeofpair2),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4ss, aes(x=Typeofpair2, y=0.8,label=count), inherit.aes = F) + 
  geom_signif(comparisons = list(c("Dominantpair","DomHelp")),annotations= c("*"),y_position = c(0)) + scale_x_discrete(labels= c("Dom-Dom","Dom-Help","Dom-Sub","Help-Sub","Sub-Sub"))+ 
  xlab("Groups") +
  ylab("Shannon diversity similarity") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")

# pair beta
physeqAdistance4grouplmmultidata <- ggpredict(physeqAdistance4grouplmmulti, terms="Pairs")

physeqAdistance4grouplmmultidatadf <- data.frame(physeqAdistance4grouplmmultidata) %>% dplyr::rename(value = predicted, Pairs = x) 

pairplot2 <- ggplot(physeqAdistance4grouplmmultidatadf, aes(x=Pairs, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4, aes(x=Pairs,y=value,fill=Pairs),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4sss, aes(x=Pairs, y=-20,label=count), inherit.aes = F) + 
  geom_signif(comparisons = list(c("diffGroup","sameGroup")),annotations= c("***"),y_position = c(-35))  + scale_x_discrete(labels= c("Between","Within"))+ 
  xlab("Groups") +
  ylab("Microbiome similarity") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")

# beta type
physeqAdistance4grouplm2multidata <- ggpredict(physeqAdistance4grouplm2multi, terms="Typeofpair2")
physeqAdistance4grouplm2multidatadf <- data.frame(physeqAdistance4grouplm2multidata) %>% dplyr::rename(value = predicted, Typeofpair2 = x) 

typeplot2 <- ggplot(physeqAdistance4grouplm2multidatadf, aes(x=Typeofpair2, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4groupsamepair, aes(x=Typeofpair2,y=value,fill=Typeofpair2),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4ss, aes(x=Typeofpair2, y=-30,label=count), inherit.aes = F)+ scale_x_discrete(labels= c("Dom-Dom","Dom-Help","Dom-Sub","Help-Sub","Sub-Sub"))+ 
  xlab("Groups") +
  ylab("Microbiome similarity") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")


# aerotolerant
physeqAdistance4groupaerolm2multidata <- ggpredict(physeqAdistance4groupaerolm2multi,terms = "Typeofpair2")
physeqAdistance4groupaerolm2multidata2 <- data.frame(physeqAdistance4groupaerolm2multidata) %>% dplyr::rename(value = predicted, Typeofpair2 = x) 

aerotypeplot2 <- ggplot(physeqAdistance4groupaerolm2multidata2, aes(x=Typeofpair2, y= value)) + geom_point() + geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + stat_halfeye(data=physeqAdistance4groupaero2, aes(x=Typeofpair2,y=value,fill=Typeofpair2),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +geom_text(data=physeqAdistance4ss, aes(x=Typeofpair2, y=-20,label=count), inherit.aes = F)+ scale_x_discrete(labels= c("Dom-Dom","Dom-Help","Dom-Sub","Help-Sub","Sub-Sub"))+ 
  xlab("Groups") +
  ylab("Aerotolerant microbiome similarity") +theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none")


# anaerobes
physeqAdistance4groupanaelm2multidata <- ggpredict(physeqAdistance4groupanaelm2multi,terms = "Typeofpair2")
physeqAdistance4groupanaelm2multidata2 <- data.frame(physeqAdistance4groupanaelm2multidata) %>% dplyr::rename(value = predicted, Typeofpair2 = x) 

anaerotypeplot2 <-ggplot(physeqAdistance4groupanaelm2multidata2, aes(x=Typeofpair2, y= value)) + 
  geom_point() + 
  geom_errorbar(aes(ymin = conf.low, ymax=conf.high),width=0.01) + 
  stat_halfeye(data=physeqAdistance4groupanae2, aes(x=Typeofpair2,y=value,fill=Typeofpair2),inherit.aes = F, adjust = 0.5,justification = -0.2, .width = 0, point_colour=NA, width = 0.5) + 
  scale_fill_manual(values = c("#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE","#D4E2AE")) +
  geom_text(data=physeqAdistance4ss, aes(x=Typeofpair2, y=3,label=count), inherit.aes = F) + 
  geom_signif(comparisons = list(c("Dominantpair","DomSub"),c("Dominantpair","HelpSub"),c("Dominantpair","SubSub")),annotations= c("p = 0.051","p = 0.034","p = 0.014"),y_position = c(-9, -6,-3)) + 
  geom_signif(comparisons = list(c("Dominantpair","DomHelp"),c("DomSub","SubSub")),annotations= c("",""),y_position = c(-41.5, -41.5),tip_length=0) +
  geom_signif(comparisons = list(c("DomHelp","HelpSub")),annotations= c("p = 0.003"),y_position = c(-43),tip_length=-0.04,extend_line = 0.05, vjust = 4) +
  scale_x_discrete(labels= c("Dom-Dom","Dom-Help","Dom-Sub","Help-Sub","Sub-Sub"))+
  xlab("Groups") +
  ylab("Anaerobic microbiome similarity") +
  theme_tufte(base_size = 15, base_family = "Arial") + 
  theme(axis.line = element_line(colour = "black", linetype=1),legend.position="none") + coord_cartesian(ylim=c(-43,3))
anaerotypeplot2

# all plots
shannontype2
pairplot2
ggarrange(pairplot2, typeplot2, widths = c(1.2,2), labels=c("A","B"))
ggarrange(aerotypeplot2,anaerotypeplot2, labels=c("A","B"))

tiff("Output/Figure1.tiff", res =400, units = "in", width = 8, height = 8, bg = "white")
pairplot2
dev.off()

tiff("Output/Figure2.tiff", res =400, units = "in", width = 10, height = 8, bg = "white")
shannontype2
dev.off()

tiff("Output/Figure3.tiff", res =400, units = "in", width = 10, height = 8, bg = "white")
anaerotypeplot2
dev.off()
