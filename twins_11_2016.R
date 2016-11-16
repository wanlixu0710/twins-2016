setwd("~/R_directory/MAPS/6_14_16_data/")
library(plyr)
library(dplyr)
library(ggplot2)
library(reshape)
library(vegan)
library(indicspecies)
rm(list=ls())

# read in the data
### combinning the data
map <- read.csv(file="map_6_14_16.csv",sep=",",header=T)
nseq <- read.csv(file="otu_table_stats.csv",sep=",", header=T)
map <- inner_join(map,nseq)
stool <- read.csv(file="Updated_5_23_16 MAPS_stool_MCB_processing set up.csv", sep=",", header=T)
stool <- stool[,(1:5)]
map <- inner_join(map, stool)
map <- map %>% group_by(Daily_Code) %>% top_n(1, nsequence) %>% filter(nsequence>=10000) #pick the highest number of sequencing for each daily code
map <- select(map, Sample, Subject_ID, Fecal_Sample_ID, PNA)

## read in the a-diversity
a.div <- read.csv(file="alpha_even_10000.csv",sep=",", header=T)
a.div <- inner_join(map, a.div)
a.div$PNA <- as.character(a.div$PNA)
a.div$Subject_ID <- as.character(a.div$Subject_ID)
a.div.30 <- subset(a.div, as.numeric(PNA)<=30)

## demo, nnns
demo <- read.csv(file="MAPS demographic data_3_11_16.csv")
demo <- demo[,c("Subject_ID", "BIrth_GA", "female", "Baby_Race", "vaginal", "PROM", "Twins", "Birth_weight", "SNAPEII")]
demo$Subject_ID <- as.character(demo$Subject_ID)
nnns <- read.csv(file="NNNS_final_8_25_15.csv")
nnns$Subject_ID <- as.character(nnns$Subject_ID)
demonnns <- inner_join(demo, nnns)
a.div.30 <- left_join(a.div.30, demonnns)

## feeding
count <- read.csv(file="count.for.upload.csv", sep=",", header=T, stringsAsFactors = F)
count$Subject_ID <- as.character(count$Subject_ID)
count <- select(count, Subject_ID, day, ten.days.fed)
a.div.30$day[a.div.30$PNA < 11] <- "1st 10 days"
a.div.30$day[a.div.30$PNA> 10 & a.div.30$PNA <= 20] <- "2nd 10 days"
a.div.30$day[a.div.30$PNA> 20 & a.div.30$PNA <= 30] <- "3rd 10 days"
a.div.30 <- inner_join(a.div.30, count)

a.div.30$subject.recode <- revalue(as.character(a.div.30$Subject_ID),  c("1"="A", "2"="B","3"="C","4"="D","5"="E","6"="F", "7"="G","8"="H","9"="I","10"="J","11"="K","12"="L", "13"="M","14"="N","15"="O","16"="P","18"="Q","19"="R","20"="S","21"="T","22"="U","23"="V","24"="W","25"="X", "26"="Y","27"="Z","28"="ZA","29"="ZB","30"="ZC", "35"="ZD","36"="ZE","40"="ZF", "41"="ZG"))

a.div.30$mo <- revalue(as.character(a.div.30$ten.days.fed), c("1"="1.MOM", "2"="4. HDM", "3"="5.Formula", "4"='2.MOM+HDM', "5"="3.MOM+Formula", "6"="6.HDM+Formula"))

## antibiotics use, numbers of days of antibiotics uses during the first 10 days
anti <- read.csv(file="Daily_antibiotic.40.50.csv",stringsAsFactors = F)
anti$Subject_ID <- as.character(anti$Subject_ID)
anti$PNA <- as.character(anti$PNA)
anti.10 <- anti%>% filter(as.numeric(PNA)<11) %>%count(Subject_ID)%>% plyr::rename(c("n"= "days.of.antibiotics.f10"))
a.div.30 <- left_join(a.div.30, anti.10)
a.div.30$days.of.antibiotics.f10[is.na(a.div.30$days.of.antibiotics.f10)] <- 0
a.div.30 <- left_join(a.div.30, anti[,c("Subject_ID", "PNA", "antibiotic.use")])
a.div.30$antibiotic.use[is.na(a.div.30$antibiotic.use)] <- 0

## niss
niss <- read.csv(file="niss_first_46.csv")
niss$Subject_ID <- as.character(niss$Subject_ID)
niss$PNA <- as.character(niss$PNA)
# calculate niss acute, chronic and contact scores.
niss.daily  <- niss%>% mutate(freq.acute.daily=Acute5+Acute4+Acute3+Acute2, freq.chronic.daily=Chronic4+Chronic3+Chronic2)%>%select(Subject_ID,PNA, Acute_Daily, Chronic_Daily, freq.acute.daily, freq.chronic.daily, Kangaroo_Care_Daily, Breastfeeding_Daily, Contact_Total_Daily, Weight) %>%rename(c("Acute_Daily"="weighted.acute.daily", "Chronic_Daily"= "weighted.chronic.daily"))
a.div.30 <- left_join(a.div.30,niss.daily)
rm(nseq, nnns, niss, count, anti, anti.10, demo, demonnns, a.div)


#### twins
a.div.twins <- a.div.30 %>% filter(Twins!="0")
niss.daily$day[niss.daily$PNA < 11] <- "1st 10 days"
niss.daily$day[niss.daily$PNA> 10 & niss.daily$PNA <= 20] <- "2nd 10 days"
niss.daily$day[niss.daily$PNA> 20 & niss.daily$PNA <= 30] <- "3rd 10 days"
## calculate moving average for the weighted NISS value by 3
library(zoo)
niss.daily <- niss.daily%>% arrange(Subject_ID,PNA) %>% group_by(Subject_ID)%>% mutate(mova.acute=rollmean(weighted.acute.daily,3,na.pad=TRUE, align="right"))%>% mutate(mova.chronic=rollmean(weighted.chronic.daily,3,na.pad=TRUE, align="right"))

niss.daily <- inner_join(niss.daily, unique(a.div.twins[,c("Subject_ID", "Twins", "day")]))
niss.day <- niss.daily %>% dplyr::group_by(Twins,day) %>% dplyr::summarize(acute.niss.mean = mean(weighted.acute.daily), chronic.niss.mean=mean(weighted.chronic.daily))
niss.day <- inner_join (niss.day, unique(niss.daily[,c("Subject_ID", "Twins", "day")]))




# plot the acute niss
ggplot (niss.daily, aes(x=as.factor(Twins),y=weighted.acute.daily))+geom_boxplot()
ggplot (niss.daily, aes(x=as.factor(Twins),y=weighted.chronic.daily))+geom_boxplot()

# plot the simpson diversity
ggplot (a.div.twins, aes(x=as.factor(Twins),y=PNA,colour=simpson))+geom_point(size=5)+scale_colour_gradient(low="green", high="blue")+labs(x="Infant", y="Postnatal age (day)", colour="Simpson", title="Gini-Simpson Index for Twins")+ theme_bw() 
library(tidyr)

niss.day <- separate(niss.day, Twins, c("twin", "pair"))
niss.day1 <- niss.day %>% group_by(twin,day) %>% top_n(1, acute.niss.mean) %>% mutate(acute.niss="A")
niss.day <- full_join(niss.day, niss.day1) 
niss.day$acute.niss[is.na(niss.day$acute.niss)] <-" B"
niss.day$acute.niss[niss.day$Subject_ID==31] <- "C"

niss.day1 <- niss.day %>% group_by(twin,day) %>% top_n(1, chronic.niss.mean) %>% mutate(chronic.niss="A")
niss.day <- full_join(niss.day, niss.day1) 
niss.day$chronic.niss[is.na(niss.day$chronic.niss)] <-" B"
niss.day$chronic.niss[niss.day$Subject_ID==31] <- "C"

### a-diversity (average for every 10 days, check the difference of simpson between group high and group low group)
a.div.twins <- left_join(a.div.twins, niss.day) 
niss.ave <- a.div.twins %>% select(twin, pair,day, simpson) %>% group_by (twin,pair, day) %>% summarize (ave.simp=mean(simpson))
niss.ave <- inner_join(niss.day,niss.ave)
summary(aov(data=niss.ave[!niss.ave$acute.niss=="C",],ave.simp~acute.niss+day)) # didn't show any difference
summary(aov(data=niss.ave[!niss.ave$chronic.niss=="C",],ave.simp~chronic.niss+day)) # didn't show any difference

a <- niss.ave%>% select(twin, day, acute.niss, ave.simp) %>% tidyr::spread( acute.niss, ave.simp)
t.test(a$B, a$A, paired = T, na.action=T) # acute doesn't yeild due the small sample sise
a <- niss.ave%>% select(twin, day, chronic.niss, ave.simp) %>% tidyr::spread( chronic.niss, ave.simp) # not working, since equal values of chronic pain between the twins
t.test(a$B, a$A, paired = T, na.action=T)



### b-diversity
## look at the b-diversity of all the twins
otu <- read.csv(file="even_table_10000_L6.csv",  header=T, sep=",")
otu <- inner_join(map,otu)
otu$Subject_ID <- as.character(otu$Subject_ID)
otu <- inner_join(niss.daily, otu)
otu <- inner_join(niss.day,otu)
row.names(otu) <- otu$Fecal_Sample_ID

nms <- metaMDS(otu[,(24:85)], dist="bray", k=2, trymax=250, wascores=TRUE, trymin=50)
stressplot(nms)
stat <- data.frame(nms$points)
stat$Fecal_Sample_ID <- rownames(otu)
stat <- inner_join(stat, otu)
ggplot(stat, aes(x=MDS1, y=MDS2, colour=factor(twin), shape=factor(pair)))+geom_point()+facet_wrap(twin~day)

# use weighted niss daily
adonis(otu[,(24:85)]~otu$twin+otu$PNA+otu$weighted.acute.daily+otu$weighted.chronic.daily, method="bray", strata = otu$Subject_ID)
adonis(otu[,(24:85)]~otu$twin+otu$PNA+otu$weighted.acute.daily+otu$weighted.chronic.daily, method="jaccard", strata = otu$Subject_ID)

# use average weighted niss by each 10 days ** not a factor
adonis(otu[,(24:85)]~otu$twin+otu$PNA+otu$acute.niss.mean+otu$chronic.niss.mean, method="bray", strata = otu$Subject_ID)
adonis(otu[,(24:85)]~otu$twin+otu$PNA+otu$acute.niss.mean+otu$chronic.niss.mean, method="jaccard", strata = otu$Subject_ID)

# use moving average of the prior 3 days of weighted niss ** not a factor
otu$mova.acute[is.na(otu$mova.acute)] <- otu$weighted.acute.daily # replace the NA with the value of the that day
otu$mova.chronic[is.na(otu$mova.chronic)] <- otu$weighted.chronic.daily # replace the NA with the value of the that day
adonis(otu[,(24:85)]~otu$twin+otu$PNA+otu$mova.acute+otu$mova.chronic, method="bray", strata = otu$Subject_ID)
adonis(otu[,(24:85)]~otu$twin+otu$PNA+otu$mova.acute+otu$mova.chronic, method="jaccard", strata = otu$Subject_ID)

## calculate moving average for the weighted NISS value by 2
niss.daily <- niss.daily%>% arrange(Subject_ID,PNA) %>% group_by(Subject_ID)%>% mutate(mova.acute=rollmean(weighted.acute.daily,2,na.pad=TRUE, align="right"))%>% mutate(mova.chronic=rollmean(weighted.chronic.daily,2,na.pad=TRUE, align="right"))
otu2 <- otu%>%select(Subject_ID, PNA, twin, pair, 24:85)
otu2 <- inner_join(niss.daily, otu2)

otu2$mova.acute[is.na(otu2$mova.acute)] <- otu2$weighted.acute.daily # replace the NA with the value of the that day
otu2$mova.chronic[is.na(otu2$mova.chronic)] <- otu2$weighted.chronic.daily # replace the NA with the value of the that day
adonis(otu2[,(17:78)]~otu2$twin+otu2$PNA+otu2$mova.acute+otu2$mova.chronic, method="bray", strata = otu2$Subject_ID)
adonis(otu2[,(17:78)]~otu2$twin+otu2$PNA+otu2$mova.acute+otu2$mova.chronic, method="jaccard", strata = otu2$Subject_ID)


