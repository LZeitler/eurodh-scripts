source('source_me.R')

#### Supplementary Figure S10: Haplotype diversity boxplot 50k vs 600k
#### supfig6.pdf

#### prepare data
## get haplotype data computed in hapex
dat <- data.frame()
for (race in races){
    fr1 <- unique(fread(paste0("output/hapex/hapseq_600k_v4_",race,".vcf.txt"),data.table=F,drop=c(1:2,8)))
    fr2 <- unique(fread(paste0("output/hapex/hapseq_55k_v4_",race,".txt"),data.table=F,drop=c(1:2,8)))
    dat <- rbind(dat, data.frame(fr1,race,set='600k'))
    dat <- rbind(dat, data.frame(fr2,race,set='50k'))
}
dat$type <- sapply(as.character(dat$race), function(x) unlist (strsplit(x, "_"))[2])
dat$pop <- sapply(as.character(dat$race), function(x) unlist (strsplit(x, "_"))[1])

s6 <- ggplot(filter(dat, haps>1))
s6 <- s6 + geom_boxplot(
               aes(x=pop,
                   y=hapdiv,
                   fill=paste(set,type)),
               alpha=.75)
s6 <- s6 + scale_fill_manual(values=c("coral4", "coral1","steelblue4", "steelblue1"))
s6 <- s6 + labs(x="Population",
                y="Haplotype diversity",
                fill="Set")

saveplot(s6,'supfig6')

