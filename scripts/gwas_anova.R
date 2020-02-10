source('source_me.R')
library(broom)

setwd('output/gwas')
races <- c("GB_LR","RT_LR","SF_LR","GB_DH","RT_DH","SF_DH") # overlapping lines with core dataset
racess <- c("GB","RT","SF")

## get allele frequecies of GWAS dataset
ca <- data.frame()
for (race in racess){
    dhb <- fread(paste0("../counts/221118_step6_red_", race, "_DH.frq.counts"), data.table=F)
    lrb <- fread(paste0("../counts/221118_step6_red_", race, "_LR.frq.counts"), data.table=F)
    ca <- rbind(ca,data.frame(inner_join(dhb,lrb,suffix = c('.DH','.LR'),by='SNP'),race))
}                                       # the tracked allele is A1.DH==A1.LR and is in all races the same


## read GCTB output/ effect sizes
files <- list.files(pattern="*.snpRes") # loads all
ef <- lapply(files,function(x) data.frame(fread(x,data.table = F),trait=substr(x,6,16)))
ef <- do.call(rbind,ef) %>%
    rename(effect=Effect,SNP=Name) %>%
    inner_join(unique(select(ca,SNP,A1.LR)),by=c('SNP')) %>%
    mutate(effect=ifelse(Allele==A1.LR,effect,-1*effect))

# combine frequecies, effect sizes and outlier class
comb <- get.freq('../freq_221118.txt') %>% filter(!race=='WA') %>%
    inner_join(select(ef,SNP,effect,trait,PIP),by=c('snp'='SNP')) %>%
    mutate(lrfreqn=l/(2*sL))
comb$lrfreq <- cut(comb$lrfreqn,10)     # binning

comb <- comb %>%
    mutate(SNP=ifelse(outl,'outlier','non-outlier'),
           abseffect=abs(effect),
           lrfreqbin=as.numeric(lrfreq)/10-.05)  %>%  # mean allele freq in bin
    group_by(race,trait,SNP,lrfreq,lrfreqbin) %>%
    mutate(meanabseffect=mean(abseffect))

## fit model
traits <- unique(comb$trait)
traitlabels <- c('shoot vigor', 'female flowering','fusarium','grain yield','oil content','plant height','protein content')
names(traitlabels) <- traits

res1 <- comb %>%
    group_by(trait) %>%
    do(tidy(aov(abseffect~outl*lrfreq,data=.))) %>%
    data.frame
res1[,3:7] <- apply(res5[,3:7],2,function(x) signif(x,3))
names(res1) <- c('Trait','Term','df','SumSq','MeanSq','F-value','p-value')

## print table
print(tibble(res1),n=28)

## output trait specific
lapply(traits, function(x) filter(res1,Trait==x))
           
fwrite(res1,'gwas_aov_table.txt')

#### explain interaction of frequency and outlier

g <- ggplot(comb)+
    geom_point(aes(lrfreqbin,meanabseffect,
                             color=SNP,
                             shape=SNP))+
    geom_smooth(aes(lrfreqbin,abseffect,color=SNP),method='lm')+
    facet_wrap(~trait,scales='free',nrow=3,labeller=as_labeller(traitlabels))+
    labs(x='allele frequency bin',y='mean absolute effect')
g

## ggsave('~/ma/r/plot/gwas-abseffect-freq-interaction.pdf',g,width = 10,height = 6)
setwd('../..')
saveplot(g,'gwas-freq-vs-effects',12,7)

