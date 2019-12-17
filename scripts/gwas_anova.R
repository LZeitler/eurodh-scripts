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
    mutate(lrfreq=factor(round(l/(2*sL),1))) # this makes 11 allele frequency bins

## fit model
traits <- unique(comb$trait)

res1 <- comb %>% group_by(trait) %>% 
    do(tidy(aov(effect~outl*lrfreq,data=.))) %>%
    data.frame
res1[,4:7] <- apply(res1[,4:7],2,function(x) signif(x,3))
names(res1) <- c('Trait','Term','df','SumSq','MeanSq','F-value','p-value')

## print table
print(tibble(res1),n=28)

## output trait specific
lapply(traits, function(x) filter(res1,Trait==x))
           
fwrite(res1,'gwas_aov_table.txt')