source('source_me.R')

#### generates TABLE S4

## get probability output from probability_farm_260419.R
freq <- get.freq('output/freq_260419.txt')                    # aSFS

## prepare aSFS data overview
inout <- data.frame()
for (r in unique(freq$race)){
    inout <- rbind(inout,data.frame(race=r,
                              inCI=nrow(filter(freq,race==r,in_CI == T)),
                              outCI=nrow(filter(freq,race==r,in_CI == F)),
                              max.y=unique (filter(freq,race==r)$sD)))
}
inout <- inout %>% mutate(outperc=round (outCI/(outCI+inCI)*100, 2))
freq$race <- as.factor(freq$race)

## occurences in populations
occur <- data.frame(table(select(freq,snp,outl))) %>% mutate(outl=as.logical(outl))
fo <- inner_join(filter(freq,outl==T),filter(occur,outl==T))
tab <- as.data.frame.matrix(addmargins(t(table(select(filter(fo),race,Freq)))))

## unique outliers
fo2 <- unique(select(fo,chr,pos,snp,Freq))
tabu <- c(table(fo2$Freq))
tabu <- c(tabu,sum(tabu))

## percentages
perc <- c('Outlier %', paste(inout$outperc, '%'), '','')

## stitch everything together
data.frame(rbind(as.matrix(cbind('Overlap'=c(1:5,'','SUMS'),
                                 rbind(cbind(tab[1:5,]),rep('',5),tab[6,]),
                                 'SUM unique'=c(unname(tabu),'SUMS'))),perc))
