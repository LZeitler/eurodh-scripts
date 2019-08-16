source('source_me.R')

#### Supplementary Figure S8: shared joint probabilities boxplot
#### prob-shared.pdf

## get probabilities
prodf <- get.prob('output/pro_260419.txt') %>%
    filter(outl==T)

## prepare data 
basepop.top <- left_join(prodf,data.frame(table(prodf$snp)),by=c('snp'='Var1'))

annotation <- data.frame(table(basepop.top$Freq))

## plot
c <- ggplot(basepop.top)
c <- c+geom_boxplot(aes(as.factor(Freq),logprob),fill='grey90')
c <- c+labs(x='Significant outlier SNP shared in populations',
            y='-log10(p)')
c <- c+scale_x_discrete(labels=c('1'=paste(annotation[1,],collapse=', n = '),
                                   '2'=paste(annotation[2,],collapse=', n = '),
                                   '3'=paste(annotation[3,],collapse=', n = '),
                                   '4'=paste(annotation[4,],collapse=', n = '),
                                   '5'=paste(annotation[5,],collapse=', n = ')
                                   ))

saveplot(c,'prob-shared')
