source('source_me.R')

#### Supplementary figure S5: Boxplots with outlier overlaps vs. allele frequencies (LR, DH, diff)
#### freq-boxplot-combined.pdf

## get aSFS output from probability_farm_260419.R
freq <- get.freq('output/freq_260419.txt')        # aSFS

## prepare plotting data
occur <- data.frame(table(select(freq,snp,outl))) # count overlaps
freq$dh_freq <- round(freq$d/freq$sD,2)           # DH freq
freq$lr_freq <- round(freq$l/(2*freq$sL),2)       # LR freq
ou <- inner_join(filter(freq,outl==T),filter(occur,outl==T,Freq>0),by='snp')
ou$delta <- ou$dh_freq-ou$l/(2*ou$sL)

## prepare plot
freq.box <- function(plt.data,y,ann.data,x.label,y.label){
    c2 <- ggplot(plt.data,aes_string('factor(Freq)',y))
    c2 <- c2+geom_boxplot(fill='grey90')
    c2 <- c2+labs(x=x.label,y=y.label)
    c2 <- c2+stat_summary(fun.y=mean,geom='point',shape=18,size=5)
    c2 <- c2+scale_x_discrete(labels=c('1'=paste(ann.data[1,],collapse=', n = '),
                                       '2'=paste(ann.data[2,],collapse=', n = '),
                                       '3'=paste(ann.data[3,],collapse=', n = '),
                                       '4'=paste(ann.data[4,],collapse=', n = '),
                                       '5'=paste(ann.data[5,],collapse=', n = ')
                                       ))
    c2
}

## data for panel A: DH freq boxplot
ou.dh <- aggregate(dh_freq~snp*Freq,ou,mean)

## data for panel B: LR freq boxplot
ou.lr <- aggregate(lr_freq~snp*Freq,ou,mean)

## data for panel C: freq diff boxplot
ou.de <- aggregate(delta~snp*Freq,ou,mean)

## overlapping outliers
ou.t <- data.frame(table(ou.dh$Freq))

## draw plots
freq.box.dh <- freq.box(ou.dh,'dh_freq',ou.t,'Number of shared outliers','allele frequency of outlier SNPs in DH')
freq.box.lr <- freq.box(ou.lr,'lr_freq',ou.t,'Number of shared outliers','allele frequency of outlier SNPs in LR')
freq.box.de <- freq.box(ou.de,'delta',ou.t,'Number of shared outliers',expression('Frequency difference ('*Delta*'p)'))

## stitch together
cplot <- plot_grid(freq.box.lr+theme(axis.title.x=element_blank(),axis.text.x=element_blank()),
                   freq.box.dh+theme(axis.title.x=element_blank(),axis.text.x=element_blank()),
                   freq.box.de,
                   labels = 'AUTO',
                   ncol = 1)

## save
saveplot(cplot,'freq-boxplot-combined',height = 12)
