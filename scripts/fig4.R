source('source_me.R')

freq <- get.freq('output/freq_260419.txt')                    # aSFS
prodf <- get.prob('output/pro_260419.txt')                    # jProb

## SNP LD pruning
system('plink --vcf 260419_step6_red_LR.vcf --indep-pairphase 500 500 0.2')
prune <- as.vector(fread('plink.prune.in',data.table=F,header=F)[,1])

## downsample for comparison landrace LR
stratify <- function(data){
    fr <- data
    out <- data.frame()
    for (ra in racess){
        for (bin in sort(unique(fr$lr_freq.bin))){
            t <- filter(fr,race==ra,outl==T,lr_freq.bin==bin)
            t1 <- filter(fr,race==ra,outl==F,lr_freq.bin==bin)
            n <- min(nrow(t),nrow(t1))
            out <- bind_rows(out,sample_n(filter(fr,race==ra,outl==T,lr_freq.bin==bin),size=n),
                             sample_n(filter(fr,race==ra,outl==F,lr_freq.bin==bin),size=n))
        }
    }
    return(out)
}

## prepare plots
fig.het <- function(data){
    c <- ggplot(data)
    c <- c+geom_violin(aes(race,hetfreq,fill=outl),alpha=.8)
    c <- c+stat_summary(aes(race,hetfreq,group=outl),alpha=.8,
                        position=position_dodge(.9),fun.y='mean',geom='point',shape=18,size=3)
    c <- c+scale_fill_viridis(discrete = T, begin = .2, end = .8, direction = -1,labels=c('non-outlier','outlier'))
    c <- c+labs(x='Accession',y='Heterozygote frequency',fill='SNP')
    c
}
