source('source_me.R')

#### prepare data
## get haplotype data computed in hapex
dat <- do.call(rbind,lapply(races, function(race){
    data.frame(fread(paste0("output/hapex/hapseq_600k_v4_",race,".vcf.txt"),data.table=F),race)}))
dat$type <- sapply(as.character(dat$race), function(x) unlist(strsplit(x, "_"))[2])
dat$pop <- sapply(as.character(dat$race), function(x) unlist(strsplit(x, "_"))[1])
dat <- dat %>% mutate(hapname=paste0(chrompos, ':', seq))

## identify major haplotype in LR and look it up in DH
hapdi <- data.frame()
for (ra in racess){
    lr <- filter(dat, pop==ra, type=='LR', haps>1) # only windows with mutliple haplotypes
    lr <- lr %>% group_by(chrompos) %>% slice(which.max(freq)) %>% data.frame() # get major haplotype
    dh <- filter(dat, pop==ra, type=='DH', chrompos %in% lr$chrompos) # remove windows with only 1 hap in LR
    dh.inherit <- filter(dh, hapname %in% lr$hapname) # major haplotypes from LR that are present in DH
    dh.fixed   <- filter(dh.inherit, freq==1)
    dh.segreg  <- filter(dh.inherit, freq<1)
    dh.lost    <- filter(lr, !hapname %in% dh.inherit$hapname) %>%
        mutate(type="DH", haps=NA, count=0, freq=0, hapdiv=NA) 
    dh.all     <- rbind(mutate(dh.fixed,  majorhap='fixed'),
                        mutate(dh.segreg, majorhap='segregating'),
                        mutate(dh.lost,   majorhap='lost'))
    combined <- inner_join(lr, dh.all[,c('chrompos','freq','hapdiv','count','majorhap','haps')],
                           by='chrompos', suffix = c('.LR','.DH'))
    combined$freqdiff <- combined$freq.DH - combined$freq.LR # a LR/DH combined df
    hapdi <- rbind(hapdi, combined)                          # output
}

## add small number correction for diversity
hapdi <- hapdi %>% mutate(N.LR=ifelse(race=='BU_LR',22,
                     ifelse(race=='GB_LR',46,23)),
                N.DH=ifelse(race=='BU_LR',36,
                     ifelse(race=='GB_LR',59,
                     ifelse(race=='RT_LR',44,
                     ifelse(race=='SC_LR',58,
                     ifelse(race=='SF_LR',69,NA))))),
                hapdiv.DH=ifelse(is.na(hapdiv.DH),NA,hapdiv.DH*(N.DH/(N.DH-1))),
                hapdiv.LR=hapdiv.LR*(N.LR/(N.LR-1)))

## count haplotype classes, numbers
hapdi <- filter(hapdi, snps>5) # filter out short haplotypes
tab <- data.frame(addmargins(table(select(hapdi, pop, majorhap))))
tab2 <- spread(tab,majorhap,Freq)
tab2$fixed.p <- round(tab2$fixed/tab2$Sum,4)*100 
tab2$lost.p <- round(tab2$lost/tab2$Sum,4)*100
tab2$seg.p <- round(tab2$segregating/tab2$Sum,4)*100




#### Fig 3 A: Haplotype diversity 600k

f3a <- unique(select(dat, chrompos, chr, pos, hapdiv, type, pop, haps)) %>%
    filter(haps>1) %>%
    ggplot()+
    geom_boxplot(aes(x=pop,y=hapdiv,fill=type),alpha=.75)+      
    scale_fill_viridis(discrete = T, begin = .2, end = .8)+
    labs(x="Population",y="Haplotype diversity",fill="Type")
                                                           



#### Fig 3 B: Haplotype fate in BU, Chrom 3

f3b <- ggplot(filter(hapdi, chr==3, pop=='BU'))
f3b <- f3b + geom_jitter(aes(pos/1000000, freq.LR,
                             color=factor(majorhap,
                                 labels = c(paste0('fixed:\n', filter(tab2, pop==ra)$fixed.p,' %'),
                                            paste0('lost:\n', filter(tab2, pop==ra)$lost.p,' %'),
                                            paste0('segregating:\n', filter(tab2, pop==ra)$seg.p,' %')))),
                         alpha=.8, size=.85)
f3b <- f3b + geom_vline(data=filter(centro, chr==3),
                    aes(xintercept=pos/1000000),
                    linetype="dashed",
                    size=.3)
f3b <- f3b + labs(x="Position on chromosome 3 [Mbp]",
                  y="Haplotype frequency (LR)",
                  ## color="Fate of most common\nhaplotype among\nall haplotypes\nin genome")
                  color="Fate of most common\nhaplotype among\nanalysed windows\nin genome") 
f3b <- f3b + theme(legend.key.height=unit(2.5,"line"))
f3b <- f3b + scale_color_manual(values = c('#481567FF','red','#3CBB75FF'))
f3b <- f3b + guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))
f3b <- f3b + coord_cartesian(ylim = c(0,1))


#### Stitch everything together
f3 <- plot_grid(f3a,
                f3b,
                align = 'vh',
                labels = c("A", "B"),
                hjust = -1, vjust = c(1.5,.5),
                nrow = 2)

saveplot(f3, "fig3")
