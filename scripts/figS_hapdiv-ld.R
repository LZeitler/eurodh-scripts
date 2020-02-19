source('source_me.R')

#### prepare data
source('scripts/make_hapdiv_ld.R')

## count haplotype classes, numbers
tab <- data.frame(addmargins(table(select(hapdi, pop, majorhap))))
tab2 <- spread(tab,majorhap,Freq)
tab2$fixed.p <- round(tab2$fixed/tab2$Sum,4)*100 
tab2$lost.p <- round(tab2$lost/tab2$Sum,4)*100
tab2$seg.p <- round(tab2$segregating/tab2$Sum,4)*100
mns <- hapdi %>% group_by(type,pop) %>% summarise(mean(snps)) %>% data.frame
mean(as.numeric(mns[,3]))


#### Fig 3 A: Haplotype diversity 600k

lvls <- c('LR','DH')
f3a <- unique(select(dat, chrompos, chr, pos, hapdiv, type, pop, haps)) %>%
    filter(haps>1) %>%
    ggplot()+
    geom_boxplot(aes(x=pop,y=hapdiv,alpha=factor(type,levels=lvls),fill=pop))+
    scale_alpha_discrete(range = c(.3, .75))+
    scale_fill_viridis(discrete = T)+
    labs(x="Population",y="Haplotype diversity",alpha="Type")+
    guides(alpha = guide_legend(override.aes = list(fill = 'black')), fill = 'none')
                                                           



#### Fig 3 B: Haplotype fate in BU, Chrom 3

f3b <- ggplot(filter(hapdi, chr==3, pop=='BU'))
f3b <- f3b + geom_jitter(aes(pos/1000000, freq.LR,
                             color=factor(majorhap,
                                 labels = c(paste0('fixed:\n', filter(tab2, pop==ra)$fixed.p,' %'),
                                            paste0('lost:\n', filter(tab2, pop==ra)$lost.p,' %'),
                                            paste0('segregating:\n', filter(tab2,
                                                                            pop==ra)$seg.p,' %')))),
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

saveplot(f3, "fig3-ld-50k")
