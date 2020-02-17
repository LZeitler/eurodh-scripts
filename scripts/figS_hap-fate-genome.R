source('source_me.R')

#### Supplementary Figure S11: Fate of most common haplotype genomewide
#### s3a_all.pdf

## prepare data
source('scripts/make_hapdiv_ld_600.R')

## plot
s3a <- ggplot(hapdi)
s3a <- s3a + geom_jitter(aes(pos/1000000, freq.LR,
                             color=factor(majorhap,
                                 labels = c(paste0('fixed'),
                                            paste0('lost'),
                                            paste0('segregating')))),
                         alpha=.8, size=.35)
s3a <- s3a + geom_vline(data=filter(centro),
                    aes(xintercept=pos/1000000),
                    linetype="dashed",
                    size=.3)
s3a <- s3a + facet_grid(pop~chr,scales = 'free_x',space="free_x",switch="x")
s3a <- s3a + labs(x="Chromosome, Position [Mbp]",
                  y="Haplotype frequency (LR)",
                  color="Fate in DH") 
s3a <- s3a + theme(legend.key.height=unit(2.5,"line"),
                   panel.spacing.x=unit(0, "lines"))
s3a <- s3a + scale_color_manual(values = c('#481567FF','red','#3CBB75FF'))
s3a <- s3a + guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))
s3a <- s3a + coord_cartesian(ylim = c(0,1))
s3a <- s3a + scale_x_continuous(breaks = seq(0,500,100))
s3a

saveplot(s3a,'s3a_all-ld',10*1.5,6*1.5)
