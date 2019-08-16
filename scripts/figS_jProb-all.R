source('source_me.R')

#### Supplementary Figure S7: joint probabilities genome wide
#### supfig2.pdf

## get probabilities
prodf <- get.prob('output/pro_260419.txt') %>% mutate(para = ifelse(chr %in% c(1,3,5,7,9), T,F)) # jProb


## plot
m <- ggplot(prodf)
m <- m + geom_point(aes(x=pos,
                        y=-log10(prob),
                        color=para),
                    size=.3,
                    alpha=.5,
                    show.legend=F)
m <- m + facet_grid(race~chr, space="free_x", scales="free",switch="x")
m <- m + theme(panel.spacing.x=unit(0, "lines"),
               axis.ticks.x = element_blank(),
               axis.text.x = element_blank())
m <- m + scale_color_manual(values = c("grey30","black"))
m <- m + labs(y="-log10(p)",x='Chromosome')

## save plot
saveplot(m, 'supfig2')
