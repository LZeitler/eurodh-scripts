source('source_me.R')

## get probability output from probability_farm_260419.R
freq <- get.freq('output/freq_260419.txt')                    # aSFS
prodf <- get.prob('output/pro_260419.txt') %>% filter(chr==3) # jProb

## prepare aSFS data overview
inout <- data.frame()
for (r in unique(freq$race)){
    inout <- rbind(inout,
                   data.frame(race=r,
                              inCI=nrow(filter(freq,race==r,
                                               in_CI == T)),
                              outCI=nrow(filter(freq,race==r,
                                                in_CI == F)),
                              outCItop=nrow(filter(freq,race==r,
                                                   out.top==T)),
                              outCIbot=nrow(filter(freq,race==r,
                                                   out.bot==T)),
                              max.y=unique (filter(freq,race==r)$sD)))
}
inout <- inout %>% mutate(outperc=round (outCI/(outCI+inCI)*100, 2),
                          out.top.perc=round(outCItop/(outCItop+inCI)*100,2),
                          out.bot.perc=round(outCIbot/(outCIbot+inCI)*100,2))
freq$race <- as.factor(freq$race)

## plot Fig. 2a: aSFS
c <- ggplot(freq)
c <- c + geom_bin2d(aes(x=anc_freq,
                        y=d),
                    size=.3,
                    alpha=.5)
c <- c + geom_ribbon(aes(ymin=bot, ymax=top,
                         x=anc_freq),
                     color = "steelblue",
                     ## fill="steelblue2",
                     fill=NA,
                     alpha=.5)
c <- c + facet_wrap(~race, scales="free", nrow = 1,
                    labeller = as_labeller(racelabels))
c <- c + scale_fill_viridis(name = "Count",
                             trans = "log",
                             guide="legend",
                             ## low = "Yellow",
                             ## high = "blue4",
                             breaks = c(1, 10, 100, 1000, 10000),
                             labels = c(1, 10, 100, 1000, 10000))
c <- c + labs(x = "Ancestral frequency",
              y = "Count (DH)")
c <- c + geom_text(data=inout,
                   aes(x = .15,
                       y = max.y*0.75,
                       label = paste(out.top.perc, "%")))
c <- c + geom_text(data=inout,
                   aes(x = .85,
                       y = max.y*0.25,
                       label = paste(out.bot.perc, "%")))
c <- c + geom_text(data=inout,
                   aes(x = .5,
                       y = max.y/2,
                       label = paste(100-outperc, "%")))
c <- c + theme(legend.key = element_blank(),
               strip.background = element_blank(),
               ## legend.text = element_text(size = 10),
               ## legend.title = element_text(size = 12)
               axis.text = element_text(size = 16),
               axis.title = element_text(size = 20),
               legend.text = element_text(size = 14),
               legend.title = element_text(size = 16))+
    scale_x_continuous(breaks = c(0,.5,1))
f2a <- c

#### plot fig 2b: jProb BU chrom 3

## prepare data 
pop.base <- 'BU'                        # which population should be compared to
basepop <- filter(prodf, race==pop.base)   # baseline population
prodf.top <- prodf %>% group_by(race) %>% filter(logprob>=cutoff) %>% data.frame()
basepop.top <- left_join(filter(basepop,logprob>=cutoff),data.frame(table(prodf.top$snp)),by=c('snp'='Var1'))
basepop.bot <- filter(basepop, logprob<cutoff)

## plot
m <- ggplot()
m <- m + geom_point(data=basepop.bot,
                    aes(x=pos/1000000,
                        y=logprob),
                    color='black',
                    size=.5,
                    alpha=.5,
                    show.legend=F)
m <- m + geom_point(data=basepop.top,
                    aes(x=pos/1000000,
                        y=logprob,
                        color=Freq),
                    size=2.5,
                    alpha=.75,
                    show.legend=T)
m <- m + scale_color_gradientn(colours = rev(heat_hcl(5)))
m <- m + guides(colour = guide_legend(override.aes = list(size=5,alpha=1)))
m <- m + theme(legend.key = element_blank(),
               strip.background = element_blank(),
               axis.text = element_text(size = 16),
               axis.title = element_text(size = 20),
               legend.text = element_text(size = 14),
               legend.title = element_text(size = 16))
m <- m + geom_vline(data=filter(centro, chr==3),
                    aes(xintercept=pos/1000000),
                    linetype="dashed",
                    size=.3)
f2b <- m + labs(y=paste("-log10(p) of", pop.base),x='Position on chromosome 3 [Mbp]', color='Overlap')


## save plot
f2 <- plot_grid(f2a,f2b,
                align = 'vh',
                labels = c("A", "B"),
                hjust = -1,
                nrow = 2)

saveplot(f2, 'fig2-n', height = 7, width = 15)
