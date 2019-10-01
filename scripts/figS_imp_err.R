source('source_me.R')

#### Supplementary figures S1 (supfig4.pdf), S2 (supfig5.pdf): Imputational error rates
#### Fig S1: Genomewide imputation error

## data prep
all <- fread("data/260419_step6_red_v4_BU_DH.vcf",
             data.table = F, select = 1:3,col.names = c('chr','pos','id')) %>%
    mutate(Type='Illumina')

## load diagnostic files generated during drop imputation runs
ac <- lapply(racess, function(ra){
    data.frame(fread(paste0('data/imperr/eval.out_full4_',ra,'_plot.diag'), select = c(1,4)), race=ra)
})
ac <- do.call(rbind,ac)
all <- left_join(ac,all,by='id')

## mean errors
aggregate(((1-accur)*100)~race,all,function(x) round(mean(x),1))

## plotting S1 top panels
q <- ggplot(all)
q <- q + stat_summary_bin(aes(pos/1000000, (1-accur)),
                          binwidth = 1.500000,
                          fun.y = mean,
                          geom="line",
                          color = "darkred",
                          size = .5)
q <- q + stat_summary_bin(aes(pos/1000000, (1-accur)),
                          binwidth = 4.500000,
                          fun.y = max,
                          geom="line",
                          color = "DeepSkyBlue3",
                          size = .5)
q <- q + theme(axis.text.x = element_text(size=10, hjust=.06, vjust=.5, angle = -90))
q <- q + labs(x='Chromosome, Position [Mbp]',y="Error rate per window")
q <- q + facet_grid(race~chr, space = 'free_x', scales = 'free_x', switch = 'both')

## plotting S1 bottom panels
q1 <- ggplot(filter(all, Type=='Illumina', race=='BU'))
q1 <- q1 + geom_histogram(binwidth = 1.500000, aes(x = pos/1000000), alpha = .5)
q1 <- q1 + theme(axis.text.x = element_text(size=10, hjust=.06, vjust=.5, angle = -90))
q1 <- q1 + labs(x='Chromosome, Position [Mbp]',y="SNP count")
q1 <- q1 + facet_grid(~chr, space = 'free_x', scales = 'free_x', switch = 'x')

s1 <- plot_grid(q+theme(axis.title.x = element_blank(),
                        axis.text.x = element_blank(),
                        axis.ticks.x = element_blank(),
                        strip.text.x = element_blank()),
                q1,ncol=1,align = 'v', axis = 'lr',
                rel_heights = c(4,1))

## write plot
saveplot(s1,'supfig4', width = 16, height = 9)




#### Fig S2: Imputation error vs recombination

## calculate genetic distances and join with BU error rates
source('scripts/interpolate_map.R')

dt <- filter(all,race=='BU')
dt <- full_join(dt,filter(ac,race=='BU'),by = c("id", "accur", "race")) %>% arrange(chr,pos)

dist <- data.frame()
for (chrom in 1:10){
    t <- filter(dt, chr==chrom)
    tdist <- mapply(function(p1,p2,ch) gendist(p1,p2,ch),
                   t$pos[1:(length(t$pos)-1)],
                   t$pos[2:length(t$pos)],
                   chrom)
    dist <- rbind(dist, data.frame(dist=tdist, chr=chrom, pos=t$pos[1:(length(t$pos)-1)],
                                   accur=t$accur[1:(length(t$pos)-1)]))
}

## calculate Rsq
m <- lm((1-accur) ~ dist, data = dist)
rsq <- format(summary(m)$r.squared, digits = 3)

## plot
s5 <- ggplot(dist,aes(dist,1-accur))
s5 <- s5 + geom_jitter(alpha =.5)
s5 <- s5 + labs(x='Genetic distance from next known SNP [cM]', y='Error rate')
s5 <- s5 + geom_smooth(method = 'lm',se=F)
s5 <- s5 + annotate(geom = 'text', label=paste("R^2==",rsq),
                    x=Inf, y=Inf, vjust=5, hjust=3, parse = T, size=5)

## write plot
saveplot(s5,'supfig5')








