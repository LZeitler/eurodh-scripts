source('source_me.R')

#### Supplementary Figure S12: genetic load estimation using GERP scores, recessive and additive models
#### genespace-box-combined.pdf

## read data generated in genespace_euler.R
sum.add <- fread('output/genload/genespace_260419v4_600_additive1.txt',data.table=F)
sum.rec <- fread('output/genload/genespace_260419v4_600_recessive1.txt',data.table=F)

## prepare plot
p.agr.boxplot2 <- function(data){
    agr <- aggregate(sum~race*outl*ind*type*model,data,mean)
    agr$type <- factor(agr$type,levels=c('LR','DH'))
    p <- ggplot(agr)+
        geom_boxplot(aes(race,sum,fill=outl),outlier.alpha=1,alpha=.8)+
        stat_summary(aes(race,sum,group=paste(type,outl)),
                     position=position_dodge(.75),
                     fun.y='mean',geom='point',shape=18,size=3)+
        scale_alpha_discrete(range = c(.1, .95))+
        guides(alpha = guide_legend(override.aes = list(fill = 'black')))+
        labs(x='Accession',y='Mean genetic load', # y='mean GERP sum in 1 cM region'
             fill='SNP type', alpha='Population')+
        theme(strip.background = element_rect(fill='grey95'),
              strip.text = element_text(margin=margin(.3,.3,.3,.3,'lines')),
              strip.text.x = element_text(size=13))+
        scale_fill_viridis(discrete = T, begin = .2, end = .8, direction = -1)
    p
}

## plot everything
p.agr.grid <- plot_grid(
    p.agr.boxplot2(rbind(data.frame(filter(sum.rec,type=='LR'),model='Recessive Model'),
                         data.frame(filter(sum.rec,type=='DH'),model='Recessive Model'),
                         data.frame(filter(sum.add,type=='LR'),model='Additive Model'),
                         data.frame(filter(sum.add,type=='DH'),model='Additive Model')))+
    facet_grid(model~type,scales='free_y',switch='y')
) 

## save
saveplot(p.agr.grid,'genespace-box-combined')

