source('source_me.R')

#### Supplementary Figure S12: genetic load estimation using GERP scores, recessive and additive models
#### genespace-box-combined.pdf

## read data generated in genespace_euler.R
sum.add <- fread('output/genload/genespace_260419v4_600_additive_seed123.txt',data.table=F)
sum.rec <- fread('output/genload/genespace_260419v4_600_recessive_seed123.txt',data.table=F)

## rename
sum.add[which(sum.add[,'outl']=='no outlier'),'outl'] <- 'non-outlier'
sum.rec[which(sum.rec[,'outl']=='no outlier'),'outl'] <- 'non-outlier'

## prepare test
spacetest <- function(data,ptype='DH',pmodel){
    t <- do.call(rbind,lapply(racess,function(r) {
        tt <- t.test(filter(data,race==r,type==ptype,outl=='outlier',model==pmodel)$sum,
                      filter(data,race==r,type==ptype,outl=='non-outlier',model==pmodel)$sum)
        c(ifelse(tt$p.value<0.05,'*',' '),mean(tt$estimate))}))
    data.frame(race=racess,type=ptype,model=pmodel,code=t[,1],emean=as.numeric(t[,2]))
}

## prepare plot
p.agr.boxplot2 <- function(data){
    agr <- aggregate(sum~race*outl*ind*type*model,data,mean)
    tests <- rbind(spacetest(agr,'LR','Recessive Model'),spacetest(agr,'DH','Recessive Model'),
                   spacetest(agr,'LR','Additive Model'),spacetest(agr,'DH','Additive Model'))
    agr$type <- factor(agr$type,levels=c('LR','DH'))
    p <- ggplot(agr)+
        geom_boxplot(aes(race,sum,fill=outl),outlier.alpha=1,alpha=.8)+
        geom_text(data=tests,aes(race,emean+1,label=code),size=8,color='orange')+
        stat_summary(aes(race,sum,group=paste(type,outl)),
                     position=position_dodge(.75),
                     fun.y='mean',geom='point',shape=18,size=3)+
        scale_alpha_discrete(range = c(.1, .95))+
        guides(alpha = guide_legend(override.aes = list(fill = 'black')))+
        labs(x='Accession',y='Mean genetic load', # y='mean GERP sum in 1 cM region'
             fill='SNP', alpha='Population')+
        theme(strip.background = element_rect(fill='grey95'),
              strip.text = element_text(margin=margin(.3,.3,.3,.3,'lines')),
              strip.text.x = element_text(size=13))
    return(p)
}

## plot everything
p.agr.grid <- plot_grid(
    p.agr.boxplot2(rbind(data.frame(filter(sum.rec,type=='LR'),model='Recessive Model'),
                         data.frame(filter(sum.rec,type=='DH'),model='Recessive Model'),
                         data.frame(filter(sum.add,type=='LR'),model='Additive Model'),
                         data.frame(filter(sum.add,type=='DH'),model='Additive Model')))+
    facet_grid(model~type,scales='free_y',switch='y')
) 
p.agr.grid

## save
saveplot(p.agr.grid,'genespace-box-combined')

