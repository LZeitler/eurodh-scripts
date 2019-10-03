source('source_me.R')

#### Panel A: PCA

## calculate PC
gdsfile <- 'data/260419_step6_red_v4.gds'
co <- snpgdsOpen(gdsfile)
pca.co <- snpgdsPCA(co)
pca.co.tab <- data.frame(sample.id = pca.co$sample.id,
                         Set = as.factor(substr(pca.co$sample.id,1,5)),
                         Accession = as.factor(substr(pca.co$sample.id,4,5)),
                         Type = as.factor(substr(pca.co$sample.id,1,2)),
                         ev1 = pca.co$eigenvect[,1],
                         ev2 = pca.co$eigenvect[,2],
                         ev3 = pca.co$eigenvect[,3],
                         stringsAsFactors = F)

## prepare plot
pcaplot <- function(data,ev.x,ev.y,...){
    p <- ggplot(data)
    p <- p + geom_point(size=2,alpha=.7,aes(...)) + labs(x = paste0('PC ', ev.x, ' (',
                                                 round(pca.co$varprop[ev.x]*100,2),"%)"),
                                       y = paste0('PC ', ev.y, ' (',
                                                 round(pca.co$varprop[ev.y]*100,2),"%)"))
    p <- p + scale_alpha_discrete(range = c(.4, .9))+
        scale_shape_manual(values = c(2,19,21))
    p <- p + scale_color_viridis(discrete = T)
    return(p)
}
p12 <- pcaplot(pca.co.tab,1,2,x=ev1,y=ev2,color=Accession,shape=Type)
p13 <- pcaplot(pca.co.tab,1,3,x=ev1,y=ev3,color=Accession,shape=Type)
p23 <- pcaplot(pca.co.tab,2,3,x=ev2,y=ev3,color=Accession,shape=Type)



#### Panel B: nucleotide diversity

## calculate site pi for 50k dataset and prepare data
pidf <- data.frame()
for (race in races){
    command <- paste0('vcftools --vcf data/260419_step6_red_v4_',
                      race,
                      '.vcf --site-pi --out output/pi/260419_step6_red_v4_',
                      race,
                      '_55_pi')
    system(command)
    pidf <- rbind(pidf,
                  data.frame(fread(paste0('output/pi/260419_step6_red_v4_',race,'_55_pi.sites.pi')),race))
}
pidf$type <- sapply(as.character(pidf$race), function(x) unlist(strsplit(x, "_"))[2])
pidf$pop <- sapply(as.character(pidf$race), function(x) unlist(strsplit(x, "_"))[1])
pidf <- pidf %>% mutate(chrompos=paste0(CHROM,':',POS))

trim <- data.frame()                    # remove monomorphic sites in LR-DH pairs
for (r in unique(pidf$pop)){
    w1 <- filter(pidf, pop==r, type=="LR", PI==0)$chrompos
    w2 <- filter(pidf, pop==r, type=="DH", PI==0)$chrompos
    w <- w1[w1%in%w2]
    trim <- bind_rows(trim, filter(pidf, pop==r, !chrompos %in% w))
}

## prepare panel B plot
fig1B <- function(x){
    lvls <- c('LR','DH')
    f1 <- ggplot(x)
    f1 <- f1 + geom_boxplot(
                   aes(x=pop,
                       y=PI,
                       alpha=factor(type,levels=lvls),
                       fill=pop))
    f1 <- f1 + stat_summary(aes(pop,PI,group=factor(type,levels=lvls)),
                            alpha=.8,position=position_dodge(.75),
                            fun.y='mean',geom='point',shape=18,size=3)
    f1 <- f1 + scale_x_discrete(labels = gsub(' \\(','\n\\(',racelabels))
    f1 <- f1 + labs(x="Population",
                    y=expression("Nucleotide diversity ("*pi*")"),
                    alpha="Type")
    f1 <- f1 + scale_fill_viridis(discrete = T)
    f1 <- f1 + scale_alpha_discrete(range = c(.3, .75))
    f1 <- f1 + guides(alpha = guide_legend(override.aes = list(fill = 'black')),
                      fill = 'none')
    return(f1)
}
f1b <- fig1B(trim)

## Wilcox tests
tests.wil <- function(data,by=racess){
    t(sapply(by,function(x){
        test <-
            wilcox.test(filter(data,pop==x,type=='LR')$PI,filter(data,pop==x,type=='DH')$PI,'two.sided')
        means <-
            c('LR'=mean(filter(data,pop==x,type=='LR')$PI),'DH'=mean(filter(data,pop==x,type=='DH')$PI))
        c('p'=test$p.value,'pop'=x,round(means,5))}))
}
tests.wil(trim)



#### Panel C: jSFS with original data

## calculate allele frequecies
for (race in races){
    command <- paste0('plink --vcf data/melch_', race,
                      '.vcf --freq counts --out output/counts/melch_', race)
    system(command)
}
mafdata <- data.frame()
for (race in racess){
    mafdata <- rbind(mafdata,data.frame(merge(fread(paste0('output/counts/melch_',race,"_DH.frq.counts"), data.table = F),
                                              fread(paste0('output/counts/melch_',race,"_LR.frq.counts"), data.table = F),
                                              by="SNP", suffixes=c(".DH", ".LR")), race))
}

## plot jSFS
c <- ggplot(mafdata,
            aes(C1.LR/(C1.LR+C2.LR),
                C1.DH/(C1.DH+C2.DH),
                color=race))
c <- c + geom_point(size=.35, alpha = .5)
c <- c + labs(x= "Allele frequency in LR",
              y= "Allele frequency in DH")
c <- c + geom_abline(aes(slope=1, intercept = 0))
c <- c + scale_color_viridis(discrete = T)
c <- c + facet_wrap(~race, labeller = as_labeller(racelabels), nrow = 1, ncol = 5)
fig1big.bot <- c + theme(legend.key = element_blank(), strip.background = element_blank(),
               legend.position = 'none') 



#### Stitch everything together
fig1big.top <- plot_grid(plot_grid(p12+
                                   guides(color=F)+
                                   theme(legend.position = c(.06,.17),
                                         legend.background = element_rect(fill='white',
                                                                          size=.3,
                                                                          linetype='solid',
                                                                          color='black'),
                                         legend.margin = margin(r=.8,l=.8,t=.2,b=.2,unit='lines'),
                                         legend.title = element_blank(),
                                         legend.text = element_text(size=9)),
                                   p13+theme(legend.position = "none")),
                         f1b+theme(axis.title.x = element_blank(),
                                   legend.position = c(.06,.17),
                                   legend.background = element_rect(fill='white',
                                                                    size=.3,
                                                                    linetype='solid',
                                                                    color='black'),
                                   legend.margin = margin(r=.8,l=.8,t=.2,b=.2,unit='lines'),
                                   legend.title = element_blank(),
                                   legend.text = element_text(size=9)),
                         labels = 'AUTO',
                         nrow = 1)

fig1big <- plot_grid(fig1big.top,
                     plot_grid(fig1big.bot,labels='C'),
                     nrow = 2)
fig1big <- plot_grid(fig1big,
                     get_legend(p12+guides(shape=F)),
                     rel_widths = c(.93,.07))
fig1big

## Save the plot
saveplot(fig1big,'fig1_big',16,7)
