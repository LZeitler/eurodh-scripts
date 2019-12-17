source('source_me.R')

## devtools::install_github("petrelharp/local_pca/lostruct")
library(lostruct)


## local windowed PCA
wstart <- 6.5E7
wstop <- 9.5E7
winfilename <- 'data/260419_step6_v4_BU_LR_chr3_mywindow.vcf'
winsize <- 500

system(paste0("vcftools --vcf ~/ma/data/dh/260419_step6_v4_BU_LR.vcf --chr 3 --from-bp ",
              format(wstart,scientific=F)," --to-bp ",
              format(wstop,scientific=F), " --recode --out data/260419_step6_v4_BU_LR_chr3_mywindow"))
system("mv data/260419_step6_v4_BU_LR_chr3_mywindow.recode.vcf data/260419_step6_v4_BU_LR_chr3_mywindow.vcf")
system("ex -sc '1i|##fileformat=VCFv4.0' -cx data/260419_step6_v4_BU_LR_chr3_mywindow.vcf")

snps <- read_vcf(winfilename)
pos <- as.vector(fread(winfilename,data.table=F)[,2])
pos <- pos[seq(1,length(pos)-winsize,winsize)]
pcs <- eigen_windows(snps,k=2,win=winsize)

pdt <- data.frame(pcs) %>%
    mutate(windows=1:nrow(pcs)) %>% 
    select(-total,-lam_1,-lam_2) %>%
    melt('windows')
pdt$Axis <- colsplit(pdt$variable,'_',c('PC','Axis','ind'))[,2]
pdt$ind <- colsplit(pdt$variable,'_',c('PC','Axis','ind'))[,3]
pdt$variable <- NULL
pdt <- inner_join(filter(pdt,Axis==1),filter(pdt,Axis==2),by=c('ind','windows'))
pdt <- inner_join(pdt,data.frame(pos,windows=1:nrow(pcs)))
pdt <- full_join(pdt,data.frame(ind=filter(pdt,windows==4)$ind, # add cluster variable based on #4
                                cluster=cut(filter(pdt,windows==4)$value.x,3)),by='ind') 

a <- ggplot(pdt)+
    geom_point(aes(value.x,value.y,color=cluster))+
    labs(x='Coord1',y='Coord2')+
    theme(legend.position='none')+
    facet_wrap(~pos)
a

saveplot(a,'bu_chr3_localpca-60100-pc12-50120')


## genomewide PCA
gdsfile <- 'data/260419_step6_v4_BU_LR.gds'
## snpgdsVCF2GDS('~/ma/data/dh/260419_step6_v4_BU_LR.vcf',gdsfile)
co <- snpgdsOpen(gdsfile)
pca.co <- snpgdsPCA(co)
snpgdsClose(co)
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
    p <- p + scale_alpha_discrete(range = c(.4, .9))
    p <- p + scale_color_viridis(discrete = T)
    return(p)
}
p12 <- pcaplot(pca.co.tab,1,2,x=ev1,y=ev2,color=Accession,shape=Type)

saveplot(p12,'pca-snprelate-bu-lr-genomewide')


(both <- plot_grid(a,p12+theme(legend.position = 'none'),ncol=1,rel_heights = c(.7,.3),labels = 'AUTO'))

saveplot(both,'pca-inversion-and-genomewide',8,13)
