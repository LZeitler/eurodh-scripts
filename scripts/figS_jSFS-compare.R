source('source_me.R')

require(magick)                         # to combine ggplot with existing png


#### Supplementary Plot S6: jSFS with original data compared to published figure
#### supfig1.pdf

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

## plot using original data
mafplot <- function(data){
    mafdata <- data
    c <- ggplot(mafdata,
                aes(C1.LR/(C1.LR+C2.LR),
                    C1.DH/(C1.DH+C2.DH)))
    c <- c + geom_point(size=.35, alpha = .5)
    c <- c + labs(x= "Allele frequency in LR",
              y= "Allele frequency in DH")
    c <- c + geom_abline(aes(slope=1, intercept = 0))
    c <- c + facet_wrap(~race, labeller = as_labeller(racelabels), nrow = 3, ncol = 2)
    c <- c + theme(legend.key = element_blank(), strip.background = element_blank()) 
    c
}

## plot together
s6 <- plot_grid(mafplot(mafdata),
                ggplot()+draw_image('plots/melch.sfs.png'),
                nrow=1,labels = 'AUTO')

saveplot(s6, 'supfig1', width = 14, height = 10)

