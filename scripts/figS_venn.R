source('source_me.R')

#### Supplementary Figure S9: outlier venn diagram
#### venn_freq_prob.pdf

## load and prepare data
freq <- get.freq('output/freq_260419.txt') %>% select(snp,race,outl.p=outl) 
prodf <- get.prob('output/pro_260419.txt') %>% select(snp,race,outl.f=outl)   
bo <- inner_join(prodf,freq, by = c("snp", "race"))

## draw plot
require(VennDiagram)

vetab <- data.frame(table(select(bo,-snp,-race)))
vetab <- vetab[-1,]

dev.off()
venn.plot <- draw.pairwise.venn(
    area1 = vetab$Freq[1]+vetab$Freq[3],
    area2 = vetab$Freq[2]+vetab$Freq[3],
    cross.area = vetab$Freq[3],
    category = c('jProb\noutlier','aSFS\noutlier'),
    scaled = T,
    fill = c(viridis(5)[4],viridis(5)[2]),
    main.fontfamily="sans",
    fontfamily = 'sans',
    cat.fontfamily = 'sans',
    cex = c(1.5,1.5,1.5),               # font size
    cat.cex = c(1.5,1.5),               # font size
    cat.dist = c(.05,.05),              # label distance
    cat.pos = c(-135,45),               # label position
    ext.dis = c(0,-.05)
)

## save plot pdf
pdf('plots/venn_freq_prob.pdf',width = 7,height = 7)
grid.draw(venn.plot)
dev.off()

## optional: convert to png
system('pdftocairo -png plots/venn_freq_prob.pdf plots/venn_freq_prob.png')
