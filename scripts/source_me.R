require(data.table)
require(tidyr)
require(dplyr)
require(cowplot)
require(ggplot2)
require(grid)
require(gridExtra)
require(viridis)
require(colorspace)
require(reshape2)
require(SNPRelate)
options(max.print=999)

setwd('..')

## a function to save plots quickly in pdf and png
saveplot <- function(x, name, width=10, height=6){
    ggsave(paste0("plots/", name, ".pdf"),
           x, width = width, height = height, limitsize = F)
    ggsave(paste0("plots/", name, ".png"),
           x, width = width, height = height, limitsize = F)
}

## a function to get complements
complement <- function(x){ 
    bases=c("A","C","G","T")
    xx <- unlist(strsplit(toupper(x),NULL))
    paste(unlist(lapply(xx,function(bbb){
    if(bbb=="A") compString<-"T"
    if(bbb=="C") compString<-"G"
    if(bbb=="G") compString<-"C"
    if(bbb=="T") compString<-"A"
    if(!bbb %in% bases) compString<-"N"
    return(compString)
})),collapse="")
}

## supply labels and accession lists
races <- c("BU_LR", "GB_LR", "RT_LR", "SC_LR", "SF_LR", "BU_DH", "GB_DH", "RT_DH", "SC_DH", "SF_DH")
racess <- c("BU", "GB", "RT", "SC", "SF")
racelabels <- c("BU" = "Bugard (BU)",
                "GB" = "Gelber Badischer (GB)",
                "RT" = "Rheintaler (RT)",
                "SC" = "Schindelmeiser (SC)",
                "SF" = "Strenzfelder (SF)",
                "WA" = "Wantzenau (WA)")

## centromere positions in V4
centro <- data.frame(chr=1:10,
                     pos=c(mean(c(136.77, 137.12)),
                           mean(c(95.51, 97.49)),
                           mean(c(85.78, 86.93)),
                           mean(c(109.07, 110.5)),
                           mean(c(104.54, 106.82)),
                           mean(c(52.3, 53.11)),
                           mean(c(56.38, 56.68)),
                           mean(c(50.53, 52.07)),
                           mean(c(53.75, 55.39)),
                           mean(c(51.39, 52.78)))*1000000)
