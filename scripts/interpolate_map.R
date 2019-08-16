map <- fread("data/map/ogut_fifthcM_map_agpv4_INCLUDE.txt", data.table=F)
map <- select(map, -V1, -V2) %>% rename('chr'='V4', 'pos'='V5', 'cm'='V3') # corrected version V4

interpolate <- list()
chroms <- 1:10
for (ch in chroms){
    sel <- map %>% select(chr, pos, cm) %>% filter(chr == ch)
    helper <- data.frame(chr=c(ch,ch),
                         pos=c(0,max(sel$pos)),
                         cm=c(min(subset(sel$cm,sel$chr==ch), na.rm=T),
                              max(subset(sel$cm,sel$chr==ch), na.rm=T)))
    temp <- rbind(helper,filter(sel, chr==ch))
    temp <- temp[order(temp$pos),]
    interpolate[[ch]] <- approxfun(temp$pos,temp$cm)
}

gendist <- function(pos1, pos2, chrom){
    x <- ifelse(pos1>pos2, pos1, pos2)
    y <- ifelse(pos1>pos2, pos2, pos1)
    interpolate[[chrom]](x) - interpolate[[chrom]](y)
}
