source('source_me.R')
library(multcompView);library(emmeans)  # for posthoc test and letters

#### Fig 4 A (recessive load on haplotypes), Fig 4 B (Heterozygosity for aSFS),
#### Supplementary Fig (additive load on haplotypes), Supplementary Fig (Heterozygosity for jProb)

#### load data

freq <- get.freq('output/freq_260419.txt')                    # aSFS
prodf <- get.prob('output/pro_260419.txt')                    # jProb

hapsums <- fread('output/genload/genespace_260419v4_600_hapsums.txt',data.table=F)

#### prepare data

## calculate heterozygosity
hetcount <- function(prefix){
    plink <- 'plink'                        # plink command for local
    for (race in races){
        command <- paste0(plink, ' --vcf ', prefix,
                          race,'.vcf --freqx --out output/het/260419_step6_red_freqx_', race)
        system(command)
    }
    ## read in data
    counts <- data.frame()
    for (race in racess){
        dhb <- fread(paste0("output/het/260419_step6_red_freqx_", race, "_DH.frqx"), data.table=F)
        lrb <- fread(paste0("output/het/260419_step6_red_freqx_", race, "_LR.frqx"), data.table=F)
        temp <- rbind(data.frame(dhb, Type = "DH", race), data.frame(lrb, Type = "LR", race))
        counts <- rbind(counts, temp)
    }
    counts$total <- rowSums(counts[,5:10])
    counts$hetfreq <- counts$C.HET./counts$total
    
    counts <- filter(counts,Type=='LR')
    return(counts)
}

counts <- hetcount('data/260419_step6_red_v4_') # positions v4
## fwrite(counts,'output/het/260419_step6_red_v4_freqx.txt',sep='\t')

counts <- fread('output/het/260419_step6_red_v4_freqx.txt',data.table=F)

## prepare SNP LD pruning
system('plink --vcf data/260419_step6_red_v4_LR.vcf --indep-pairphase 500 500 0.2 --out output/prune/prune')
prune <- as.vector(fread('output/prune/prune.prune.in',data.table=F,header=F)[,1])

## prep aSFS data
freq.het <- inner_join(freq,counts,by=c('race'='race','snp'='SNP'))
freq.het <- rbind(filter(freq.het,outl==F,snp%in%prune),filter(freq.het,outl==T)) # apply pruning
freq.het$lr_freq.bin <- round(freq.het$l/(2*freq.het$sL),1)

## prep probability data
prob.het <- inner_join(prodf,counts,by=c('race'='race','snp'='SNP'))
prob.het <- rbind(filter(prob.het,outl==F,snp%in%prune),filter(prob.het,outl==T)) # apply pruning
prob.het <- inner_join(prob.het,                                                  # get allele frequencies
                       select(freq.het,race,snp,lr_freq.bin),by=c('race'='race','snp'='snp'))

## downsample function
stratify <- function(data){
    out <- data.frame()
    for (ra in racess){
        for (bin in sort(unique(data$lr_freq.bin))){
            t <- filter(data,race==ra,outl==T,lr_freq.bin==bin)
            t1 <- filter(data,race==ra,outl==F,lr_freq.bin==bin)
            n <- min(nrow(t),nrow(t1))
            out <- bind_rows(out,
                             sample_n(filter(data,race==ra,outl==T,lr_freq.bin==bin),size=n),
                             sample_n(filter(data,race==ra,outl==F,lr_freq.bin==bin),size=n))
        }
    }
    return(out)
}

#### plotting

## prepare test for heterozygosity
het.test <- function(sample){
    out <- data.frame()
    for (ra in racess){
        test <-
            ifelse(mean(filter(sample,race==ra,outl==F)$hetfreq)<mean(filter(sample,race==ra,outl==T)$hetfreq),T,F)
        out <- rbind(out,
                     data.frame(t=test,race=ra))
    }
    return(out)
}
shuffl.test <- function(data,n=1000){
    te <- do.call(rbind,lapply(1:n,function(x){
        cat(x/n*100,'% \r'); flush.console()
        return(het.test(stratify(data)))
    }))
    te$sig <- te$t
    tests <- aggregate(sig~race,filter(te,sig==T),function(x) 1-length(x)/(nrow(te)/5)) %>%
        mutate(code=ifelse(sig<0.05,'*',' '))
    tests
}

## run tests
tests.freq <- shuffl.test(freq.het)
tests.prob <- shuffl.test(prob.het)

save(list=c("tests.freq","tests.prob"),file="output/heterozygosity_tests.Rdata")
load("output/heterozygosity_tests.Rdata")

## prepare plots
fig.het <- function(data){
    whichtest <- if (names(data)[1]=='prob') 'jProb test' else 'aSFS test'
    tests <- if (names(data)[1]=='prob') tests.prob else tests.freq
    c <- ggplot(data)
    c <- c+geom_violin(aes(race,hetfreq,
                           ## alpha=outl,
                           color=race,
                           fill=outl),
                       ## fill='black',
                       lwd=1)
    c <- c+geom_text(data=tests,aes(race,.8,label=code),size=7)
    c <- c+stat_summary(aes(race,hetfreq,group=outl),alpha=.8,
                        position=position_dodge(.9),fun.y='mean',geom='point',shape=18,size=3)
    ## c <- c+scale_fill_viridis(discrete = T, begin = .2, end = .8, direction =
    ## -1,labels=c('non-outlier','outlier'))
    c <- c+scale_alpha_discrete(range=c(.2,.6),labels=c('non-outlier','outlier'))+
        scale_color_viridis(discrete = T)+
        guides(alpha = guide_legend(override.aes = list(fill = 'black', lwd = 0)), color = 'none')
    c <- c+labs(x='Accession',y=paste0('Heterozygote frequency\n',whichtest),alpha='SNP')
    c
}
fig.hap <- function(data,type='recessive'){
    if (type=='recessive') {
        var <- 'gerpr'
        lstr <- labs(y='Load of haplotypes\nrecessive model', alpha='Type', x='Accession')
        agr <- aggregate(gerpr~race*type,hapsums,max) %>% rename(max=gerpr)
    }
    if (type=='additive') {
        var <- 'gerpa'
        lstr <- labs(y='Load of haplotype\nadditive model', alpha='Type', x='Accession')
        agr <- aggregate(gerpa~race*type,hapsums,max) %>% rename(max=gerpa)
    }
    adj <- lm(hapsums[,var]~hapsums$race+hapsums$type+hapsums$race*hapsums$type)
    cld <- cld.emmGrid(emmeans(adj,~race*type),Letters=letters,adjust="tukey") %>% mutate(lttr=gsub(" ","",.group,fixed=T))
    cld <- inner_join(cld,agr)
    ggplot(hapsums)+
        geom_boxplot(aes_string('race',var,alpha='factor(type,levels=c("LR","DH"))',fill='race'),outlier.alpha=1)+
        ## geom_text(data=cld,aes(race,max+120,group=factor(type,levels=c("LR","DH")),label=lttr),position=position_dodge(.75))+
        scale_fill_viridis(discrete = T)+
        scale_alpha_discrete(range = c(.2, .8))+
        guides(alpha = guide_legend(override.aes = list(fill = 'black')),
                      fill = 'none')+
        lstr
}



## draw plots
sampled <- stratify(freq.het)
c1 <- fig.het(sampled)

sampled <- stratify(prob.het)
c2 <- fig.het(sampled)

h1 <- fig.hap(hapsums,'recessive')

h2 <- fig.hap(hapsums,'additive')

## stitch together Fig 4
f4 <- plot_grid(h1+theme(legend.position = 'none'),
                c1+theme(legend.position = 'none'),ncol=1,labels='AUTO')
f4 <- plot_grid(f4,plot_grid(get_legend(h1),get_legend(c1),ncol=1),rel_widths=c(.8,.1))
f4

saveplot(f4,'fig4-a-recload-b-hetsfs')

## stitch together Fig 4 (alternative)
sourced <- T
source('scripts/gwas_anova.R')
mytraits <- c('grain_yield','Early_vigor','plant_heigh')

f4g <- gwasplt(filter(comb,trait%in%mytraits),nrow=1)
f4a <- plot_grid(
    c1+theme(legend.position = 'none'),
    f4g+theme(legend.position = 'none'),ncol=1,labels='AUTO')
f4a <- plot_grid(f4a,plot_grid(get_legend(f4g),ncol=1),rel_widths=c(.8,.1))
f4a

saveplot(f4a,'fig4-b-effects-a-hetsfs')

## stitch together Fig 4 Supplementary
f4s <- plot_grid(h2+theme(legend.position = 'none'),
                 c2+theme(legend.position = 'none'),ncol=1,labels='AUTO')
f4s <- plot_grid(f4s,plot_grid(get_legend(h2),get_legend(f4g),ncol=1),rel_widths=c(.8,.1))
f4s

saveplot(f4s,'fig4s-a-addload-b-hetprob')

## Fig Supplementary (Hets jProb only)
f4_c2 <- plot_grid(c2+theme(legend.position = 'none'),
                   get_legend(f4g),ncol=2,rel_widths=c(.8,.1))

saveplot(f4_c2,'fig4s-hetprob',height = 4)

#### GERP ANOVA and post hoc
a1 <- aov(hapsums$gerpr~hapsums$race+hapsums$type+hapsums$race*hapsums$type)
adj <- lm(hapsums$gerpr~hapsums$race+hapsums$type+hapsums$race*hapsums$type)
anova(adj)
TukeyHSD(a1)



