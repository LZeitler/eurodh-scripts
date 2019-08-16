source('source_me.R')

#### Supplementary figure S4: PCA with PC1, PC2, PC3

gdsfile.co.red <- 'data/260419_step6_red_v4.gds'

co <- snpgdsOpen(gdsfile.co.red)        # reduced

pca.co <- snpgdsPCA(co)

pca.co.tab <- data.frame(sample.id = pca.co$sample.id,
                         Set = as.factor(substr(pca.co$sample.id,1,5)),
                         Accession = as.factor(substr(pca.co$sample.id,4,5)),
                         Type = as.factor(substr(pca.co$sample.id,1,2)),
                         ev1 = pca.co$eigenvect[,1],
                         ev2 = pca.co$eigenvect[,2],
                         ev3 = pca.co$eigenvect[,3],
                         stringsAsFactors = F)

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

p12 <- pcaplot(pca.co.tab,1,2,x = ev1, y = ev2, color = Accession, shape=Type)
p13 <- pcaplot(pca.co.tab,1,3,x = ev1, y = ev3, color = Accession, shape=Type)
p23 <- pcaplot(pca.co.tab,2,3,x = ev2, y = ev3, color = Accession, shape=Type)
pcabox <- plot_grid(p12+theme(legend.position = "none"),
                    p13+theme(legend.position = "none"),
                    p23+theme(legend.position = "none"),
                    labels='AUTO')

pcaleg <- get_legend(p12)

pcaboxes <- plot_grid(pcabox,pcaleg,rel_widths = c(.9,.1))

saveplot(pcaboxes,'pca_boxes')

snpgdsClose(co)

