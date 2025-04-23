
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(ggridges)
library(viridis)

library(ComplexHeatmap)
library(circlize)
library(viridis)
library(ggbeeswarm)


theme_Publication<-function(base_size=9, base_family='sans') {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_family=base_family)
    + theme(plot.title = element_text(face = "plain",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            axis.title = element_text(face = "plain",size = 9),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(),
            axis.text = element_text(size=9), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing = unit(0, "cm"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="plain"),
            legend.text=element_text(size=9),
            legend.title=element_blank(),
            panel.border = element_rect(colour = "black", fill=NA, size=.5)
    ))
  
} 
rmback<-
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank()) +
  theme(strip.background =element_rect(fill="white"))

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


fix_plot = theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank())


## change me!
setwd('~/Desktop/Research/Postdoc_Krogan/Projects/Collaborations/COMBO/2023_Plasma/COMBO_plasma/figures/fig1')

##### Fig1B
id_table <- as.data.frame(cbind(c('Peptide', 'Protein'), c(7102, 860)))
id_table$V2=as.numeric(id_table$V2)

pl = ggplot(id_table, aes(x=V1, y=V2))+geom_bar(stat='identity')+
  geom_text(data = id_table, aes(x = V1, y = V2+200, label = V2), size=1.5)
pl = pl + ylab('# peptides/proteins') + 
  scale_y_continuous(expand=c(0,0))+
  theme_Publication()+
  rmback+ 
  scale_colour_viridis_d(option = 'viridis')+fix_plot+
  theme(axis.title.x=element_blank(), axis.text.x = element_text(size=7))+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
pl
ggsave('Fig1_B.pdf', pl, width = 1.5, height = 1.5)


### Fig1C
dd = read.csv('data/info_pep_level_run.csv')
dd = subset(dd, dd$type!='inte')
dd2= subset(dd,dd$type=='pep')
dd3 = subset(dd, dd$type!='pep')
dd = rbind(dd3, dd2)
pepdd = subset(dd, dd$type=='pep') 
median(pepdd$max)
pl = ggplot(pepdd, aes(x=rid, y=max, col=type, group=type))+geom_line(size=0.5, alpha=0.75)
pl = pl + xlab('Injection nr')+ylab('# Peptides')
pl = pl + theme_Publication()+rmback
pl =pl +  scale_x_continuous(expand = c(0,0)) +theme(axis.line = element_line(colour = "black"),
                                                                                                 panel.grid.major = element_blank(),
                                                                                                 panel.grid.minor = element_blank(),
                                                                                                 panel.border = element_blank(),
                                                                                                 panel.background = element_blank())
pl = pl+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+theme(legend.position='none') + geom_hline(yintercept = 1708)
pl
ggsave(plot=pl, filename = 'Fig1_C.pdf', width = 1.5, height = 1.5)

## Fig1D
pepdd = subset(dd, dd$type=='prot') 
mean(pepdd$max)

pl = ggplot(pepdd, aes(x=rid, y=max, col=type, group=type))+geom_line(size=0.5, alpha=0.75)
pl = pl + xlab('Injection nr')+ylab('# Proteins')
pl = pl + theme_Publication()+rmback
pl =pl +  scale_x_continuous(expand = c(0,0)) +theme(axis.line = element_line(colour = "black"),
                                                     panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(),
                                                     panel.border = element_blank(),
                                                     panel.background = element_blank())
pl = pl+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+theme(legend.position='none')
ggsave(plot=pl, filename = 'Fig1_D.pdf', width = 1.5, height = 1.5)

## FIg1E
df = read.csv('data/completness.csv')
pl3 = ggplot(df, aes(x=idx, y=cmplt))+geom_line(aes(group=dummy), col='#3C5488B2')+
  ylab('Data completeness')+
  xlab('Protein number')+theme_Publication()+
  rmback+fix_plot+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
pl3

ggsave(plot = pl3, filename = 'Fig1_E.pdf', width = 1.5, height = 1.5)


## Fig1 F
dd = read.csv('test_rank_mol.csv')
dd$cc = dd$conc * dd$unit
dd$cc = log10(dd$cc)
pl = ggplot(dd, aes(x=cc))+geom_density(aes(fill=detected), alpha=0.75)+
  scale_fill_viridis_d()+theme_Publication()+rmback+fix_plot+
  scale_x_continuous(expand=c(0,0))+
  scale_y_continuous(expand=c(0,0))+
  ylab('Density')+xlab('Log10(ng/L)')

ggsave('Fig1_Fpdf', height = 2, width = 2, plot=pl)


##FIg2E
df = read.csv('data/cv_all.csv')
df$X0=df$X0*100
df$loc <- factor(df$loc, levels=c('Gambia ', 'Peru ', 'PRO South Africa ', 'Uganda ', 'Total'))

pl= ggplot(df, aes(x=loc, y=X0, group=loc, color=loc), alpha=0.25)+
  geom_beeswarm(cex = 0.25, alpha=0.25) + ylab('Protein-level CV %')+fix_plot+
  scale_color_manual(values = c(
    "Gambia " = "#8FAADC",     # original blue from AAAS palette
    "Peru " = "#A20056FF",            # changed to pink
    "PRO South Africa " = "#008B45FF", # original red
    "Uganda " = "#631879FF",     # original purple
    "Total " = "#808180FF"            # changed to grey
  ))+
  theme(legend.position = 'none')+theme(axis.title.x=element_blank())+ theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

pl+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
ggsave(plot = pl, filename = 'Fig2_E.pdf', width = 2.5, height = 2)



## Fig2D
dd = read.csv('data/pca_input.csv')
pl = ggplot(dd, aes(x=PC1, y=PC2, col=Loc))+geom_point(size=.5, alpha=.75)+
  scale_color_manual(values = c(
    "Gambia " = "#8FAADC",     # original blue from AAAS palette
    "Peru " = "#A20056FF",            # changed to pink
    "South Africa " = "#008B45FF", # original red
    "Uganda " = "#631879FF"     # original purple
              # changed to grey
  ))+
  theme_Publication()+fix_plot+rmback+theme(legend.direction = 'vertical', legend.position = 'right')
ggsave(plot= pl,filename =  'Fig2_D.pdf', height = 3, width = 4)



## revision


### order per run
dd = read.csv('../../meta/info_pep_level_run_no_dedup.csv')
dd = subset(dd, dd$type!='inte')

pepdd = subset(dd, dd$type=='pep') 
median(pepdd$X0)
pl = ggplot(pepdd, aes(x=run_nr, y=X0, col=type, group=type))+geom_line(size=0.5, alpha=0.75)
pl = pl + xlab('Injection nr')+ylab('# Peptides')
pl = pl + theme_Publication()+rmback
pl =pl +  scale_x_continuous(expand = c(0,0)) +theme(axis.line = element_line(colour = "black"),
                                                     panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(),
                                                     panel.border = element_blank(),
                                                     panel.background = element_blank())
pl1 = pl+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+theme(legend.position='none')
pl1

## proteins
pepdd = subset(dd, dd$type=='prot') 
mean(pepdd$max)

pl = ggplot(pepdd, aes(x=run_nr, y=X0, col=type, group=type))+geom_line(size=0.5, alpha=0.75)
pl = pl + xlab('Injection nr')+ylab('# Proteins')
pl = pl + theme_Publication()+rmback
pl =pl +  scale_x_continuous(expand = c(0,0)) +theme(axis.line = element_line(colour = "black"),
                                                     panel.grid.major = element_blank(),
                                                     panel.grid.minor = element_blank(),
                                                     panel.border = element_blank(),
                                                     panel.background = element_blank())
pl2 = pl+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+theme(legend.position='none')
pl2