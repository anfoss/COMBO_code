
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


setwd("~/Desktop/Research/Postdoc_Krogan/Projects/Collaborations/COMBO/2023_Plasma/COMBO_plasma/ML_only")




### roc figure
dd = read.csv('output/roc_df.csv')
pl=ggplot(dd, aes(x=fpr, y=tpr, group=ids, col=ids))+geom_line(size=0.5, alpha=0.9)+
  geom_point(aes(x=0.3, y=0.9), colour="grey20", size=0.5)+
  scale_colour_viridis_d(option='turbo')+
  theme_Publication()+rmback+
  xlab('1-Specificity')+ylab('Sensitivity')+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+fix_plot+
  theme(legend.direction = 'vertical')+
  theme(legend.position = c(0.8,0.3))
pl
ggsave(pl= pl, filename = 'figures/roc_auc.pdf', height = 2.5, width = 2.5)


#### lassocoef


dd = read.csv('output/lassocoef.csv')
dd$kk = seq(1, length(dd$Pr_gn), by=1)
dd = subset(dd,dd$abs>0)
dd$gn = gsub(".*@","",dd$Pr_gn)
dd$sss = dd$abs
dd = head(dd, 10)
pl =ggplot(dd, aes(x=sss, y=reorder(gn, -sss)))+geom_bar(stat='identity', fill='#3C5488B2')+
  xlab('Abs(feature importance)LASSO')+theme_Publication()+rmback+fix_plot+
  scale_x_continuous(expand=c(0,0))+ylab('Top10 features')


ggsave(pl= pl, filename = 'figures/lassimpo.pdf', height = 2.5, width = 2.5)



#### sens at 90 spec

dd = read.csv('output/sens_70_spec.csv')
dd = subset(dd,dd$X2=='BestN')

pl =ggplot(dd, aes(x=X0, y=reorder(X1, -X0)))+geom_bar(stat='identity', fill='#3C5488B2')+
  xlab('Sensitivity at 70% specificity')+theme_Publication()+rmback+fix_plot+
  scale_x_continuous(expand=c(0,0))+ylab('Nr features')+geom_vline(xintercept = 0.9, linetype='dashed', col='grey40')
pl
ggsave(pl= pl, filename = 'figures/sens_spec.pdf', height = 2.5, width = 2)




#### pointrange for top 3-4-5 protein combination
df2 = read.csv('output/biosignature_test.csv')
p = ggplot(df2, aes(x=Pr_gn, y=mean, group=tb, col=tb), alpha=0.8) + 
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_pointrange(aes(ymin=mean-std, ymax=mean+std), position = position_dodge(width = 0.4), size=0.3)+
  theme_Publication()+ rmback+theme(axis.title.x = element_blank())+ylab('Unlikely TB-normalized intensity')+scale_colour_lancet()+theme(axis.text.x=element_text(size=6))
p = p+fix_plot+ theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
p
ggsave(pl= p, filename = 'figures/biosign.pdf', height = 2.5, width = 5)





#### test differences

df = read.csv('output/biosignature_test_no_grouping.csv')

my_comparisons <- list( c("Confirmed TB", "Unlikely TB"), c("Confirmed TB", "Unconfirmed TB"), c("Confirmed TB", "Healthy"), c('Confirmed TB', 'Latent TB') )

ggboxplot(df, "tb", "value", color = "tb", 
          add = "jitter", add.params = list(size = 0.2, alpha = 0.5), 
          facet.by = "Pr_gn", short.panel.labs = FALSE)+stat_compare_means(comparisons = my_comparisons)


