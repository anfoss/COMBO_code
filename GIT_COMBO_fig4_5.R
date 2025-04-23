
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



#### Fig 4A
dd = read.csv('new_biosignature/imp_coef_lasso.csv')
dd$kk = seq(1, length(dd$Pr_gn), by=1)
dd = subset(dd,dd$abs>0)
dd$gn = gsub(".*@","",dd$Pr_gn)
dd$sss = dd$abs
dd = head(dd, 10)
pl =ggplot(dd, aes(x=sss, y=reorder(gn, -sss)))+geom_bar(stat='identity', fill='grey50')+
  xlab('Abs(feature importance)LASSO')+theme_Publication()+rmback+fix_plot+
  scale_x_continuous(expand=c(0,0))+ylab('Top10 features')

pl
ggsave(pl= pl, filename = 'revision_figures/Fig4_A.pdf', height = 2.5, width = 2.5)


### Fig 4 B
dd = read.csv('new_biosignature/20250423_ROC_CV.csv')
pl=ggplot(dd, aes(x=fpr, y=tpr, group=clf, col=clf))+geom_line(size=0.5, alpha=0.9)+
  geom_point(aes(x=0.3, y=0.9), colour="grey20", size=0.5)+
  scale_colour_viridis_d(option='turbo')+
  theme_Publication()+rmback+
  xlab('1-Specificity')+ylab('Sensitivity')+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+fix_plot+
  theme(legend.direction = 'vertical')+
  theme(legend.position = c(0.8,0.3))
pl
ggsave(pl= pl, filename = 'revision_figures/Fig4_B.pdf', height = 3, width = 3)




#### Fig 4 C
dd = read.csv('new_biosignature/biosignature_output.csv')
dd$nfeat = as.character(dd$nfeat)
pl =ggplot(dd, aes(x=prec70sens, y=nfeat))+geom_bar(stat='identity', fill='grey50')+
  xlab('% Sensitivity at 70% specificity')+theme_Publication()+rmback+fix_plot+
  scale_x_continuous(expand=c(0,0))+ylab('Nr features')+geom_vline(xintercept = 0.9, linetype='dashed', col='darkred')
pl
ggsave(pl= pl, filename = 'revision_figures/Fig4_C.pdf', height = 3, width = 3)




#### Fig4E
df2 = read.csv('new_biosignature/biosignature_clinical_sites.csv')
p = ggplot(df2, aes(x=Pr_gn, y=mean, group=tb, col=tb), alpha=0.8) + 
  geom_hline(yintercept = 0, linetype='dashed')+
  geom_pointrange(aes(ymin=mean-std, ymax=mean+std), position = position_dodge(width = 0.4), size=0.3)+
  theme_Publication()+ rmback+theme(axis.title.x = element_blank())+ylab('Unlikely TB-normalized intensity')+scale_colour_lancet()+theme(axis.text.x=element_text(size=6))
p = p+fix_plot+ theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
p
ggsave(pl= p, filename = 'revision_figures/Fig4_E.pdf', height = 2.5, width = 5)


## SupFig


dd = read.csv('new_biosignature/topN_vs_bestN.csv')
dd = subset(dd,dd$n!=1)
pl=ggplot(dd, aes(x=fpr, y=tpr, group=t, col=t))+geom_line(size=0.5, alpha=0.9)+
  geom_point(aes(x=0.3, y=0.9), colour="grey20", size=0.5)+
  scale_color_aaas()+
  theme_Publication()+rmback+
  xlab('1-Specificity')+ylab('Sensitivity')+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+fix_plot+
  facet_wrap(~n, nrow = 2, scales = 'free_x')+
  theme(legend.position=c(0.8,0.3), legend.direction = 'vertical')+
  theme(panel.spacing.x=unit(2, "lines"))
pl
ggsave(pl= pl, filename = 'revision_figures/SupFig_topn_best_n.pdf', height = 5, width = 6)


## supfig individual AUCs
dd = read.csv('new_biosignature/individual_AUC.csv')
pl=ggplot(dd, aes(x=fpr, y=tpr, group=clf, col=clf))+geom_line(size=0.5, alpha=0.9)+
  geom_point(aes(x=0.3, y=0.9), colour="grey20", size=0.5)+
  scale_colour_aaas()+
  geom_abline(slope = 1, intercept = 0, linetype='dashed', col='grey70')+
  theme_Publication()+rmback+
  xlab('1-Specificity')+ylab('Sensitivity')+
  scale_y_continuous(expand=c(0,0))+
  scale_x_continuous(expand=c(0,0))+fix_plot+
  theme(legend.direction = 'vertical')+
  theme(legend.position = c(0.8,0.3))
pl
ggsave(pl= pl, filename = 'revision_figures/SupFig_individual_AUCs.pdf', height = 4.5, width = 4.5)


### Fig 5A
dd = read.csv('new_biosignature/pred_count.csv')
dd$cc = gsub('model', '', dd$cc)
dd$X0 = gsub('0', 'Negative', dd$X0)
dd$X0 = gsub('1', 'Positive', dd$X0)
pl =ggplot(dd, aes(x=cc, y=nr, fill=X0))+geom_bar(stat='identity', alpha=0.75)+
  ylab('# Unconfirmed TB')+theme_Publication()+rmback+fix_plot+
  scale_y_continuous(expand=c(0,0))+xlab('N-protein models')+scale_fill_viridis_d()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  geom_text(aes(label = nr), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 3) +
  theme(legend.position = 'right', legend.direction = 'vertical')+theme(axis.title.x = element_blank())
pl
ggsave(pl= pl, filename = 'revision_figures/Fig5A.pdf', height = 2.5, width = 3)

### Fig 5C


dd = read.csv('new_biosignature/lda_projection.csv')

dd$class = as.character(dd$class)
pl = ggplot(dd, aes(x=X0, y=X1, col=class))+geom_point(size=0.5, alpha=0.75)+xlab('PC1')+ylab('PC2')+theme_Publication()+rmback+
  scale_color_colorblind()+guides(col=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_colorblind()+
  stat_ellipse(geom = "polygon",
               aes(fill = class), 
               alpha = 0.1)

pl
ggsave(filename = 'revision_figures/Fig5C_ellipses.pdf', pl, height = 3, width = 3)


dd = read.csv('new_biosignature/lda_projection.csv')
dd$class = as.character(dd$class)
pl = ggplot(dd, aes(x=X0, y=X1, col=class))+geom_point(size=0.5, alpha=0.75)+xlab('PC1')+ylab('PC2')+theme_Publication()+rmback+
  scale_color_colorblind()+guides(col=guide_legend(nrow=2,byrow=TRUE))+
  scale_fill_colorblind()

pl
ggsave(filename = 'revision_figures/Fig5C.pdf', pl, height = 3, width = 3)


### SupFig XXX
dd = read.csv('new_biosignature/lda_loadings.csv')
dd$X0= abs(dd$X0)
p = ggplot()+geom_point(data = subset(dd,dd$biosignature=='False'), aes(x=rank, y=X0), col='#808080', size=0.5, alpha=0.30)
p = p+geom_point(data=subset(dd,dd$biosignature=='True'), aes(x=rank, y=X0), col='#8B0000', size=1, alpha=1)+
  geom_text_repel(data=subset(dd,dd$biosignature=='True'),aes(x=rank, y=X0,label=as.character(Gene)), hjust = 2, max.overlaps = 100, col = 'black',size = 3)
p
p4 = p+ theme_Publication()+
  rmback+fix_plot+theme(legend.position='none')+ylab('Absolute Feature Loading PC1')+xlab('Rank')
ggsave(filename = 'revision_figures/SupFig_lda_loadings.pdf', p4, height = 3, width = 3)


## SupFig XXX
dd = read.csv('new_biosignature/healthy_latent_pca.csv')
pl = ggplot(dd, aes(x=X0, y=X1, col=class))+geom_point(size=1, alpha=1)+
  facet_grid(~nrfreat)+xlab('PC1')+ylab('PC2')+theme_Publication()+rmback+
  scale_color_aaas()+scale_fill_aaas()+
  stat_ellipse(geom = "polygon",
               aes(fill = class), 
               alpha = 0.1)
pl
ggsave(filename = 'revision_figures/SupFig_pca_healthy_TB.pdf', pl, height = 3, width = 6)




