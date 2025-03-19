# complicated plot

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
library(ggExtra)
library(corrplot)


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



fix_plot = theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank())

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


data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}



fix_plot = theme(axis.line = element_line(colour = "black"),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.border = element_blank(),
                 panel.background = element_blank())



new_lines_adder = function(test.string, interval) {
  #split at spaces
  string.split = strsplit(test.string," ")[[1]]
  # get length of snippets, add one for space
  lens <- nchar(string.split) + 1
  # now the trick: split the text into lines with
  # length of at most interval + 1 (including the spaces)
  lines <- cumsum(lens) %/% (interval + 1)
  # construct the lines
  test.lines <- tapply(string.split,lines,function(line)
    paste0(paste(line,collapse=" "),"\n"),simplify = TRUE)
  # put everything into a single string
  result <- paste(test.lines,collapse="")
  return(result)
}
add_newlines = function(x, interval) {
  
  # make sure, x is a character array   
  x = as.character(x)
  # apply splitter to each
  t = sapply(x, FUN = new_lines_adder, interval = interval,USE.NAMES=FALSE)
  return(t)
}


setwd("~Desktop/Research/Postdoc_Krogan/Projects/Collaborations/COMBO/2023_Plasma/COMBO_plasma/figures/fig3")




###### volcano plot
dd = read.csv('data/Confirmed_vs_Unlikely.csv')

mycol = c('F'='#808080','T'='#21908c')

dd$issign = ifelse(dd$q<0.05, 'T','F')
p = ggplot(dd, aes(y=logq, x=Log2FC, col=issign))+geom_point(size=1)+xlab('Log2(Confirmed TB/Unlikely TB)')+
  ylab('-Log10(BH-adjusted) p')+
  geom_hline(yintercept = -log10(0.05), linetype='dashed', col='grey10', size=.5)+
  scale_colour_manual(values=mycol)+
  theme_Publication()+geom_vline(xintercept = 0,  linetype='dashed', col='grey10', size=.5)+
  rmback+fix_plot+theme(legend.position='none')+theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))
p

ggsave(filename = 'volcano.pdf', plot = p, height = 1.5, width = 1.5)





#### barplot up and down nr

### up 
mycol = c('Up'='#a02c2cff','Down'='#0088aaff')


dd = as.data.frame(c(18,30))
dd$dd = c('Up', 'Down')
colnames(dd) = c('a', 'b')
p = ggplot(dd, aes(y=a, x=b))+geom_bar(stat='identity', aes(fill=b))+theme_Publication()+
  rmback+fix_plot+ rmback+ fix_plot+
  scale_y_continuous(expand=c(0,0))+theme(axis.title.y=element_blank(), legend.title = element_blank()) + xlab('# DEP')
p =  p + theme(plot.margin=grid::unit(c(0,0,0,0), "mm"))+scale_fill_manual(values=mycol)+theme(legend.position='none')
ggsave(plot=p,filename =  'dep_count.pdf', width = 1.25, height = 1.25)
