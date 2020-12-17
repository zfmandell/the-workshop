###Imports
require(ggplot2)
require(gridExtra)
require(grid)
require(plyr)
require(Matching)

WT_FPUTR<-merged$WT_total_FPUTR
WT_CDS<-merged$WT_total_CDS
deaD_FPUTR<-merged$deaD_total_FPUTR
deaD_CDS<-merged$deaD_total_CDS

avg_frame <- c(WT_FPUTR,deaD_FPUTR)
treatment_average <- factor( c(rep( 'a', length(WT_FPUTR)),rep( 'b', length(deaD_FPUTR))))
average_frame <- data.frame(treatment_average,avg_frame)

average<- ggplot(average_frame, aes(factor(treatment_average), avg_frame, fill=treatment_average))+
  theme(legend.position='none',
        axis.text.x = element_text(colour="grey20",size=22,face="plain"),
        axis.text.y = element_text(colour="grey20",size=22,face="plain"),
        axis.title.x = element_blank(),#element_text(colour="grey20",size=14,face="plain"),
        axis.title.y = element_text(colour="grey20",size=22,face="plain"),
        plot.title = element_text(hjust = 0.5,colour="grey20",size=24,face="plain"),
        panel.background = element_rect(fill = 'white'),
        axis.line = element_line(colour = "black", size = 1))+
  geom_boxplot(width = 0.2,outlier.shape=NA)+ylab('BP Per Feature')+
  ggtitle("Total Base Pairs in WT vs deaD KO Per Feature")+
  scale_x_discrete(labels = c("a" = "FPUTR WT","b" = "FPUTR deaD KO"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))





