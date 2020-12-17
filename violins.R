
###Imports
require(ggplot2)
require(gridExtra)
require(grid)
require(plyr)
require(Matching)


total<-Total$log2FC
#dependent<-log2(Dependent$`WT Score`)
#independent<-log2(Independent$`WT Score`)
strongest<-Strongest$log2FC

avg_frame <- c(strongest,total)
treatment_average <- factor( c(rep( 'a', length(strongest)),rep( 'b', length(total))))
average_frame <- data.frame(treatment_average,avg_frame)
#
ggplot(average_frame, aes(factor(treatment_average), avg_frame, fill=treatment_average))+
           theme_bw()+
           scale_fill_manual(values=c("#6D89D5","#FFE573"))+
           geom_boxplot()+ylab('Log2(FC)')+
           ggtitle('Impact of Mutation')+
           xlab('')+
           theme(legend.position="none")+
           scale_x_discrete(labels = c("a" = "100 Strongest in WT","b" = "Total"))


