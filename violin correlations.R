require(ggplot2)
require(gridExtra)
require(grid)
require(plyr)
require(Matching)


A_B<-compiled_bottom_ind_A$hairpin_deltaG
A_T<-compiled_top_ind_A$hairpin_deltaG
A_S<-compiled_strong_A$hairpin_deltaG
G_B<compiled_bottom_ind_G$hairpin_deltaG
G_T<-compiled_top_ind_G$hairpin_deltaG
G_S<-compiled_strong_G$hairpin_deltaG

avg_frame <- c(A_B,A_T,A_S,G_B,G_T,G_S)
treatment_average <- factor( c(rep( 'a',50),rep( 'b',50) ,rep( 'c',50),rep( 'd',50),rep( 'e',50),rep( 'f',50))
average_frame <- data.frame(treatment_average,avg_frame)


ggplot(average_frame, aes(factor(treatment_average), avg_frame, fill=treatment_average))+
    geom_boxplot(alpha=0.5,outlier.size = 0, coef = 0)+
    theme(legend.position='none',
        axis.text.x = element_text(colour="grey20",size=22,face="plain"),
        axis.text.y = element_text(colour="grey20",size=22,face="plain"),
        axis.title.x = element_blank(),#element_text(colour="grey20",size=14,face="plain"),
        axis.title.y = element_text(colour="grey20",size=22,face="plain"),
        plot.title = element_text(hjust = 0.5,colour="grey20",size=24,face="plain"),
        panel.background = element_rect(fill = 'white'),
        axis.line = element_line(colour = "black", size = 1))+
    ylab('ΔG of Hairpin')+
    ggtitle('Energetic Stability of Hairpin')+
    scale_x_discrete(labels = c("a" = "Bottom NusA Independent","b" = "Top NusA Independent","c"='Most NusA Dependent',"d"='Bottom NusG Independent','e'='Top NusG Independent','f'='Most NusG Dependent'))



