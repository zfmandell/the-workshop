
###Imports
require(ggplot2)
require(gridExtra)
require(grid)
require(plyr)
require(Matching)

SIz <- SI_build$`dG-Hairpin (kcal/mol)`
Az <- A_built[A_built$`d%T` >= 25,]$`dG-Hairpin (kcal/mol)`
Gz <- G_built[G_built$`d%T` >= 25,]$`dG-Hairpin (kcal/mol)`
Rz <- R_built[R_built$`d%T` >= 25,]$`dG-Hairpin (kcal/mol)`

  
avg_frame <- c(SIz,Az,Gz,Rz)
treatment_average <- factor(c(rep( 'a', length(SIz)),rep( 'b', length(Az)),
                              rep( 'c', length(Gz)),rep( 'd', length(Rz))))
                            
average_frame <- data.frame(treatment_average,avg_frame)

ggplot(average_frame, aes(factor(treatment_average), avg_frame,fill=treatment_average))+
           theme_bw()+scale_y_reverse()+
           geom_violin()+
           geom_boxplot(alpha=0.2)+
            xlab('')+
            theme_classic()+
          ylab('ΔG (kcal/mol)')+
           theme(legend.position="none")+
  scale_x_discrete(labels = c("a" = 'SI',"b" = 'A',
                              'c'='G','d'='R'))
                            
                             



  

