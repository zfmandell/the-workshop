
###Imports
require(ggplot2)
require(gridExtra)
require(grid)
require(plyr)
require(Matching)

Az <- dA_WT_nonthresholded$`log2FC-Pause`
Bz <- dG_WT_nonthresholded$`log2FC-Pause`
Cz <- dA_dG_nonthresholded$`log2FC-Pause`

avg_frame <- c(Az,Bz,Cz)
treatment_average <- factor(c(rep( 'a', length(Az)),rep( 'b', length(Bz)),
                              rep( 'c', length(Cz))))
                              #rep( 'd', length(Dz)),
                              #rep( 'e', length(Ez)),rep( 'f', length(Fz)),
                              #rep( 'g', length(Gz)),rep( 'h', length(Hz)),rep( 'i', length(ctrl))))
                      
average_frame <- data.frame(treatment_average,avg_frame)

ggplot(average_frame, aes(factor(treatment_average), avg_frame))+
           theme_bw()+geom_hline(yintercept=2)+geom_hline(yintercept=-2)+
           geom_violin(fill='grey')+
           geom_point(average_frame[average_frame$avg_frame >= 2 | average_frame$avg_frame <= -2,],mapping = aes(y= avg_frame,col=threshold))+
           scale_color_manual(values = "red")+
            xlab('')+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
          ylab('log2 (WT Score)/(Mutant Score)')+
           theme(legend.position="none")+
  scale_x_discrete(labels = c("a" = 'NusA depletion\nvs.\nWT',"b" = 'nusG deletion\nvs.\nWT',
                              'c'='NusA depletion\nvs.\nnusG deletion'))
                             



                               #'d'='BG1030',
                              #'e'='BG546\nBG1030','f'='BG546\nBG838\nBG1030',
                              #'g'='BG1\nBG546\nBG838\nBG1030','h'='BG1\nBG838',
                              #'i'='all'))scale_fill_brewer(palette="Dark2")


