library(ggplot2)
DE <- WT_dG_minusAr_dropped_results
DE$threshold = as.factor(abs(DE$log2FoldChange) >= 2 & DE$padj <= 0.005)

ggplot(DE, aes(log2FoldChange, -log10(padj)))+
geom_point(aes(col=threshold,alpha = 0.005))+
scale_color_manual(values=c("black", "red"))+
ylab('-log10(FDR)')+
xlab('log2(FC)')+
theme(legend.position="none")+
geom_hline(yintercept=2.3) +
geom_vline(xintercept = -2) +
geom_vline(xintercept = 2) +
theme_minimal()+
ylim(0,10)+
xlim(-10,10)+
theme(legend.position="none")+
ggtitle("WT vs dG + Ar, abs(FC) >= 4, padj<=0.005")




