DE <- condition_treated_results
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
  xlim(-20,20)+
  ylim(0,350)+
  theme_minimal()+
  theme(legend.position="none")
  
