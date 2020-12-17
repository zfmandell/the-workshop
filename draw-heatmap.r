library(gplots)

# Read normalized counts
data = deaD_DE_max
gene = data[,1]
vals = as.matrix(data[,3:ncol(data)])

# Adds a little noise to each element
# To avoid the clusteing function failing on zero
# variance datalines.
vals = jitter(vals, factor = 1, amount=0.00001)


# Calculate zscore
score = NULL
for (i in 1:nrow(vals)) {
    row=vals[i,]
    zscore=(row-mean(row))/sd(row)
    score =rbind(score,zscore)
}

row.names(score) = gene$id
zscore=score

# Generate heatmap
mat = as.matrix(zscore)
pheatmap(mat)
# Opent the drawing device.
#pdf('|cat')

#colors = colorRampPalette(c("green","black","red"),space="rgb")(256)
#heatmap.2(mat,col=colors,density.info="none",trace="none", margins=c(7,7),lhei=c(1,5))
#dev.off()
