d = read.table('atac_20cell.sig.od.txt', header=F)

set.seed(2018)

ds = d[sample(dim(d)[1], 20000),c(-1,-2)]

colname_file = read.table('ct_list.18.colnames.txt', header=F)

ds = as.matrix(ds)

colnames(ds) = colname_file[,1]

library(pheatmap)

my_colorbar=colorRampPalette(c('white', 'dodgerblue'))(n = 128)
col_breaks = c(seq(0, 2000,length=33))

#png('atac_ccRE_sig.png')
pdf('atac_ccRE_sig.pdf')
pheatmap(ds, color=my_colorbar, cluster_cols = TRUE,cluster_rows=TRUE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=FALSE,show_colnames=TRUE)
dev.off()


