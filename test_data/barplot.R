counts = read.table('new_ct_pk_num.txt', header=F)


rownames_counts = counts[,2]

counts_c = as.matrix(counts[,1])

rownames(counts_c) = rownames_counts
#counts_c = counts_c[order(counts_c),]

png('barplot.counts.new.png', width=2000, height=2000)
barplot(t(counts_c), ylim=c(0,165000))
dev.off()




counts = read.table('old_ct_pk_num.txt', header=F)
rownames_counts = counts[,2]

counts_c = as.matrix(counts[,1])

rownames(counts_c) = rownames_counts
#counts_c = counts_c[order(counts_c),]

png('barplot.counts.old.png', width=2000, height=2000)
barplot(t(counts_c), ylim=c(0,165000))
dev.off()
