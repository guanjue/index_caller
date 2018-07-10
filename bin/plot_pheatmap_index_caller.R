####################################################
### get parameters
args = commandArgs(trailingOnly=TRUE)
signal_matrix_file = args[1]
output_filename = args[2]
signal_input_list = args[3]
signal_matrix_start_col = args[4]
high_color = args[5]
low_color = args[6]
log2 = args[7]
smallnum = as.numeric(args[8])

####################################################
### use pheatmap to plot heatmaps
color_heatmap = function(color_matrix, high_color, low_color, format, outputname){
	library(pheatmap)
	format(outputname, width = 1000, height = 1000) ### output name
	par() ### set heatmap margins
	### plot pheatmap
	my_colorbar=colorRampPalette(c(low_color, high_color))(n = 128)
	col_breaks = c(seq(0, 2000,length=33))
	pheatmap(color_matrix, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
	dev.off()
}

#################################################### 
############ read input files
####################################################
### read signal matrix file
signal_matrix_od = as.matrix(read.table(signal_matrix_file, header=FALSE))
### extract signal matrix without info
signal_matrix = signal_matrix_od[ , signal_matrix_start_col:dim(signal_matrix_od)[2] ]
rownames(signal_matrix) = signal_matrix_od[,1]
### convert to numeric matrix
class(signal_matrix) = 'numeric'

### log2 transform
if (log2=='T'){
	signal_matrix = log2(signal_matrix+smallnum)
}

###### read colnames file
colname_file = read.table(signal_input_list, header=F)
### add colnames
print(dim(signal_matrix))
colnames(signal_matrix) = colname_file[,1]


format = png
color_heatmap(signal_matrix, high_color, low_color, format, output_filename)


signal_vec = as.vector(signal_matrix)

library(mixtools)

mixmdl = normalmixEM(signal_vec,k = 3)
plot(mixmdl,which=2)
lines(density(signal_vec), lty=2, lwd=2)


post = apply(mixmdl$posterior, 1, function(x) which.max(x)-1)


post_matrix = t(matrix(post, nrow =  dim(signal_matrix)[2], byrow = TRUE))

level1_lim = 2.9
level2_lim = 6

post_matrix = (signal_matrix >=level1_lim )*1

post_matrix[signal_matrix>=level2_lim] =2


new_label = apply(post_matrix, 1, function(x) paste(x, collapse="_"))


signal_matrix_relabel = signal_matrix
rownames(signal_matrix_relabel) = new_label

signal_matrix_relabel = signal_matrix_relabel[order(rownames(signal_matrix_relabel)),]
color_heatmap(signal_matrix_relabel, high_color, low_color, format, paste(output_filename, '.newlabel.png', sep=''))

png(paste(output_filename, 'density.png', sep=''))
plot(density(signal_vec, bw=0.2))
abline(v = level1_lim, col='orange', lty = 2, lwd = 1.5)
abline(v = level2_lim, col='red', lty = 2, lwd = 1.5)
dev.off()


label2newlabel = cbind(signal_matrix_od[,1], new_label)
write.table(label2newlabel, paste(output_filename, '.label2newlabel.txt', sep=''), quote=FALSE, sep='\t', row.names = FALSE, col.names = FALSE)




