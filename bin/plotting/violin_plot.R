#library(ggplot2)
#library(vioplot)
library(mvtnorm)

### get parameters
args = commandArgs(trailingOnly=TRUE)
index_set_inputfile = args[1]
regin_index_inputfile = args[2]
regin_signal_inputfile = args[3]

index_set_all_heatmap = args[4]
index_set_thresh_heatmap = args[5]
index_all_heatmap = args[6]


### read index set signal matrix
read_enriched_index_set_matrix = function(inputfile){
	data_index_set = as.matrix(read.table(inputfile, header=T))
	data_index_set_binary_patter = data_index_set[,1]
	rownames(data_index_set) = data_index_set_binary_patter
	data_index_set = data_index_set[,-1]
	class(data_index_set) = "numeric" 
	return(c(data_index_set, data_index_set_binary_patter))
}


### read index matrix
read_index_matrix = function(inputfile){
	### read matrix as string matrix
	data_index_matrix = read.table(inputfile, header=T, sep='\t', colClasses = "character")
	rownames(data_index_matrix) = data_index_matrix[,1]
	data_index_matrix = data_index_matrix[,-1]
	### collapse string vector
	data_index_matrix = apply(data_index_matrix, MARGIN=1, FUN=function(x) paste(x, collapse='_') )
	return(data_index_matrix)
}

### read index signal matrix
read_signal_matrix = function(inputfile){
	### read matrix as string matrix
	data_signal_matrix = as.matrix(read.table(inputfile, header=T, sep='\t'))
	rownames(data_signal_matrix) = data_signal_matrix[,1]
	data_signal_matrix = data_signal_matrix[,-1]
	return(data_signal_matrix)
}

### get index binary label matrix
index_label_matrix = read_index_matrix('celltype.index.sorted.txt')

### get index signal matrix
index_sig_matrix = read_signal_matrix('celltype.index.tpm.sorted.txt') 

### get index_set binary labels
index_set_filtered_matrix = read_signal_matrix('celltype.index_set_filtered.sorted.txt')
### to numeric matrix
class(index_set_filtered_matrix) = "numeric"
index_set_filtered_label = rownames(index_set_filtered_matrix)

### get multiple variable normal distribution from previous peak calling result
### get new index calling probability

### Previous matrix (plus 0.01 then log2 transform)
index_sig_matrix_log2 = log2(index_sig_matrix+0.01)

index_sig_matrix_p_call = c()
enriched_index_set_position = rep(0, dim(index_sig_matrix_log2)[1])

### enriched index set:
for (i in c(1: length(index_set_filtered_label))){
	###
	print(paste('index set', toString(i), index_set_filtered_label[i], sep=':'))
	### get index_sig_i data (plus 0.01 then log2 transform)
	index_sig_i_position = index_label_matrix==index_set_filtered_label[i]
	index_sig_i = log2(index_sig_matrix[index_sig_i_position, ] + 0.01)

	### save enriched id positions
	enriched_index_set_position = enriched_index_set_position + index_sig_i_position

	### boxplot index_set cell type signals
	png(paste('index_set_boxplot/index_set_boxplot.', toString(i), '.', index_set_filtered_label[i], '.png', sep=''))
	boxplot(index_sig_i, ylim=c(min(index_sig_matrix_log2), quantile(index_sig_matrix_log2, probs=0.99)))
	dev.off()

	### get index_set mean vector & correlation matrix 
	x_mean = colMeans(index_sig_i)
	x_cor = cor(index_sig_i)

	### read mvnorm function for matrix
	pmvnorm_all = function(vector){
		x_mean_used <<- x_mean
		x_cor_used <<- x_cor
		p = pmvnorm(mean = x_mean_used, corr = x_cor_used, upper=vector, lower=rep(-Inf, length(x_mean)))
		return(p[1])
	}

	### calculate the probability based on pmvnorm
	print('start pmvnorm')
	index_sig_matrix_i_p = apply(index_sig_matrix_log2, MARGIN=1, FUN=function(x) pmvnorm_all(x))

	###
	index_sig_matrix_p_call = cbind(index_sig_matrix_p_call, index_sig_matrix_i_p)
}

### NOT enriched index set:
for (i in c(1)){
	###
	print('NOT enriched index set:')
	### get NOT enriched index set: enriched_index_set_position==0 
	### (plus 0.01 then log2 transform)
	index_sig_i = log2(index_sig_matrix[enriched_index_set_position==0, ] + 0.01)

	### boxplot index_set cell type signals
	png(paste('index_set_boxplot/index_set_boxplot.', toString(i), '.', 'NOT_enriched', '.png', sep=''))
	boxplot(index_sig_i, ylim=c(min(index_sig_matrix_log2), quantile(index_sig_matrix_log2, probs=0.99)))
	dev.off()

	### get index_set mean vector & correlation matrix 
	x_mean = colMeans(index_sig_i)
	x_cor = cor(index_sig_i)

	### 
	index_sig_matrix_i_p = c()
	for ( j in c(1: dim(index_sig_matrix_log2)[1]) ){
		if (i%%1000==0){
			print(j)
		}
		p = pmvnorm_all(index_sig_matrix_log2[j,], x_mean, x_cor)
		index_sig_matrix_i_p[j] = p[1]
	}

	###
	index_sig_matrix_p_call = cbind(index_sig_matrix_p_call, index_sig_matrix_i_p)
}


print(dim(index_sig_matrix_p_call))
#print(head(index_sig_matrix_p_call))

colnames(index_sig_matrix_p_call) = c(1: dim(index_sig_matrix_p_call)[2])
rownames(index_sig_matrix_p_call) = rownames(index_sig_matrix)

### the sum of probability is standardized to 1
print('get index call probability matrix')
index_sig_matrix_p_call = t( apply(index_sig_matrix_p_call, MARGIN=1, FUN=function(x) x/sum(x)) )
### write all index set call & probability 
write.table(index_sig_matrix_index_set_call, 'index_sig_matrix_index_set_call.txt', quote=F, sep='\t', row.names = TRUE, col.names = NA)

### get belonging index set id
print('get index call index_id & probability')
index_sig_matrix_index_set_call = t( apply((index_sig_matrix_p_call), MARGIN=1, FUN=function(x) c(which.max(x), max(x))) )
### write all index set probility matrix 
write.table(index_sig_matrix_p_call, 'index_sig_matrix_p_call.txt', quote=F, sep='\t', row.names = TRUE, col.names = NA)



############################################
print('sort index based on new index call')
### sort input signal matrix based on the new index set call
index_sig_matrix_sort = index_sig_matrix[order(index_sig_matrix_index_set_call[i,1]), ]
### write all index set call & probability 
write.table(index_sig_matrix_sort, 'index_sig_matrix_sort.txt', quote=F, sep='\t', row.names = TRUE, col.names = NA)



############################################
print('get new index set numbers:')
### recalculate index set DNA region number
index_set_num = table(index_sig_matrix_index_set_call[i,1])

### change previous index set number to new called number
#index_set_filtered_matrix = index_set_filtered_matrix[c(1:15),]
index_set_filtered_matrix_recalnum = c()
for (i in c(1: dim(index_set_filtered_matrix)[1])){
	index_set_filtered_matrix_vec = index_set_filtered_matrix[i, ]

	### recalculate the number of DNA region in index set
	index_set_filtered_matrix_vec_new = index_set_filtered_matrix_vec / max(index_set_filtered_matrix_vec) * index_set_num[i]

	index_set_filtered_matrix_recalnum = cbind(index_set_filtered_matrix_recalnum, index_set_filtered_matrix_vec_new)
}

### add colnames & rownames
colnames(index_set_filtered_matrix_recalnum) = colnames(index_set_filtered_matrix)
rownames(index_set_filtered_matrix_recalnum) = rownames(index_set_filtered_matrix)

### write index_set_filtered_matrix_recalnum
write.table(index_set_filtered_matrix_recalnum, 'celltype.index_set_filtered.sorted.recalnum.txt', quote=F, sep='\t', row.names = TRUE, col.names = NA)







