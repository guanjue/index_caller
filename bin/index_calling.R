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

### read index signal matrix
read_wg_signal_matrix = function(inputfile){
	### read matrix as string matrix
	data_signal_matrix = (read.table(inputfile, header=F, sep='\t'))
	#rownames(data_signal_matrix) = data_signal_matrix[,1]
	data_signal_matrix = as.matrix(data_signal_matrix[,c(-1,-2,-3)])
	return(data_signal_matrix)
}

### read DNA region matrix
read_DNAregion_matrix = function(inputfile){
	### read matrix as string matrix
	data_signal_matrix = as.matrix(read.table(inputfile, header=F, sep='\t'))
	rownames(data_signal_matrix) = data_signal_matrix[,1]
	data_signal_matrix = data_signal_matrix[,-1]
	return(data_signal_matrix)
}

### get DNA region intervals
index_region_matrix = read_DNAregion_matrix('DNA_regin_210k_indexsort_onlyinterval.txt')#[c(1:100000),]
### index region length
index_region_length = as.numeric(index_region_matrix[,3]) - as.numeric(index_region_matrix[,2])#[c(1:100000)]
total_len = sum(index_region_length)

### get index binary label matrix
index_label_matrix = read_index_matrix('celltype.index.sorted.txt')#[c(1:100000)]

### get index signal matrix
index_sig_matrix = read_signal_matrix('celltype.index.signal.sorted.txt')#[c(1:100000),]
NOcREs_sig_matrix = read_wg_signal_matrix('celltype.tpm.NOcRE.txt')#[c(1:100000),]
wg_sig_matrix_log2 = log2(read_wg_signal_matrix('celltype.tpm.sorted.txt')+0.01)#[c(1:100000),]


### get index_set binary labels
index_set_filtered_matrix = read_signal_matrix('celltype.index_set_filtered.sorted.txt')#[c(1,2,3),]
### to numeric matrix
class(index_set_filtered_matrix) = "numeric"
index_set_filtered_label = rownames(index_set_filtered_matrix)

### get multiple variable normal distribution from previous peak calling result
### get new index calling probability

### Previous matrix (plus 0.01 then log2 transform)
index_sig_matrix_log2 = log2(index_sig_matrix+0.01)

### initialize index call p matrix
index_sig_matrix_p_call = c()
enriched_index_set_position = rep(0, dim(index_sig_matrix_log2)[1])

### enriched index set:
for (i in c(1: length(index_set_filtered_label))){
	###
	print(paste('index set', toString(i), index_set_filtered_label[i], sep=':'))
	### get index set index position
	index_sig_i_position = index_label_matrix==index_set_filtered_label[i]
	### get each index's signal (plus 0.01 then log2 transform)
	index_sig_i_log2 = log2(index_sig_matrix[index_sig_i_position, ] + 0.01)
	### get index set peak proportion
	index_region_length_i = -log(sum(index_region_length[index_sig_i_position]/total_len))

	### save enriched_index_set DNA region id positions
	enriched_index_set_position = enriched_index_set_position + index_sig_i_position

	### boxplot index_set cell type signals
	png(paste('index_set_boxplot/index_set_boxplot.', toString(i), '.', index_set_filtered_label[i], '.png', sep=''))
	boxplot(index_sig_i_log2, ylim=c(min(index_sig_matrix_log2), quantile(index_sig_matrix_log2, probs=0.99)))
	dev.off()

	### get index_set mean vector & variance-covariance matrix 
	x_mean = colMeans(index_sig_i_log2)
	x_cov = cov(index_sig_i_log2)

	### read mvnorm function for matrix
	dmvnorm_all_slow = function(vector){
		library(mvtnorm)
		d = dmvnorm(vector, x_mean, x_cov)
		return(d)
	}

	### dnorm function (fast version)
	dnorm_fast = function(matrix, x_mean, x_cov){
		#data_r = t(apply(matrix, 1, FUN=function(x) x-x_mean))#t(matrix) - x_mean
		data_r = t(t(matrix) - x_mean)
		p = length(x_mean)
		### 
		x2 = apply(data_r  %*% solve(x_cov) * (data_r), 1, sum) 
		#print(dim(x2))
		d = det(x_cov) * 2 * pi 
		lp = (log(d) * p + x2)/2
		return(lp)
	}


	### calculate the likelihood density 
	print('start dmvnorm')
	index_sig_matrix_i_p = dnorm_fast(wg_sig_matrix_log2, x_mean, x_cov)
	#print(length(index_sig_matrix_i_p))
	### consider overall index set region length
	index_sig_matrix_i_p_p = index_sig_matrix_i_p + index_region_length_i

	### append index call matrix (-log(p))
	index_sig_matrix_p_call = cbind(index_sig_matrix_p_call, index_sig_matrix_i_p_p)
}



### NOT enriched index set:
for (i in c(1)){
	###
	print(paste('index set', 'NOT enriched', sep=':'))
	### get index set index position
	index_sig_i_position = enriched_index_set_position==0
	### get each index's signal (plus 0.01 then log2 transform)
	index_sig_i_log2 = log2(index_sig_matrix[index_sig_i_position, ] + 0.01)
	### get index set peak proportion
	index_region_length_i = -log(sum(index_region_length[index_sig_i_position]/total_len))

	### save enriched_index_set DNA region id positions
	enriched_index_set_position = enriched_index_set_position + index_sig_i_position

	### boxplot index_set cell type signals
	png(paste('index_set_boxplot/index_set_boxplot.', toString(i), '.', index_set_filtered_label[i], '.png', sep=''))
	boxplot(index_sig_i_log2, ylim=c(min(index_sig_matrix_log2), quantile(index_sig_matrix_log2, probs=0.99)))
	dev.off()

	### get index_set mean vector & variance-covariance matrix 
	x_mean = colMeans(index_sig_i_log2)
	x_cov = cov(index_sig_i_log2)

	### read mvnorm function for matrix
	dmvnorm_all_slow = function(vector){
		library(mvtnorm)
		d = dmvnorm(vector, x_mean, x_cov)
		return(d)
	}

	### dnorm function (fast version)
	dnorm_fast = function(matrix, x_mean, x_cov){
		#data_r = t(apply(matrix, 1, FUN=function(x) x-x_mean))#t(matrix) - x_mean
		data_r = t(t(matrix) - x_mean)
		p = length(x_mean)
		### 
		x2 = apply(data_r  %*% solve(x_cov) * (data_r), 1, sum) 
		#print(dim(x2))
		d = det(x_cov) * 2 * pi 
		lp = (log(d) * p + x2)/2
		return(lp)
	}


	### calculate the likelihood density 
	print('start dmvnorm')
	index_sig_matrix_i_p = dnorm_fast(wg_sig_matrix_log2, x_mean, x_cov)
	#print(length(index_sig_matrix_i_p))
	### consider overall index set region length
	index_sig_matrix_i_p_p = index_sig_matrix_i_p + index_region_length_i

	### append index call matrix (-log(p))
	index_sig_matrix_p_call = cbind(index_sig_matrix_p_call, index_sig_matrix_i_p_p)
}



### NOT enriched index set:
for (i in c(1)){
	###
	print(paste('index set', 'NOT enriched whole genome', sep=':'))
	### get each index's signal (plus 0.01 then log2 transform)
	index_sig_i_log2 = log2(NOcREs_sig_matrix + 0.01)
	### get index set peak proportion
	index_region_length_i = -log(sum(index_region_length[index_sig_i_position]/total_len))

	### save enriched_index_set DNA region id positions
	enriched_index_set_position = enriched_index_set_position + index_sig_i_position

	### boxplot index_set cell type signals
	png(paste('index_set_boxplot/index_set_boxplot.', toString(i), '.', index_set_filtered_label[i], '.png', sep=''))
	boxplot(index_sig_i_log2, ylim=c(min(index_sig_matrix_log2), quantile(index_sig_matrix_log2, probs=0.99)))
	dev.off()

	### get index_set mean vector & variance-covariance matrix 
	x_mean = colMeans(index_sig_i_log2)
	x_cov = cov(index_sig_i_log2)

	### read mvnorm function for matrix
	dmvnorm_all_slow = function(vector){
		library(mvtnorm)
		d = dmvnorm(vector, x_mean, x_cov)
		return(d)
	}

	### dnorm function (fast version)
	dnorm_fast = function(matrix, x_mean, x_cov){
		#data_r = t(apply(matrix, 1, FUN=function(x) x-x_mean))#t(matrix) - x_mean
		data_r = t(t(matrix) - x_mean)
		p = length(x_mean)
		### 
		x2 = apply(data_r  %*% solve(x_cov) * (data_r), 1, sum) 
		#print(dim(x2))
		d = det(x_cov) * 2 * pi 
		lp = (log(d) * p + x2)/2
		return(lp)
	}


	### calculate the likelihood density 
	print('start dmvnorm')
	index_sig_matrix_i_p = dnorm_fast(wg_sig_matrix_log2, x_mean, x_cov)
	#print(length(index_sig_matrix_i_p))
	### consider overall index set region length
	index_sig_matrix_i_p_p = index_sig_matrix_i_p + index_region_length_i

	### append index call matrix (-log(p))
	index_sig_matrix_p_call = cbind(index_sig_matrix_p_call, index_sig_matrix_i_p_p)
}


#print(dim(index_sig_matrix_p_call))
#print(head(index_sig_matrix_p_call))

colnames(index_sig_matrix_p_call) = c(1: dim(index_sig_matrix_p_call)[2])
#rownames(index_sig_matrix_p_call) = rownames(index_sig_matrix)

### the sum of probability is standardized to 1
#print('get index call probability matrix')
#index_sig_matrix_p_call = t( apply(index_sig_matrix_p_call, MARGIN=1, FUN=function(x) x/sum(x)) )

### write all index set call & probability 
write.table(index_sig_matrix_p_call, 'index_sig_matrix_index_set_p.txt', quote=F, sep='\t', row.names = TRUE, col.names = NA)


### get belonging index set id
print('get index call index_set_id & probability')
index_sig_matrix_index_set_call = t( apply((index_sig_matrix_p_call), MARGIN=1, FUN=function(x) c(which.max(x), max(x))) )
colnames(index_sig_matrix_index_set_call) = c('index_set_id', 'probability')
### write all index set probility matrix 
write.table(index_sig_matrix_index_set_call, 'index_sig_matrix_index_set_call.txt', quote=F, sep='\t', row.names = TRUE, col.names = )


############################################
print('sort index based on new index call')
### sort input signal matrix based on the new index set call
index_sig_matrix_sort = index_sig_matrix[order(index_sig_matrix_index_set_call[,1]), ]
### write all index set call & probability 
write.table(index_sig_matrix_sort, 'index_sig_matrix_index_set_call_sort.txt', quote=F, sep='\t', row.names = TRUE, col.names = )

############################################
print('get new index set numbers:')
### recalculate index set DNA region number
index_set_num = table(index_sig_matrix_index_set_call[,1])

### change previous index set number to new called number
index_set_filtered_matrix_recalnum = c()
for (i in c(1: dim(index_set_filtered_matrix)[1])){
	index_set_filtered_matrix_vec = index_set_filtered_matrix[i, ]
	### recalculate the number of DNA region in index set
	if (i %in% names(index_set_num)){
		index_set_filtered_matrix_vec_new = index_set_filtered_matrix_vec / max(index_set_filtered_matrix_vec) * index_set_num[names(index_set_num) == i]
		index_set_filtered_matrix_recalnum = cbind(index_set_filtered_matrix_recalnum, index_set_filtered_matrix_vec_new)
	} else{
		index_set_filtered_matrix_vec_new = index_set_filtered_matrix_vec - index_set_filtered_matrix_vec
		index_set_filtered_matrix_recalnum = cbind(index_set_filtered_matrix_recalnum, index_set_filtered_matrix_vec_new)		
	}
}
### transpose
index_set_filtered_matrix_recalnum = t( index_set_filtered_matrix_recalnum )

### add colnames & rownames
colnames(index_set_filtered_matrix_recalnum) = colnames(index_set_filtered_matrix)
rownames(index_set_filtered_matrix_recalnum) = rownames(index_set_filtered_matrix)

### write index_set_filtered_matrix_recalnum
write.table(index_set_filtered_matrix_recalnum, 'celltype.index_set_filtered.sorted.recalnum.txt', quote=F, sep='\t', row.names = TRUE, col.names = NA)

png('probability_hist.png')
hist(index_sig_matrix_index_set_call[,2], breaks = 50)#, xlim=c(0, 150))
dev.off()

############################################
print('get new index set numbers (threshold 0.95):')
### recalculate index set DNA region number
index_set_num = table(index_sig_matrix_index_set_call[index_sig_matrix_index_set_call[,2]>=0,1])

### change previous index set number to new called number
index_set_filtered_matrix_recalnum = c()
for (i in c(1: dim(index_set_filtered_matrix)[1])){
	index_set_filtered_matrix_vec = index_set_filtered_matrix[i, ]
	### recalculate the number of DNA region in index set
	if (i %in% names(index_set_num)){
		index_set_filtered_matrix_vec_new = index_set_filtered_matrix_vec / max(index_set_filtered_matrix_vec) * index_set_num[names(index_set_num) == i]
		index_set_filtered_matrix_recalnum = cbind(index_set_filtered_matrix_recalnum, index_set_filtered_matrix_vec_new)
	} else{
		index_set_filtered_matrix_vec_new = index_set_filtered_matrix_vec - index_set_filtered_matrix_vec
		index_set_filtered_matrix_recalnum = cbind(index_set_filtered_matrix_recalnum, index_set_filtered_matrix_vec_new)		
	}
}
### transpose
index_set_filtered_matrix_recalnum = t( index_set_filtered_matrix_recalnum )

### add colnames & rownames
colnames(index_set_filtered_matrix_recalnum) = colnames(index_set_filtered_matrix)
rownames(index_set_filtered_matrix_recalnum) = rownames(index_set_filtered_matrix)

### write index_set_filtered_matrix_recalnum
write.table(index_set_filtered_matrix_recalnum, 'celltype.index_set_filtered.sorted.recalnum.thresh.txt', quote=F, sep='\t', row.names = TRUE, col.names = NA)

