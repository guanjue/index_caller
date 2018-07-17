library(mixtools)

split_label = function(x, ct_num){
	xs = unlist(strsplit(toString(x), "_"))
	if (xs[1] == 'N'){
		xs =  rep('N', ct_num)
	} else if (xs[1] == 'X'){
		xs =  rep('X', ct_num)
	}
	return(xs)
}


data_od = read.table('atac_20cell.sig.18.txt.ic.txt', header = F)
ct_names = read.table('ct_list.18.txt')
ct_num = dim(ct_names)[1]

data_sig = data_od[,c(-1,-2)]
data_old_label = as.matrix(data_od[,2])
data_old_label_unique = unique(data_old_label)

data_old_labe_mat = t(apply(data_old_label, 1, function(x) split_label(x, ct_num) ))
data_new_labe_mat = c()

iteration_num = 100
gmm_k = 2

for (i in 1:dim(data_sig)[2]){
	print(paste('cell type ', toString(i), ': ', toString(ct_names[i,]), sep=''))
	### get each column (cell type) signal & label
	data_sig_ct = data_sig[,i]
	data_label_ct = data_old_labe_mat[,i]
	###### iteration update data_sig_ct_0 & data_sig_ct_1
	for (j in 1:iteration_num){
		### initializing dif_num for each column for each iteration
		dif_num = 0
		### get peak & background
		data_sig_ct_0 = data_sig_ct[data_label_ct=='0']
		data_sig_ct_1 = data_sig_ct[data_label_ct!='0']

		for (index_set in data_old_label_unique){
			### get each index-set signal vector in a column
			select_index_peak_id = (data_old_label == index_set)
			each_index_set_ctsig = data_sig_ct[select_index_peak_id]
			old_label = data_label_ct[select_index_peak_id][1]
			### two sample t-test: t_test_greater_than_bg
			t_test_greater_than_bg = t.test(each_index_set_ctsig, data_sig_ct_0, alternative = "greater", var.equal = FALSE)
			t_test_greater_than_bg_p = t_test_greater_than_bg$p.value
			if (t_test_greater_than_bg_p < (1e-3/ct_num)){
				### relabel as level 1 peak
				new_label = '1'
				if (old_label != new_label){
					data_label_ct[data_old_label == index_set] = new_label
				}
			} else {
				### relabel as bg
				new_label = '0'
				if (old_label != new_label){
					data_label_ct[data_old_label == index_set] = new_label
				}
			}
			#### count dif num of the index-set
			if (old_label != new_label){
				dif_num = dif_num + 1
			}		
		}
		print(paste('change label number:', toString(dif_num), sep = ' '))
		if (dif_num==0){
			print('converged!')
			break
		}
	}

	### last iteration: use gmm to cluster 1 & 2
	### fit gmm
	png(paste('ct_sig.', toString(i), '.density.0.png', sep=''))
	#plot(density(data_sig_ct_1), lty=2, lwd=2)
	hist(data_sig_ct_0, breaks=100)
	dev.off()
	png(paste('ct_sig.', toString(i), '.density.1.png', sep=''))
	#plot(density(data_sig_ct_1), lty=2, lwd=2)
	hist(data_sig_ct_1, breaks=100)
	dev.off()
	mixmdl = normalmixEM(data_sig_ct_1,k = gmm_k)
	png(paste('ct_sig.', toString(i), '.gmm.png', sep=''))
	plot(mixmdl,which=2)
	lines(density(data_sig_ct_1), lty=2, lwd=2)
	dev.off()
	###### get post mean
	post = apply(mixmdl$posterior, 1, function(x) which.max(x)-1)
	c1_a = signal_vec[post == 0]
	c2_a = signal_vec[post == 1]
	c1_m = mean(c1_a)
	c2_m = mean(c2_a)
	if (min(c(c1_m, c2_m))==c1_m){
		c1 = c1_a
		c2 = c2_a
	} else {
		c1 = c2_a
		c2 = c2_a
	}
	level1_lim = min(c1)
	level2_lim = min(c2)
	print('level1_lim: ')
	print(level1_lim)
	print('level2_lim: ')
	print(level2_lim)
	### relabel 1 as 1 or 2
	data_label_ct[data_sig_ct_1 >= level1_lim] = '1'
	data_label_ct[data_sig_ct_1 >= level2_lim] = '2'
	### add new label vector
	data_new_labe_mat = cbind(data_new_labe_mat, data_label_ct)

}

data_new_labe_vec = apply(data_new_labe_mat, 1, function(x) paste(x, collapse="_"))

### add id for each new index label to avoid same new label
l = 0
for (index_set in data_old_label_unique){
	l = l+1
	check_if_N_X = split_label(index_set, ct_num)[1]
	if (check_if_N_X == 'X'){
		new_label_vec = data_new_labe_vec[data_old_label == index_set]
		new_label_vec_1 = paste('X', new_label_vec[1], toString(l), sep=":")
		data_new_labe_vec[data_old_label == index_set] = new_label_vec_1		
	#} else if (check_if_N_X == 'N'){
	#	new_label_vec = (data_new_labe_vec[data_old_label == index_set])
	#	new_label_vec_1 = paste('N', new_label_vec[1], toString(l), sep=":")
	#	data_new_labe_vec[data_old_label == index_set] = new_label_vec_1			
	} else {
		new_label_vec = data_new_labe_vec[data_old_label == index_set]
		new_label_vec_1 = paste(new_label_vec[1], toString(l), sep=":")
		data_new_labe_vec[data_old_label == index_set] = new_label_vec_1	
	}
}


data_new = cbind(data_od[,1], data_new_labe_vec, data_od[,3:(dim(data_od)[2])])

write.table(data_new, 'atac_20cell.sig.18.txt.ic.Nlabel.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')


data_new_labe_vec_uniq =  unique(data_new_labe_vec)
new_index_sig_mean_mat = c()

for (new_index in data_new_labe_vec_uniq){
	new_index_sig = as.matrix(data_new[data_new_labe_vec==new_index,c(-1,-2)])
	#print(head(new_index_sig))
	new_index_sig_mean = colMeans(new_index_sig)
	new_index_sig_mean_mat = rbind(new_index_sig_mean_mat, new_index_sig_mean)
}

new_index_sig_mean_mat = cbind(data_new_labe_vec_uniq, new_index_sig_mean_mat)
write.table(new_index_sig_mean_mat, 'atac_20cell.sig.18.txt.ic.Nlabel.meansig.txt', quote=FALSE, col.names=FALSE, row.names=FALSE, sep='\t')



####################################################
### use pheatmap to plot heatmaps
color_heatmap = function(color_matrix, high_color, low_color, format, outputname){
	library(pheatmap)
	format(outputname, width = 200, height = 200) ### output name
	par() ### set heatmap margins
	### plot pheatmap
	my_colorbar=colorRampPalette(c(low_color, high_color))(n = 128)
	col_breaks = c(seq(0, 2000,length=33))
	pheatmap(color_matrix, color=my_colorbar, cluster_cols = FALSE,cluster_rows=FALSE,annotation_names_row=FALSE,annotation_names_col=TRUE,show_rownames=TRUE,show_colnames=TRUE)
	dev.off()
}

new_index_sig_mean_mat = as.matrix(new_index_sig_mean_mat[order(data_new_labe_vec_uniq),-1])
class(new_index_sig_mean_mat) = 'numeric'
rownames(new_index_sig_mean_mat) = new_index_sig_mean_mat[,1]

color_heatmap(new_index_sig_mean_mat, 'blue', 'white', pdf, 'atac_20cell.sig.18.txt.ic.Nlabel.meansig.pdf')




















