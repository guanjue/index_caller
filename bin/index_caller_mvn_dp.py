import os
import numpy as np

################################################################################################
### read 2d array
def read2d_array(filename,dtype_used):
	data=open(filename,'r')
	data0=[]
	for records in data:
		tmp = [x.strip() for x in records.split('\t')]
		data0.append(tmp)
	data0 = np.array(data0,dtype=dtype_used)
	data.close()
	return data0

################################################################################################
### write 2d matrix
def write2d_array(array,output):
	r1=open(output,'w')
	for records in array:
		for i in range(0,len(records)-1):
			r1.write(str(records[i])+'\t')
		r1.write(str(records[len(records)-1])+'\n')
	r1.close()

################################################################################################
### cov2corr
def cov2corr(cov, return_std=False):
	cov = np.asanyarray(cov)
	std_ = np.sqrt(np.diag(cov))
	#print(std_)
	std_[std_<=0.1]=0.1
	#print(std_)
	corr = cov / np.outer(std_, std_)
	if return_std:
		return corr, std_
	else:
		return corr

################################################################################################
### corr2cov
def corr2cov(corr, std):
	corr = np.asanyarray(corr)
	std_ = np.asanyarray(std)
	cov = corr * np.outer(std_, std_)
	return cov

################################################################################################
### get convert bedtools window output to matrix of pk and intersect function label info
#data_new = read2d_array('atac_20cell_wg.indexcaller.txt', str)
data_new = read2d_array('atac_20cell.sig.5.txt', str)



##################
###### randomly select # bins from whole genome as index-set 0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0
np.random.seed(2018)
'''
random_sample_num = 100000
data_new_bg_id = np.random.randint(data_new.shape[0], size=random_sample_num)
data_new_bg = data_new[data_new_bg_id,:]
data_new_bg_label = np.repeat('0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0', random_sample_num).reshape(random_sample_num, 1)
data_new_bg = np.concatenate((data_new_bg[:,0].reshape(random_sample_num,1), data_new_bg_label, data_new_bg[:,1:]), axis=1)

data_new_sig = data_new[:,1:].astype(str)
data_New_sig_matrix_bed = data_new[:,0]
'''

data = read2d_array('atac_20cell.sig.5.txt', str)


##################
###### add randomly selected signals to the train index matrix based on peak calling results
#data = np.concatenate((data, data_new_bg), axis=0)

##################
###### set std upper limit for each index-set xi
std_upper_lim = np.sqrt(3.0)
#std_lower_lim = np.sqrt(1.0)


data_sig_mean_dict = {}
data_sig_cov_dict = {}
data_sig_p_dict = {}
data_sig_cov_all = []

change_num_all = []
for iteration_num in range(0,100):
	print('iteration_num: ' + str(iteration_num) + '!!!')
	#print(data[0:3])
	data_sig_dict = {}
	data_index_vec = []
	data_index_count = {}

	for info in data:
		index = info[1]
		if not (index in data_sig_dict):
			data_sig_dict[index] = [info[2:]]
			data_index_vec.append(index)
			data_index_count[index] = 1
		else:
			data_sig_dict[index].append(info[2:])
			data_index_count[index] = data_index_count[index] + 1

	data_sig_matrix = data[:,2:].astype(float)
	data_sig_matrix_bed = data[:,0:2].astype(str)

	###### get mean vector & covariance matrix & prior vector of index-sets based peak calling index-sets
	data_sig_cov_all = []
	total_row_num = data.shape[0]
	ct_num = data.shape[1]-2

	if iteration_num == 0:
		for index in data_sig_dict:
			data_sig_matrix_i = np.array(data_sig_dict[index], dtype=float)
			data_sig_dict[index] = data_sig_matrix_i
			data_sig_p_dict[index] = float(data_sig_matrix_i.shape[0]) / total_row_num
			data_sig_mean_dict[index] = np.mean(data_sig_matrix_i, axis=0)
			cov_matrix = np.cov(data_sig_matrix_i.T)

			
			##################
			###### get min std of index-set i
			std_vector = np.sqrt(np.diag(cov_matrix))
			std_0_num = np.sum(std_vector==0)
			
			if std_0_num > 0:
				for i in np.argwhere(std_vector==0):
					mu, sigma = 0, 0.1
					peak_num = data_sig_matrix_i.shape[0]
					add_noise = np.random.normal(mu, sigma, peak_num)
					data_sig_matrix_i[:,i] = add_noise.reshape(peak_num,1)
				cov_matrix = np.cov(data_sig_matrix_i.T)
			
			##################
			###### get max std of index-set i
			std_max = np.max(std_vector)

			if std_max > std_upper_lim:
				corr_i, std_i = cov2corr(cov_matrix, return_std=True)
				std_adj_num = np.sum(std_i>std_upper_lim)
				std_upper_lim_add_noise = np.random.normal(std_upper_lim, 0.01, std_adj_num)
				std_i[std_i>std_upper_lim] = 3.0#std_upper_lim_add_noise
				#std_i = std_i / np.max(std_i) * std_upper_lim
				cov_matrix = corr2cov(corr_i, std_i)
				#print(data_sig_cov_i)
			#data_sig_cov_i[np.isnan(data_sig_cov_i)] = 1e-5
			data_sig_cov_dict[index] = cov_matrix
			data_sig_cov_all.append(data_sig_cov_dict[index])
		### for exmpty 0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0
		data_sig_mean_dict['0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0'] = np.repeat(2.0, ct_num)
		cov_matrix = np.identity(ct_num)
		data_sig_cov_dict[index] = cov_matrix
		data_sig_cov_all.append(data_sig_cov_dict[index])		

	
	else:
		###### if the sample number decrease to 0
		for index in data_sig_mean_dict:
			if index in data_sig_dict:
				data_sig_matrix_i = np.array(data_sig_dict[index], dtype=float)
				data_sig_dict[index] = data_sig_matrix_i
				data_sig_p_dict[index] = float(data_sig_matrix_i.shape[0]) / total_row_num
				data_sig_mean_dict[index] = np.mean(data_sig_matrix_i, axis=0)
				cov_matrix = np.cov(data_sig_matrix_i.T)
				##################
				###### get min std of index-set i
				std_vector = np.sqrt(np.diag(cov_matrix))
				std_0_num = np.sum(std_vector==0)
				
				if std_0_num > 0:
					for i in np.argwhere(std_vector==0):
						mu, sigma = 0, 0.1
						peak_num = data_sig_matrix_i.shape[0]
						add_noise = np.random.normal(mu, sigma, peak_num)
						data_sig_matrix_i[:,i] = add_noise.reshape(peak_num,1)
					cov_matrix = np.cov(data_sig_matrix_i.T)
				
				##################
				###### get max std of index-set i
				std_max = np.max(std_vector)	
				if std_max > std_upper_lim:
					corr_i, std_i = cov2corr(cov_matrix, return_std=True)
					std_adj_num = np.sum(std_i>std_upper_lim)
					std_upper_lim_add_noise = np.random.normal(std_upper_lim, 0.01, std_adj_num)
					std_i[std_i>std_upper_lim] = 3.0#std_upper_lim_add_noise
					#std_i = std_i / np.max(std_i) * std_upper_lim
					cov_matrix = corr2cov(corr_i, std_i)
					#print(data_sig_cov_i)
			else:
				#print(index)
				data_sig_mean_dict[index] = np.repeat(2.0, ct_num)
				cov_matrix = np.identity(ct_num)

			#data_sig_cov_i[np.isnan(data_sig_cov_i)] = 1e-5
			data_sig_cov_dict[index] = cov_matrix
			data_sig_cov_all.append(data_sig_cov_dict[index])		

	##################
	###### add non-peak prior to prior vector
	#data_sig_p_dict['0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0'] = 1.0 - 215120.0/13554672.0

	data_sig_cov_all = np.array(data_sig_cov_all)
	mvn_qda_matrix = []
	data_index_count_vec = []
	data_index_vec = np.sort(np.array(data_index_vec))




	all_cov = []
	#constant = -18/2*np.log(2*3.14)



	##################
	###### loop index-set
	for index in data_index_vec:
		data_sig_mean_i = data_sig_mean_dict[index]
		data_sig_cov_i = data_sig_cov_dict[index]
		data_sig_p_i = data_sig_p_dict[index]
		data_index_count_vec.append(data_index_count[index])
		all_cov.append(data_sig_cov_i)
		#print(index)
		#print(data_sig_cov_dict[index])
		small_num = 0.0
		data_sig_cov_i[np.isnan(data_sig_cov_i)] = 0.0
		

		try:
			np.linalg.det(data_sig_cov_i+small_num)
		except np.linalg.LinAlgError as err:
			if 'Singular matrix' in str(err):
				print(index)
				print('check problem np.linalg.det')
				print(data_sig_cov_i)
				print(np.linalg.det(data_sig_cov_i+small_num))
				print(np.log(np.linalg.det(data_sig_cov_i+small_num)))
				print('check problem ok')

		#if iteration_num == 8:
		#print(data_sig_cov_i)
		if (np.isnan(np.log(np.linalg.det(data_sig_cov_i+small_num)+small_num))):
			print(index)
			print('check problem np.linalg.det')
			print(data_sig_cov_i)
			print(np.linalg.det(data_sig_cov_i+small_num))
			print(np.log(np.linalg.det(data_sig_cov_i+small_num)))
			print('check problem ok')

		#if (np.isnan(np.sum(np.linalg.inv(data_sig_cov_i+small_num)))):
		try:
			np.linalg.inv(data_sig_cov_i+small_num)
		except np.linalg.LinAlgError as err:
			if 'Singular matrix' in str(err):
				print(index)
				print('check problem np.linalg.inv')
				print(data_sig_cov_i)
				print(np.linalg.inv(data_sig_cov_i+small_num))
				print(np.log(np.linalg.inv(data_sig_cov_i+small_num)))
				print('check problem ok')

		'''
		if np.isnan(np.sum(mvn_qda)):
			print('cov')
			print(data_sig_cov_i)
			print('corr_i')
			print(corr_i)
			#print(data_sig_cov_dict[index])
			print('sig')
			print(data_sig_dict[index])
			print('lim')
			print(np.max(data_sig_dict[index], axis=0))
			print(np.min(data_sig_dict[index], axis=0))
			print('dif')
			print(np.max(data_sig_dict[index], axis=0)-np.min(data_sig_dict[index], axis=0))
			print('mvn_qda')
			print(- np.sum(np.dot((data_sig_matrix-data_sig_mean_i), np.linalg.inv(data_sig_cov_i)) * (data_sig_matrix-data_sig_mean_i), axis=1))
			print(-1/2*np.log(np.linalg.det(data_sig_cov_i)))
		#print(index)
		#print(mvn_qda)
		'''

		mvn_qda = -1/2*np.log(np.linalg.det(data_sig_cov_i+small_num)+small_num) - np.sum(np.dot((data_sig_matrix-data_sig_mean_i), np.linalg.inv(data_sig_cov_i+small_num)) * (data_sig_matrix-data_sig_mean_i), axis=1) + np.log(data_sig_p_i) #+ constant
		mvn_qda_matrix.append(mvn_qda)

	#break
	all_cov = (np.array(all_cov))
	#print(all_cov.shape)
	all_cov = all_cov.reshape(all_cov.shape[0]*all_cov.shape[1]*all_cov.shape[2],1)
	write2d_array(all_cov, 'all_cov.txt')



	##################
	###### convert mvn density score to p-value
	mvn_qda_matrix = np.transpose(np.array(mvn_qda_matrix))
	mvn_qda_p_matrix = np.exp(mvn_qda_matrix)#-mvn_qda_matrix.max(axis=1, keepdims=True))
	
	mvn_qda_p_matrix = mvn_qda_p_matrix / mvn_qda_p_matrix.sum(axis=1, keepdims=True)

	'''
	print(mvn_qda_matrix[0:10,:])
	print(np.exp((mvn_qda_matrix[0:10,:])))
	print(mvn_qda_p_matrix[0:10,:])
	print(mvn_qda_p_matrix.shape)
	print(np.sum(mvn_qda_p_matrix, axis=1).shape)
	print(np.sum(mvn_qda_p_matrix, axis=1)[0:10])
	'''
	##################
	###### if 1-sum_mvn_p > 0.5 -> set to index-set 0
	#if iteration_num == 0:
	#	mvn_qda_p0_matrix = 1.0-np.max(mvn_qda_p_matrix, axis=1)
	#	mvn_qda_p_matrix = np.concatenate((mvn_qda_p_matrix, mvn_qda_p0_matrix.reshape(mvn_qda_p0_matrix.shape[0], 1)), axis=1)

	'''
	print('no pk:')
	print(np.sum(mvn_qda_p0_matrix>0.5))
	'''
	



	##################
	###### set the index-set with highest p as the index-set label for each bin
	#print(mvn_qda_p_matrix[0,:])
	cluster_id = np.argmax(mvn_qda_p_matrix, axis=1)
	#print(cluster_id)
	#print(mvn_qda_p_matrix[0:3,:])
	#print(cluster_id[0:10])



	##################
	###### get the highest probability
	y_pred = np.amax(mvn_qda_p_matrix, axis=1)
	#print(cluster_id.shape)
	#data_index_vec = np.append(data_index_vec, '0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0')

	### count index-set number
	data_index_pred_count = {}
	index_pred_vec = []
	for i in cluster_id:
		index_pred = data_index_vec[i]
		index_pred_vec.append(index_pred)
		if not (index_pred in data_index_pred_count):
			data_index_pred_count[index_pred] = 1
		else:
			data_index_pred_count[index_pred] = data_index_pred_count[index_pred] +1

	### replace rare index-sets
	for i in range(0,len(index_pred_vec)):
		index_pred = index_pred_vec[i]
		if data_index_pred_count[index_pred] < 100:
			index_pred_vec[i] = 'X_X_X_X_X'

	### recount index-set number
	data_index_pred_count = {}
	for i in range(0,len(index_pred_vec)):
		index_pred = index_pred_vec[i]
		if not (index_pred in data_index_pred_count):
			data_index_pred_count[index_pred] = 1
		else:
			data_index_pred_count[index_pred] = data_index_pred_count[index_pred] +1



	count_table = open('iteration_' + str(iteration_num) + '_count_table.txt','w')
	change_num = 0
	for index in data_index_vec:
		print(index)
		count_table.write(index+'\t')

		if (index in data_index_pred_count) and (index in data_index_count):
			print(str(data_index_pred_count[index])+'-'+str(data_index_count[index])+'=' +str(data_index_pred_count[index]-data_index_count[index]))
			count_table.write(str(data_index_pred_count[index])+'-'+str(data_index_count[index])+'=' +str(data_index_pred_count[index]-data_index_count[index])+'\n')
			change_num = change_num + abs(data_index_pred_count[index] - data_index_count[index])
		elif (index in data_index_pred_count) and not (index in data_index_count):
			print(str(data_index_pred_count[index])+'-'+str(0)+'=' +str(data_index_pred_count[index]))
			count_table.write(str(data_index_pred_count[index])+'-'+str(0)+'=' +str(data_index_pred_count[index])+'\n')
			change_num = change_num + abs(data_index_pred_count[index] - 0)
		elif index in data_index_count:
			print(str(0)+'-'+str(data_index_count[index])+'=' +str(0-data_index_count[index]))
			count_table.write(str(0)+'-'+str(data_index_count[index])+'=' +str(0-data_index_count[index])+'\n')
			change_num = change_num + abs(0-data_index_count[index])
		else:
			print(str(0)+'-'+str(0)+'=' +str(0))
			count_table.write(str(0)+'-'+str(0)+'=' +str(0)+'\n')
	count_table.close()

	#print(mvn_qda_matrix.shape)
	#print(mvn_qda_p_matrix.shape)
	#print(data_sig_matrix.shape)
	change_num1 = 0
	for i in range(0, len(index_pred_vec)):
		if data[i,1] != index_pred_vec[i]:
			change_num1 = change_num1+1

	change_num_all.append([change_num, change_num1])
	print('change_num:')
	print(len(data_index_vec))
	print('iteration_num: ' + str(iteration_num) + '!!!' + str(change_num))
	print('iteration_num: ' + str(iteration_num) + '!!!' + str(change_num1))
	### updata previous signal matrix



	data[:,1] = index_pred_vec

change_num_all = np.array(change_num_all)
write2d_array(change_num_all, 'change_num_all.txt')





data_sig_pred_dict = {}
data_sig_pred_vec = []

for i, info in zip(cluster_id, data):
	index = data_index_vec[i]
	data_sig_pred_vec.append(index)
	if not (index in data_sig_pred_dict):
		data_sig_pred_dict[index] = [ info[2:] ]
	else:
		data_sig_pred_dict[index].append(info[2:])



data[:,1] = data_sig_pred_vec

data = np.array(data)

write2d_array(data, 'data.adj.txt')



data_sig_pred_dict = {}
cluster_vec = []
for records in data:
	cluster = records[1]
	if not (cluster in data_sig_pred_dict):
		data_sig_pred_dict[cluster] = [ records[2:] ]
		cluster_vec.append(cluster)
	else:
		data_sig_pred_dict[cluster].append(records[2:])

cluster_vec = np.sort(np.array(cluster_vec))

pred_cluster_mean = []

for cluster in cluster_vec:
	sig_mat = np.array(data_sig_pred_dict[cluster], dtype=float)
	sig_mat_mean = np.mean(sig_mat, axis=0)
	pred_cluster_mean.append(sig_mat_mean)

pred_cluster_mean = np.array(pred_cluster_mean)

data_sig_out = np.concatenate((cluster_vec.reshape((cluster_vec).shape[0],1), pred_cluster_mean), axis=1)
write2d_array(data_sig_out, 'qda_atac_wg.mvn_check.5.txt')



#print(data_New_sig_matrix_bed.shape)
#print(y_pred.shape)
#data_sig_out_matrix = np.concatenate((data_New_sig_matrix_bed.reshape((y_pred).shape[0],1), np.array(index_pred_vec).reshape((y_pred).shape[0],1), np.array(y_pred).reshape((y_pred).shape[0],1), data_new_sig), axis=1)
#write2d_array(data_sig_out_matrix, 'qda_atac_matrix_wg_check.mvn.5.txt')





