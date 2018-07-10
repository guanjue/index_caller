import os
import numpy as np
from scipy import stats

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
def cov2corr(cov, std_lower_lim, return_std=False):
	cov = np.asanyarray(cov)
	std_ = np.sqrt(np.diag(cov))
	#print(std_)
	###### set lower limit for std vector
	std_[std_<=std_lower_lim]=std_lower_lim
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
def index_caller(input_signal_mat, output_signal_mat, index_count_lim, std_upper_lim, std_lower_lim, iteration_num_all, alpha_for_empty_index_set, random_seed):
	### get convert bedtools window output to matrix of pk and intersect function label info
	#data_new = read2d_array('atac_20cell_wg.indexcaller.txt', str)
	#data_new = read2d_array('atac_20cell.sig.10.txt', str)

	#input_signal_mat = 'atac_20cell.sig.10.txt'
	#output_signal_mat = 'data.adj.10.txt'
	#index_count_lim = 100
	#std_upper_lim = 3.0
	#std_lower_lim = 0.1
	#iteration_num_all = 100
	#alpha_for_empty_index_set = 1.0
	#random_seed = 2018

	##################
	###### randomly select # bins from whole genome as index-set 0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0
	np.random.seed(random_seed)
	'''
	random_sample_num = 100000
	data_new_bg_id = np.random.randint(data_new.shape[0], size=random_sample_num)
	data_new_bg = data_new[data_new_bg_id,:]
	data_new_bg_label = np.repeat('0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0', random_sample_num).reshape(random_sample_num, 1)
	data_new_bg = np.concatenate((data_new_bg[:,0].reshape(random_sample_num,1), data_new_bg_label, data_new_bg[:,1:]), axis=1)

	data_new_sig = data_new[:,1:].astype(str)
	data_New_sig_matrix_bed = data_new[:,0]
	'''

	###### read input signal matrix
	data = read2d_array(input_signal_mat, str)

	#data_sig = np.array(data[:,2:], dtype=float)
	#data_sig = np.add(data_sig, np.random.uniform(l_add, h_add, data.shape[0] * (data.shape[1]-2)).reshape(data.shape[0], (data.shape[1]-2)))
	#data[:,2:] = data_sig
	##################
	###### add randomly selected signals to the train index matrix based on peak calling results
	#data = np.concatenate((data, data_new_bg), axis=0)

	##################
	###### set std upper limit for each index-set xi
	std_upper_lim = np.sqrt(std_upper_lim)

	###### get mean vector & covariance matrix & prior probability vector
	data_sig_mean_dict = {}
	data_sig_cov_dict = {}
	data_sig_p_dict = {}
	#data_sig_cov_all = []

	###### get total mean vector and used it as the mean vector for new index-sets
	total_mean_vec = np.mean(np.array(data[:,2:], dtype=float), axis = 0)
	print(total_mean_vec)

	###### train index-caller
	change_num_all = []
	for iteration_num in range(0,iteration_num_all):
		print('iteration_num: ' + str(iteration_num) + '!!!')
		#print(data[0:3])
		data_sig_dict = {}
		data_index_vec = []
		data_index_count = {}

		###### get mean vector & covariance matrix & index-set peak counts (for prior calculation)
		for info in data:
			index = info[1]
			if not (index in data_sig_dict):
				data_sig_dict[index] = [info[2:]]
				data_index_vec.append(index)
				data_index_count[index] = 1
			else:
				data_sig_dict[index].append(info[2:])
				data_index_count[index] = data_index_count[index] + 1

		###### get signal matrix
		data_sig_matrix = data[:,2:].astype(float)
		###### get bed info & prior index-set label
		data_sig_matrix_bed = data[:,0:2].astype(str)

		###### get mean vector & covariance matrix & prior vector of index-sets based peak calling index-sets
		#data_sig_cov_all = []
		###### get total peak number for prior calculation
		total_row_num = data.shape[0]
		###### get colnumber for cell types
		ct_num = data.shape[1]-2
		###### add 1 for total row number for additional index-set
		if iteration_num == 0:
			###### for the first iteration
			for index in data_sig_dict:
				data_sig_matrix_i = np.array(data_sig_dict[index], dtype=float)
				data_sig_dict[index] = data_sig_matrix_i
				data_sig_p_dict[index] = float(data_sig_matrix_i.shape[0]) / (total_row_num + alpha_for_empty_index_set)
				data_sig_mean_dict[index] = np.mean(data_sig_matrix_i, axis=0)
				cov_matrix = np.cov(data_sig_matrix_i.T)
				###### add identity matrix to avoid matrix sigularity (!!! change to 1 to 0.01 to avoid singularity!!!)
				cov_matrix = np.add(cov_matrix, np.identity(cov_matrix.shape[0])/100)
				##################
				###### get max std of index-set i
				std_vector = np.sqrt(np.diag(cov_matrix))
				std_max = np.max(std_vector)
				###### if max std > std_upper_lim; then set std>std_upper_lim equal to std_upper_lim
				if std_max > std_upper_lim:
					###### first split coviance matrix to correlation matrix * std vector
					corr_i, std_i = cov2corr(cov_matrix, std_lower_lim, return_std=True)
					###### limit std
					std_i[std_i>std_upper_lim] = std_upper_lim
					###### merge correlation matrix & std vector back to covariance matrix
					cov_matrix = corr2cov(corr_i, std_i)
				###### add covariance matrix to dict
				data_sig_cov_dict[index] = cov_matrix
				#data_sig_cov_all.append(data_sig_cov_dict[index])
			### for exmpty 0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0_0
			#data_sig_mean_dict['0_0_0_0_0'] = np.repeat(2.0, ct_num)
			#cov_matrix = np.identity(ct_num)
			#data_sig_cov_dict[index] = cov_matrix
			#data_sig_cov_all.append(data_sig_cov_dict[index])				
		else:
			###### for the iterations other than 1
			###### if the sample number decrease to 0
			for index in data_sig_mean_dict:
				if index in data_sig_dict:
					###### if index-set already in the dict
					data_sig_matrix_i = np.array(data_sig_dict[index], dtype=float)
					data_sig_dict[index] = data_sig_matrix_i
					data_sig_p_dict[index] = float(data_sig_matrix_i.shape[0]) / (total_row_num + alpha_for_empty_index_set)
					data_sig_mean_dict[index] = np.mean(data_sig_matrix_i, axis=0)
					cov_matrix = np.cov(data_sig_matrix_i.T)
					###### add identity matrix to avoid matrix sigularity (!!! change to 1 to 0.01 to avoid singularity!!!)
					cov_matrix = np.add(cov_matrix, np.identity(cov_matrix.shape[0])/100)
					##################
					###### get max std of index-set i
					std_vector = np.sqrt(np.diag(cov_matrix))
					std_max = np.max(std_vector)
					###### if max std > std_upper_lim; then set std>std_upper_lim equal to std_upper_lim
					if std_max > std_upper_lim:
						###### first split coviance matrix to correlation matrix * std vector
						corr_i, std_i = cov2corr(cov_matrix, std_lower_lim, return_std=True)
						###### limit std
						std_i[std_i>std_upper_lim] = std_upper_lim
						###### merge correlation matrix & std vector back to covariance matrix
						cov_matrix = corr2cov(corr_i, std_i)
				else:
					###### if index-set NOT in the dict 
					###### (add a new index-set at each iteration for add more index-sets)
					###### use total means as mean vector
					data_sig_mean_dict[index] = np.repeat(2.0, ct_num)
					###### use identity matrix as covariance matrix
					cov_matrix = np.identity(ct_num)
				###### update covariance matrix in dict
				data_sig_cov_dict[index] = cov_matrix

		##################
		###### add non-peak prior for add new index-sets (mvn cluster)
		### add inactive cluster
		data_sig_p_dict['_'.join(['N'] * ct_num)+':'+str(iteration_num)] = 1.0 / (total_row_num+alpha_for_empty_index_set)
		data_sig_mean_dict['_'.join(['N'] * ct_num)+':'+str(iteration_num)] = total_mean_vec
		data_sig_cov_dict['_'.join(['N'] * ct_num)+':'+str(iteration_num)] = np.identity(ct_num)
		data_index_vec.append('_'.join(['N'] * ct_num)+':'+str(iteration_num))
		data_index_count['_'.join(['N'] * ct_num)+':'+str(iteration_num)] = 0.0

		###### qda score vector
		mvn_qda_matrix = []
		###### index-set peak counts vector
		data_index_count_vec = []

		###### sort index-set vector 
		data_index_vec = np.sort(np.array(data_index_vec))
		###### constant for qda score calculation
		#constant = -18/2*np.log(2*3.14)

		##################
		###### loop index-set
		for index in data_index_vec:
			###### get current iteration index-set mean vector & covariance matrix 
			###### & peak counts & prior probability 
			data_sig_mean_i = data_sig_mean_dict[index]
			data_sig_cov_i = data_sig_cov_dict[index]
			data_sig_p_i = data_sig_p_dict[index]
			data_index_count_vec.append(data_index_count[index])
			###### change na cov if there is na in covariance matrix
			#small_num = 0.0
			#data_sig_cov_i[np.isnan(data_sig_cov_i)] = 0.0
			###### if sigularity matrix then...
			#try:
			#	np.linalg.det(data_sig_cov_i+small_num)
			#except np.linalg.LinAlgError as err:
			#	if 'Singular matrix' in str(err):
			#		print(index)
			#		print('check problem np.linalg.det')
			#		print(data_sig_cov_i)
			#		print(np.linalg.det(data_sig_cov_i+small_num))
			#		print(np.log(np.linalg.det(data_sig_cov_i+small_num)))
			#		print('check problem ok')

			#if iteration_num == 8:
			#print(data_sig_cov_i)
			#if (np.isnan(np.log(np.linalg.det(data_sig_cov_i+small_num)+small_num))):
			#	print(index)
			#	print('check problem np.linalg.det')
			#	print(data_sig_cov_i)
			#	print(np.linalg.det(data_sig_cov_i+small_num))
			#	print(np.log(np.linalg.det(data_sig_cov_i+small_num)))
			#	print('check problem ok')

			#if (np.isnan(np.sum(np.linalg.inv(data_sig_cov_i+small_num)))):
			#try:
			#	np.linalg.inv(data_sig_cov_i+small_num)
			#except np.linalg.LinAlgError as err:
			#	if 'Singular matrix' in str(err):
			#		print(index)
			#		print('check problem np.linalg.inv')
			#		print(data_sig_cov_i)
			#		print(np.linalg.inv(data_sig_cov_i+small_num))
			#		print(np.log(np.linalg.inv(data_sig_cov_i+small_num)))
			#		print('check problem ok')

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

			###### calcualte qda score
			mvn_qda = -1/2*np.log(np.linalg.det(data_sig_cov_i)) - np.sum(np.dot((data_sig_matrix-data_sig_mean_i), np.linalg.inv(data_sig_cov_i)) * (data_sig_matrix-data_sig_mean_i), axis=1) + np.log(data_sig_p_i) #+ constant
			###### append qda score to mvn_qda_matrix
			mvn_qda_matrix.append(mvn_qda)

		##################
		###### convert mvn density score to p-value
		mvn_qda_matrix = np.transpose(np.array(mvn_qda_matrix))
		mvn_qda_p_matrix = np.exp(mvn_qda_matrix)#-mvn_qda_matrix.max(axis=1, keepdims=True))
		print(mvn_qda_p_matrix)
		###### calcualte posterior probability
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
		###### if 1-sum_mvn_p > 0.18.-> set to index-set 0
		#if iteration_num == 0:
		#	mvn_qda_p0_matrix = 1.0-np.max(mvn_qda_p_matrix, axis=1)
		#	mvn_qda_p_matrix = np.concatenate((mvn_qda_p_matrix, mvn_qda_p0_matrix.reshape(mvn_qda_p0_matrix.shape[0], 1)), axis=1)

		'''
		print('no pk:')
		print(np.sum(mvn_qda_p0_matrix>0.18.)
		'''
		##################
		###### set the index-set with highest p as the index-set label for each bin
		cluster_id = np.argmax(mvn_qda_p_matrix, axis=1)
		#print(cluster_id)
		#print(mvn_qda_p_matrix[0:3,:])
		#print(cluster_id[0:10])
		##################
		###### get the highest probability
		y_pred = np.amax(mvn_qda_p_matrix, axis=1)
		y_pred_mean = np.mean(mvn_qda_p_matrix, axis=1)
		#print(cluster_id.shape)
		#data_index_vec = np.append(data_index_vec, '_'.join(['0'] * ct_num))
		print(y_pred[0:20])
		print(y_pred_mean[0:20])
		print(stats.describe(y_pred))

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

		### replace rare index-sets by X_X_...
		for i in range(0,len(index_pred_vec)):
			index_pred = index_pred_vec[i]
			if data_index_pred_count[index_pred] < index_count_lim:
				print(index_pred)
				print(data_index_pred_count[index_pred])
				index_pred_vec[i] = '_'.join(['X'] * ct_num)

		### recount index-set number
		data_index_pred_count = {}
		for i in range(0,len(index_pred_vec)):
			index_pred = index_pred_vec[i]
			if not (index_pred in data_index_pred_count):
				data_index_pred_count[index_pred] = 1
			else:
				data_index_pred_count[index_pred] = data_index_pred_count[index_pred] +1


		###### get the number of peaks with different labels
		count_table = open(output_signal_mat+ '.iteration_' + str(iteration_num) + '_count_table.txt','w')
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

		###### total changed number
		change_num_0 = 0
		for i in range(0, len(index_pred_vec)):
			if data[i,1] != index_pred_vec[i]:
				change_num_0 = change_num_0+1

		change_num_all.append([change_num, change_num_0])
		print('change_num:')
		print(len(data_index_vec))
		print('iteration_num: ' + str(iteration_num) + '!!!' + str(change_num))
		print('iteration_num: ' + str(iteration_num) + '!!!' + str(change_num_0))
		
		### updata index in previous signal matrix
		data[:,1] = index_pred_vec

	###### write change_number vector
	change_num_all = np.array(change_num_all)
	write2d_array(change_num_all, output_signal_mat+'.change_num_all.txt')

	###### get signal matrix with updated index-set clusters
	data_sig_pred_vec = []
	###### get new cluster labels
	for i, info in zip(cluster_id, data):
		index = data_index_vec[i]
		data_sig_pred_vec.append(index)
	###### change the index-set labels in original matrix
	data[:,1] = data_sig_pred_vec
	###### write output
	write2d_array(data, output_signal_mat)

	###### get new index-set mean signal matrix
	data_sig_pred_dict = {}
	cluster_vec = []
	for records in data:
		cluster = records[1]
		if not (cluster in data_sig_pred_dict):
			data_sig_pred_dict[cluster] = [ records[2:] ]
			cluster_vec.append(cluster)
		else:
			data_sig_pred_dict[cluster].append(records[2:])

	###### sort index-set based on label
	cluster_vec = np.sort(np.array(cluster_vec))
	###### get mean signal matrix
	pred_cluster_mean = []
	for cluster in cluster_vec:
		sig_mat = np.array(data_sig_pred_dict[cluster], dtype=float)
		sig_mat_mean = np.mean(sig_mat, axis=0)
		pred_cluster_mean.append(sig_mat_mean)
	pred_cluster_mean = np.array(pred_cluster_mean)
	###### write output with index-set label
	data_sig_out = np.concatenate((cluster_vec.reshape((cluster_vec).shape[0],1), pred_cluster_mean), axis=1)
	write2d_array(data_sig_out, output_signal_mat+'.meansig.txt')


############################################################################

import getopt
import sys
def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:o:c:u:l:t:a:r:")
	except getopt.GetoptError:
		print 'time python oldlabel2newlabel.py -i signal_mat_file -l oldlabel2newlabel_mat_file -c old_label_col_in_signal_mat -o output_signal_mat'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'time python oldlabel2newlabel.py -i signal_mat_file -l oldlabel2newlabel_mat_file -c old_label_col_in_signal_mat -o output_signal_mat'		
		elif opt=="-i":
			input_signal_mat=str(arg.strip())				
		elif opt=="-o":
			output_signal_mat=str(arg.strip())
		elif opt=="-c":
			index_count_lim=int(arg.strip())
		elif opt=="-u":
			std_upper_lim=float(arg.strip())
		elif opt=="-l":
			std_lower_lim=float(arg.strip())
		elif opt=="-t":
			iteration_num_all=int(arg.strip())
		elif opt=="-a":
			alpha_for_empty_index_set=float(arg.strip())
		elif opt=="-r":
			random_seed=int(arg.strip())

	index_caller(input_signal_mat, output_signal_mat, index_count_lim, std_upper_lim, std_lower_lim, iteration_num_all, alpha_for_empty_index_set, random_seed)

if __name__=="__main__":
	main(sys.argv[1:])





