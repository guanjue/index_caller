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
data = read2d_array('atac_20cell.sig.5.txt', str)



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
write2d_array(data_sig_out, 'qda_atac_wg.mvn_check.5.od.txt')




