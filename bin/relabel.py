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


data = read2d_array('atac_20cell.sig.18.txt.ic.txt', 'str') 

ct_num = data.shape[1]-2

old_label_mat = []

for i in range(0, data.shape[0]):
	label_vec_tmp = data[i,1].split('_')
	if len(label_vec_tmp) == ct_num:
		old_label_mat.append(label_vec_tmp)
	else:
		old_label_mat.append(['N'+str(i)] * ct_num)

label_mat = np.array(old_label_mat)
label_vec_old = data[:,1]
label_vec_uniq_old = np.unique(label_vec_old, return_index = False)

#print(label_vec_uniq_old)

new_label_mat = []
print(label_mat.shape)
for i in range(0, ct_num):
	print('cell type: ' + str(i))
	label_vector = label_mat[:,i]
	signal_vector = (data[:,2+i]).astype(float)
	new_label = []

	for j in range(0,100):
		print('iteration: ' + str(j))
		dif_label_num = 0
		label_vector_pre = label_vector.copy()
		signal_null = signal_vector[label_vector=='0']
		###
		for index_set_old in label_vec_uniq_old:
			used_id_tmp = label_vec_old==index_set_old
			index_set_sig_vec = signal_vector[used_id_tmp]
			### get t-sample test 
			t_stat, tow_side_pval = stats.ttest_ind(signal_null, index_set_sig_vec, equal_var=False)
			one_side_pval = stats.t.cdf(t_stat, len(signal_null) + len(index_set_sig_vec) - 2) * 2
			if np.isnan(one_side_pval):
				print(signal_null[0:10])
				print(signal_null.shape)
				print(index_set_sig_vec[0:10])
				print(index_set_sig_vec.shape)
			#print(one_side_pval)
			if one_side_pval < 0.05:
				index_ct = '1'
			else:
				index_ct = '0'
			#print('check labels')
			#print(label_vector[used_id_tmp][0:10])
			label_vector[used_id_tmp] = index_ct
			#print(label_vector[used_id_tmp][0:10])
			#print(label_vector_pre[used_id_tmp][0:2])
			#print(label_vector[used_id_tmp][0:2])
			#print(label_vector_pre[used_id_tmp][0]!=label_vector[used_id_tmp][0])
			#print('check labels DONE')

			if label_vector_pre[used_id_tmp][0]!=label_vector[used_id_tmp][0]:
				dif_label_num = dif_label_num +1
		print(dif_label_num)
		if dif_label_num == 0:
			print('Converged!')
			break
			

	###### append label_vector vector to new_label_mat
	new_label_mat.append(label_vector)

new_label_mat = np.transpose(np.array(new_label_mat))
print(new_label_mat.shape)


new_label_vec = []
for i in range(0, data.shape[0]):
	new_label_vec.append('_'.join(new_label_mat[i,:]))

new_label_vec = np.array(new_label_vec)
print(new_label_vec[0:10])
data[:,1] = new_label_vec

write2d_array(data, 'atac_20cell.sig.18.txt.ic.Nrelabel.txt')

