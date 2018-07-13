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


data = read2d_array('atac_20cell.sig.18.txt.ic.txt', str) 

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
print('old_label number:')
print(label_vec_uniq_old.shape)
#print(label_vec_uniq_old)

new_label_mat = []
print(label_mat.shape)
for i in range(0, ct_num):
	print('cell type: ' + str(i))
	label_vector = label_mat[:,i]
	signal_vector = (data[:,2+i]).astype(float)
	new_label = []

	signal_0 = signal_vector[label_vector=='0']
	signal_1 = signal_vector[label_vector=='1']	
	print(np.mean(signal_0))
	print(np.var(signal_0))
	print(np.mean(signal_1))
	print(np.var(signal_1))
	for j in range(0,100):
		print('iteration: ' + str(j))
		dif_label_num = 0
		label_vector_pre = label_vector.copy()
		signal_0 = signal_vector[label_vector=='0']
		signal_1 = signal_vector[label_vector=='1']		
		###
		for index_set_old in label_vec_uniq_old:
			used_id_tmp = (label_vec_old==index_set_old)
			index_set_sig_vec = signal_vector[used_id_tmp]
			### get t-sample test 
			t_stat, tow_side_pval = stats.ttest_ind(signal_0, index_set_sig_vec, equal_var=False)
			#one_side_pval_0 = stats.t.cdf(t_stat, len(signal_0) + len(index_set_sig_vec) - 2) * 2
			if (tow_side_pval < (1e-2 / len(label_vec_uniq_old))) & (np.mean(signal_0) < np.mean(index_set_sig_vec)):
				index_ct = '1'
				#t_stat_1, tow_side_pval_1 = stats.ttest_ind(signal_1, index_set_sig_vec, equal_var=False)
				#one_side_pval_1 = stats.t.cdf(t_stat_1, len(signal_1) + len(index_set_sig_vec) - 2) * 2
				#if (tow_side_pval_1 < (1e-2 / len(label_vec_uniq_old))) & (np.mean(signal_1) < np.mean(index_set_sig_vec)):
				#	index_ct = '2'
				#else:
				#	index_ct = '1'
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
			print(index_ct)
			print(np.mean(signal_0))
			print(np.var(signal_0))
			print(np.mean(signal_1))
			print(np.var(signal_1))
			print(np.mean(index_set_sig_vec))
			break
				

	###### append label_vector vector to new_label_mat
	new_label_mat.append(label_vector)

new_label_mat = np.transpose(np.array(new_label_mat))
print(new_label_mat.shape)


new_label_vec = []
for i in range(0, data.shape[0]):
	new_label_vec.append('_'.join(new_label_mat[i,:]))

new_label_vec = np.array(new_label_vec, dtype="S50")
new_label_vec_uniq = np.unique(new_label_vec)
print('new_label number:')
print(new_label_vec_uniq.shape)
final_uniq_label = {}
l = 0
for old_label in label_vec_uniq_old:
	l = l+1
	print('label: ' + str(l))
	new_label = new_label_vec[label_vec_old==old_label][0]
	print(new_label)
	if not (new_label in final_uniq_label):
		final_uniq_label[new_label] = 0
	else:
		new_label = new_label + ':' + str(l)
		print(new_label)
		final_uniq_label[new_label] = 0
		new_label_vec[label_vec_old==old_label] = new_label
		print('after change label')
		print(new_label_vec[label_vec_old==old_label][0])



print(new_label_vec[0:10])
#data[:,1].astype('S50')
#data[:,1] = new_label_vec
new_label_vec = np.reshape(new_label_vec, (new_label_vec.shape[0],1))
print(new_label_vec.shape)
print(data.shape)
data_change_1 = np.concatenate((data, new_label_vec), axis=1)
#data_change_1[:,1] = data_change_1[:,data_change_1.shape[1]-1]
#data_change_1 = data_change_1[:,:-1]

write2d_array(data_change_1, 'atac_20cell.sig.18.txt.ic.Nrelabel1.txt')

